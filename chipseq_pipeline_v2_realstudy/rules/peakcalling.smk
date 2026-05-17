def macs2_inputs(wc):
    row = find_row(wc.run_id)
    return {
        "treat": f"results/{wc.run_id}/{row['aligner']}/treat/aligned.sorted.bam",
        "ctrl": f"results/{wc.run_id}/{row['aligner']}/con/aligned.sorted.bam",
    }


def peaks_all():
    outs = []
    for row in RUNS:
        run_id = row["run_id"]
        outs.append(f"results/{run_id}/peaks/macs2/{run_id}_peaks.bed")
    return outs


def selected_roles(study_id, prefix):
    rows = []
    for row in selected_manifest_rows():
        if row["study_id"] != study_id:
            continue
        if str(row.get("role", "")).startswith(prefix):
            rows.append(row["role"])
    return sorted(rows)


def merged_study_bam(study_id, group):
    return f"analysis_outputs/realstudy_bams/{study_id}/{group}.merged.bam"


def merged_study_bam_dependencies(wc):
    prefix = "treatment_" if wc.group == "treatment" else "control_"
    rows = selected_rows_for_prefix(wc.study_id, prefix)
    if not rows:
        raise ValueError(f"No selected rows found for {wc.study_id} / {wc.group}")
    deps = []
    for row in rows:
        if str(row.get("needs_alignment", "")).strip().lower() == "true":
            deps.append(f"data/aligned/{row['study_id']}/{row['role']}/aligned.sorted.bam")
        else:
            deps.append(f"analysis_outputs/realstudy_ingest_prep/download_markers/{row['study_id']}/{row['role']}.done")
    return deps


def merged_study_bam_paths(wc):
    prefix = "treatment_" if wc.group == "treatment" else "control_"
    rows = selected_rows_for_prefix(wc.study_id, prefix)
    if not rows:
        raise ValueError(f"No selected rows found for {wc.study_id} / {wc.group}")
    return [manifest_bam_path(row) for row in rows]


def ingested_peak_inputs(wc):
    return {
        "treat": merged_study_bam(wc.study_id, "treatment"),
        "ctrl": merged_study_bam(wc.study_id, "control"),
    }


def ingested_peak_targets():
    return [
        f"analysis_outputs/realstudy_peakcalls/{study_id}/{study_id}_peaks.bed"
        for study_id in selected_study_ids()
    ]


def ingested_peak_mode(wc):
    observed = observed_study(wc.study_id)
    mode = str(observed.get("macs2_mode", "")).strip()
    if mode:
        return mode
    for row in selected_manifest_rows():
        if row["study_id"] == wc.study_id:
            return "broad" if str(row.get("signal_class", "")).strip() == "broad" else "narrow"
    return "narrow"


rule merge_selected_realstudy_bams:
    input:
        deps=merged_study_bam_dependencies
    output:
        bam="analysis_outputs/realstudy_bams/{study_id}/{group}.merged.bam",
        bai="analysis_outputs/realstudy_bams/{study_id}/{group}.merged.bam.bai",
    threads: 4
    params:
        bams=merged_study_bam_paths
    shell:
        r"""
        mkdir -p analysis_outputs/realstudy_bams/{wildcards.study_id}
        samtools merge -f -@ {threads} {output.bam}.unsorted.bam {params.bams}
        samtools sort -@ {threads} -o {output.bam} {output.bam}.unsorted.bam
        rm -f {output.bam}.unsorted.bam
        samtools index {output.bam}
        """


rule call_peaks_macs2_realstudy:
    input: unpack(macs2_inputs)
    output:
        bed="results/{run_id}/peaks/macs2/{run_id}_peaks.bed"
    params:
        gsize=lambda wc: genome_size_for_assembly(find_row(wc.run_id).get("assembly", "unknown")),
        mode=lambda wc: find_row(wc.run_id).get("macs2_mode", "narrow"),
        mode_flag=lambda wc: "--broad" if find_row(wc.run_id).get("macs2_mode", "narrow") == "broad" else "",
        outdir=lambda wc: f"results/{wc.run_id}/peaks/macs2"
    shell:
        r"""
        macs2 callpeak \
          -t {input.treat} \
          -c {input.ctrl} \
          -g {params.gsize} \
          -n {wildcards.run_id} \
          {params.mode_flag} \
          --outdir {params.outdir}
        if [[ "{params.mode}" == "broad" ]]; then
          cp {params.outdir}/{wildcards.run_id}_peaks.broadPeak {output.bed}
        else
          cp {params.outdir}/{wildcards.run_id}_peaks.narrowPeak {output.bed}
        fi
        """


rule call_peaks_macs2_ingested_realstudy:
    input: unpack(ingested_peak_inputs)
    output:
        bed="analysis_outputs/realstudy_peakcalls/{study_id}/{study_id}_peaks.bed"
    params:
        gsize=lambda wc: genome_size_for_assembly(study_assembly(wc.study_id)),
        mode=ingested_peak_mode,
        mode_flag=lambda wc: "--broad" if ingested_peak_mode(wc) == "broad" else "",
        outdir=lambda wc: f"analysis_outputs/realstudy_peakcalls/{wc.study_id}/macs2"
    shell:
        r"""
        macs2 callpeak \
          -t {input.treat} \
          -c {input.ctrl} \
          -g {params.gsize} \
          -n {wildcards.study_id} \
          {params.mode_flag} \
          --outdir {params.outdir}
        if [[ "{params.mode}" == "broad" ]]; then
          cp {params.outdir}/{wildcards.study_id}_peaks.broadPeak {output.bed}
        else
          cp {params.outdir}/{wildcards.study_id}_peaks.narrowPeak {output.bed}
        fi
        """


rule peaks_done:
    input: peaks_all()


rule ingested_peaks_done:
    input: ingested_peak_targets()
