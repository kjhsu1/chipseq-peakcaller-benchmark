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


def ingested_peak_inputs(wc):
    treat_roles = selected_roles(wc.study_id, "treatment_")
    ctrl_roles = selected_roles(wc.study_id, "control_")
    if not treat_roles:
        raise ValueError(f"No selected treatment roles found for {wc.study_id}")
    if not ctrl_roles:
        raise ValueError(f"No selected control roles found for {wc.study_id}")
    return {
        "treat": expand(
            "data/aligned/{study_id}/{role}/aligned.sorted.bam",
            study_id=[wc.study_id] * len(treat_roles),
            role=treat_roles,
        ),
        "ctrl": expand(
            "data/aligned/{study_id}/{role}/aligned.sorted.bam",
            study_id=[wc.study_id] * len(ctrl_roles),
            role=ctrl_roles,
        ),
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


rule call_peaks_macs2_realstudy:
    input: unpack(macs2_inputs)
    output:
        bed="results/{run_id}/peaks/macs2/{run_id}_peaks.bed"
    params:
        gsize=lambda wc: "1.0e8",
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
        gsize=lambda wc: "1.0e8",
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
