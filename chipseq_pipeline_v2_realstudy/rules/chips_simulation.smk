def chips_cfg(name, default=None):
    return config.get("chips", {}).get(name, default)


def chips_run_targets(cond):
    outs = []
    for row in chips_runs():
        run_id = row["run_id"]
        outs.extend([
            f"results_chips/{run_id}/{cond}/reads_R1.fastq",
            f"results_chips/{run_id}/{cond}/reads_R2.fastq",
        ])
    return outs


def chips_alignment_targets():
    outs = []
    for row in chips_runs():
        run_id = row["run_id"]
        for cond in ("con", "treat"):
            base = f"results_chips/{run_id}/bowtie2/{cond}/aligned.sorted.bam"
            outs.extend([base, f"{base}.bai"])
    return outs


def chips_peak_targets():
    return [
        f"results_chips/{row['run_id']}/peaks/macs2/{row['run_id']}_peaks.bed"
        for row in chips_runs()
    ]


def chips_selected_study_ids():
    configured = chips_cfg("study_ids")
    if configured:
        return list(configured)
    return selected_study_ids()


def chips_runs():
    selected = set(chips_selected_study_ids())
    return [row for row in RUNS if row.get("study_id") in selected]


def chips_model_targets():
    return [f"analysis_outputs/chips_models/{study_id}/{study_id}.json" for study_id in chips_selected_study_ids()]


def chips_treatment_bam(wc):
    rows = selected_rows_for_prefix(wc.study_id, "treatment_")
    if not rows:
        raise ValueError(f"No treatment role available for ChIPs learn: {wc.study_id}")
    return merged_study_bam(wc.study_id, "treatment")


def chips_model_for_run(wc):
    return f"analysis_outputs/chips_models/{find_row(wc.run_id)['study_id']}/{find_row(wc.run_id)['study_id']}.json"


def chips_reference_for_run(wc):
    assembly = find_row(wc.run_id).get("assembly", "unknown") or "unknown"
    return reference_fasta_for_assembly(assembly)


def chips_numreads(row, cond):
    assembly = row.get("assembly", "unknown") or "unknown"
    genome_size = float(genome_size_for_assembly(assembly))
    coverage = float(row["coverage_treat"] if cond == "treat" else row["coverage_ctrl"])
    fragment_length = float(row.get("fragment_length", 150))
    return max(1, int(round((coverage * genome_size) / fragment_length)))


def chips_numcopies(row):
    mode = str(row.get("macs2_mode", "narrow"))
    if mode == "broad":
        return int(chips_cfg("numcopies_histone", 25))
    return int(chips_cfg("numcopies_tf", 1000))


def chips_seed(row, cond):
    offset = 101 if cond == "treat" else 202
    return int(row["seed"]) * 1009 + offset


def chips_bowtie2_index(wc):
    assembly = find_row(wc.run_id).get("assembly", "unknown") or "unknown"
    return bowtie2_index_for_assembly(assembly)


rule chips_learn_model:
    input:
        bam=chips_treatment_bam,
        peaks="analysis_outputs/realstudy_peakcalls/{study_id}/{study_id}_peaks.bed",
    output:
        model="analysis_outputs/chips_models/{study_id}/{study_id}.json",
    params:
        binary=lambda wc: chips_cfg("binary", "chips"),
        score_col=lambda wc: chips_cfg("peak_score_column", 4),
        outprefix=lambda wc: f"analysis_outputs/chips_models/{wc.study_id}/{wc.study_id}",
    shell:
        r"""
        {params.binary} learn \
          -b {input.bam} \
          -p {input.peaks} \
          -t bed \
          -c {params.score_col} \
          -o {params.outprefix}
        """


rule chips_simreads_treat:
    input:
        model=chips_model_for_run,
        ref=chips_reference_for_run,
        peaks=lambda wc: f"analysis_outputs/realstudy_peakcalls/{find_row(wc.run_id)['study_id']}/{find_row(wc.run_id)['study_id']}_peaks.bed",
    output:
        r1="results_chips/{run_id}/treat/reads_R1.fastq",
        r2="results_chips/{run_id}/treat/reads_R2.fastq",
    params:
        binary=lambda wc: chips_cfg("binary", "chips"),
        score_col=lambda wc: chips_cfg("peak_score_column", 4),
        outprefix=lambda wc: f"results_chips/{wc.run_id}/treat/chips_treat",
        numreads=lambda wc: chips_numreads(find_row(wc.run_id), "treat"),
        numcopies=lambda wc: chips_numcopies(find_row(wc.run_id)),
        readlen=lambda wc: int(find_row(wc.run_id).get("read_length", 38)),
        seed=lambda wc: chips_seed(find_row(wc.run_id), "treat"),
        threads=lambda wc: int(chips_cfg("threads", 4)),
    shell:
        r"""
        mkdir -p results_chips/{wildcards.run_id}/treat
        {params.binary} simreads \
          -p {input.peaks} \
          -t bed \
          -c {params.score_col} \
          -f {input.ref} \
          -o {params.outprefix} \
          --model {input.model} \
          --numcopies {params.numcopies} \
          --numreads {params.numreads} \
          --readlen {params.readlen} \
          --paired \
          --seed {params.seed} \
          --thread {params.threads}
        mv {params.outprefix}_1.fastq {output.r1}
        mv {params.outprefix}_2.fastq {output.r2}
        """


rule chips_simreads_control:
    input:
        ref=chips_reference_for_run,
    output:
        r1="results_chips/{run_id}/con/reads_R1.fastq",
        r2="results_chips/{run_id}/con/reads_R2.fastq",
    params:
        binary=lambda wc: chips_cfg("binary", "chips"),
        outprefix=lambda wc: f"results_chips/{wc.run_id}/con/chips_control",
        numreads=lambda wc: chips_numreads(find_row(wc.run_id), "con"),
        numcopies=lambda wc: chips_numcopies(find_row(wc.run_id)),
        readlen=lambda wc: int(find_row(wc.run_id).get("read_length", 38)),
        seed=lambda wc: chips_seed(find_row(wc.run_id), "con"),
        threads=lambda wc: int(chips_cfg("threads", 4)),
    shell:
        r"""
        mkdir -p results_chips/{wildcards.run_id}/con
        {params.binary} simreads \
          -t wce \
          -f {input.ref} \
          -o {params.outprefix} \
          --numcopies {params.numcopies} \
          --numreads {params.numreads} \
          --readlen {params.readlen} \
          --paired \
          --seed {params.seed} \
          --thread {params.threads}
        mv {params.outprefix}_1.fastq {output.r1}
        mv {params.outprefix}_2.fastq {output.r2}
        """


rule align_bowtie2_chips:
    input:
        r1="results_chips/{run_id}/{cond}/reads_R1.fastq",
        r2="results_chips/{run_id}/{cond}/reads_R2.fastq",
    output:
        bam="results_chips/{run_id}/bowtie2/{cond}/aligned.sorted.bam",
        bai="results_chips/{run_id}/bowtie2/{cond}/aligned.sorted.bam.bai",
    threads: 4
    params:
        index=chips_bowtie2_index,
    shell:
        r"""
        bowtie2 -p {threads} -x {params.index} -1 {input.r1} -2 {input.r2} \
          | samtools view -b - \
          | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """


rule call_peaks_macs2_chips:
    input:
        treat="results_chips/{run_id}/bowtie2/treat/aligned.sorted.bam",
        ctrl="results_chips/{run_id}/bowtie2/con/aligned.sorted.bam",
    output:
        bed="results_chips/{run_id}/peaks/macs2/{run_id}_peaks.bed",
    params:
        gsize=lambda wc: genome_size_for_assembly(find_row(wc.run_id).get("assembly", "unknown")),
        mode=lambda wc: find_row(wc.run_id).get("macs2_mode", "narrow"),
        mode_flag=lambda wc: "--broad" if find_row(wc.run_id).get("macs2_mode", "narrow") == "broad" else "",
        outdir=lambda wc: f"results_chips/{wc.run_id}/peaks/macs2",
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


rule chips_done:
    input:
        chips_model_targets(),
        chips_run_targets("treat"),
        chips_run_targets("con"),
        chips_alignment_targets(),
        chips_peak_targets()
