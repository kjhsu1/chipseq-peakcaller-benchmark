def v2_file_path(file_id):
    return V2_FILE_BY_ID[file_id]["expected_local_path"]


def v2_aligned_bam(file_id):
    return f"{v2_root()}/aligned/{file_id}/aligned.sorted.bam"


def v2_stage_marker(file_id):
    row = V2_FILE_BY_ID[file_id]
    return f"{v2_root()}/staging/{row['study_id']}/{file_id}.verified.tsv"


def v2_bowtie2_index_files():
    prefix = V2_CONFIG.get("bowtie2_index", "")
    if not prefix:
        return []
    return [f"{prefix}{suffix}" for suffix in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]]


def v2_parent_bam(study_id, role):
    return f"{v2_root()}/parents/{study_id}/{role}.parent.bam"


def v2_parent_qc(study_id, role):
    return f"{v2_root()}/qc/parents/{study_id}/{role}.csv"


def v2_study_files(study_id, role):
    return list(V2_CONFIG["studies"][study_id][f"{role}_files"])


def v2_parent_inputs(wc):
    return [v2_aligned_bam(file_id) for file_id in v2_study_files(wc.study_id, wc.role)]


def v2_parent_replicate_qc(wc):
    return [
        path
        for file_id in v2_study_files(wc.study_id, wc.role)
        for path in [
            f"{v2_root()}/qc/replicates/{file_id}.flagstat.txt",
            f"{v2_root()}/qc/replicates/{file_id}.stats.txt",
        ]
    ]


def v2_control_subset(study_id, seed, coverage):
    return f"{v2_root()}/control_subsamples/{study_id}/seed_{seed}/{coverage_label(coverage)}.bam"


def v2_peak_control(wc):
    row = V2_RUN_BY_ID[wc.run_id]
    if row["run_type"] == "full_control_anchor":
        return v2_parent_bam(row["study_id"], "control")
    return v2_control_subset(row["study_id"], row["seed"], row["control_coverage_x"])


def v2_peak_treatment(wc):
    row = V2_RUN_BY_ID[wc.run_id]
    return v2_parent_bam(row["study_id"], "treatment")


def v2_sampling_manifests():
    return [
        f"{v2_root()}/control_subsamples/{study_id}/seed_{seed}/sampling_manifest.csv"
        for study_id in V2_CONFIG.get("studies", {})
        for seed in V2_CONFIG.get("seeds", [])
    ]


rule build_realstudy_v2_design:
    input:
        config="configs/realstudy_v2.yaml",
        files=lambda wc: V2_CONFIG["manifests"]["files"]
    output:
        runs=f"{v2_root()}/design/run_table.csv",
        audit=f"{v2_root()}/design/design_audit.json"
    shell:
        r"""
        python scripts/build_realstudy_v2_design.py \
          --config {input.config} --files-manifest {input.files} \
          --run-table {output.runs} --audit-json {output.audit}
        """


rule stage_realstudy_v2_fastq:
    input:
        manifest=lambda wc: V2_CONFIG["manifests"]["files"]
    output:
        fastq="data/raw/realstudy_v2/{study_id}/{file_id}.fastq.gz",
        marker=f"{v2_root()}/staging/{{study_id}}/{{file_id}}.verified.tsv"
    params:
        local_flag=lambda wc: "--require-local" if config_bool("realstudy_require_local_inputs", True) else ""
    shell:
        r"""
        python scripts/stage_realstudy_v2_file.py \
          --manifest {input.manifest} --file-id {wildcards.file_id} \
          --marker {output.marker} {params.local_flag}
        """


rule align_realstudy_v2_fastq:
    input:
        fastq=lambda wc: v2_file_path(wc.file_id),
        marker=lambda wc: v2_stage_marker(wc.file_id),
        index=v2_bowtie2_index_files()
    output:
        bam=f"{v2_root()}/aligned/{{file_id}}/aligned.sorted.bam",
        bai=f"{v2_root()}/aligned/{{file_id}}/aligned.sorted.bam.bai",
        duplicates=f"{v2_root()}/qc/replicates/{{file_id}}.markdup.txt"
    threads: lambda wc: int(V2_CONFIG.get("threads", 8))
    params:
        index=lambda wc: V2_CONFIG["bowtie2_index"],
        sample=lambda wc: V2_FILE_BY_ID[wc.file_id]["sample_id"]
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {output.duplicates})
        bowtie2 -p {threads} -x {params.index} -U {input.fastq} \
          --rg-id {wildcards.file_id} --rg SM:{params.sample} \
          | samtools view -b - \
          | samtools sort -@ {threads} -o {output.bam}.unmarked.bam
        samtools markdup -s -f {output.duplicates} {output.bam}.unmarked.bam {output.bam}
        rm -f {output.bam}.unmarked.bam
        samtools index -@ {threads} {output.bam}
        """


rule qc_realstudy_v2_replicate:
    input:
        bam=f"{v2_root()}/aligned/{{file_id}}/aligned.sorted.bam",
        bai=f"{v2_root()}/aligned/{{file_id}}/aligned.sorted.bam.bai"
    output:
        flagstat=f"{v2_root()}/qc/replicates/{{file_id}}.flagstat.txt",
        stats=f"{v2_root()}/qc/replicates/{{file_id}}.stats.txt"
    shell:
        r"""
        mkdir -p $(dirname {output.flagstat})
        samtools flagstat {input.bam} > {output.flagstat}
        samtools stats {input.bam} > {output.stats}
        """


rule build_realstudy_v2_parent_bam:
    input:
        bams=v2_parent_inputs,
        qc=v2_parent_replicate_qc
    output:
        bam=f"{v2_root()}/parents/{{study_id}}/{{role}}.parent.bam",
        bai=f"{v2_root()}/parents/{{study_id}}/{{role}}.parent.bam.bai"
    threads: lambda wc: int(V2_CONFIG.get("threads", 8))
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        samtools merge -f -@ {threads} {output.bam} {input.bams}
        samtools index -@ {threads} {output.bam}
        """


rule summarize_realstudy_v2_parent_bam:
    input:
        bam=f"{v2_root()}/parents/{{study_id}}/{{role}}.parent.bam",
        bai=f"{v2_root()}/parents/{{study_id}}/{{role}}.parent.bam.bai"
    output:
        qc=f"{v2_root()}/qc/parents/{{study_id}}/{{role}}.csv",
        checksum=f"{v2_root()}/qc/parents/{{study_id}}/{{role}}.sha256"
    params:
        parent_id=lambda wc: f"{wc.study_id}__{wc.role}_parent",
        minimum=lambda wc: V2_CONFIG["minimum_eligible_control_fragments"] if wc.role == "control" else 0
    shell:
        r"""
        python scripts/summarize_realstudy_v2_parent_bam.py \
          --bam {input.bam} --parent-library-id {params.parent_id} \
          --study-id {wildcards.study_id} --role {wildcards.role} \
          --minimum-eligible-fragments {params.minimum} \
          --output-csv {output.qc} --checksum-output {output.checksum}
        """


rule build_realstudy_v2_control_series:
    input:
        bam=f"{v2_root()}/parents/{{study_id}}/control.parent.bam",
        bai=f"{v2_root()}/parents/{{study_id}}/control.parent.bam.bai",
        qc=f"{v2_root()}/qc/parents/{{study_id}}/control.csv",
        checksum=f"{v2_root()}/qc/parents/{{study_id}}/control.sha256"
    output:
        bam_0p5=f"{v2_root()}/control_subsamples/{{study_id}}/seed_{{seed}}/0p5x.bam",
        bai_0p5=f"{v2_root()}/control_subsamples/{{study_id}}/seed_{{seed}}/0p5x.bam.bai",
        bam_1=f"{v2_root()}/control_subsamples/{{study_id}}/seed_{{seed}}/1x.bam",
        bai_1=f"{v2_root()}/control_subsamples/{{study_id}}/seed_{{seed}}/1x.bam.bai",
        bam_2=f"{v2_root()}/control_subsamples/{{study_id}}/seed_{{seed}}/2x.bam",
        bai_2=f"{v2_root()}/control_subsamples/{{study_id}}/seed_{{seed}}/2x.bam.bai",
        bam_4=f"{v2_root()}/control_subsamples/{{study_id}}/seed_{{seed}}/4x.bam",
        bai_4=f"{v2_root()}/control_subsamples/{{study_id}}/seed_{{seed}}/4x.bam.bai",
        bam_8=f"{v2_root()}/control_subsamples/{{study_id}}/seed_{{seed}}/8x.bam",
        bai_8=f"{v2_root()}/control_subsamples/{{study_id}}/seed_{{seed}}/8x.bam.bai",
        bam_16=f"{v2_root()}/control_subsamples/{{study_id}}/seed_{{seed}}/16x.bam",
        bai_16=f"{v2_root()}/control_subsamples/{{study_id}}/seed_{{seed}}/16x.bam.bai",
        bam_32=f"{v2_root()}/control_subsamples/{{study_id}}/seed_{{seed}}/32x.bam",
        bai_32=f"{v2_root()}/control_subsamples/{{study_id}}/seed_{{seed}}/32x.bam.bai",
        ledger=f"{v2_root()}/control_subsamples/{{study_id}}/seed_{{seed}}/rank_ledger.sqlite",
        manifest=f"{v2_root()}/control_subsamples/{{study_id}}/seed_{{seed}}/sampling_manifest.csv"
    params:
        target_args=lambda wc: " ".join(
            f"--target {coverage_label(cov)}={next(row['requested_control_fragments'] for row in v2_sampled_runs() if row['study_id'] == wc.study_id and int(row['seed']) == int(wc.seed) and float(row['control_coverage_x']) == float(cov))}"
            for cov in V2_CONFIG["control_coverages_x"]
        ),
        output_args=lambda wc: " ".join(
            f"--output {coverage_label(cov)}={v2_control_subset(wc.study_id, wc.seed, cov)}"
            for cov in V2_CONFIG["control_coverages_x"]
        )
    shell:
        r"""
        python scripts/build_empirical_control_series.py \
          --parent-bam {input.bam} --source-checksum-file {input.checksum} \
          --study-id {wildcards.study_id} --control-parent-id {wildcards.study_id}__control_parent \
          --seed {wildcards.seed} --genome-size-bp {V2_CONFIG[genome_size_bp]} \
          --normalization-fragment-length-bp {V2_CONFIG[normalization_fragment_length_bp]} \
          --minimum-eligible-fragments {V2_CONFIG[minimum_eligible_control_fragments]} \
          {params.target_args} {params.output_args} \
          --ledger {output.ledger} --manifest {output.manifest} \
          --algorithm-version {V2_CONFIG[sampling_algorithm_version]}
        """


rule call_realstudy_v2_peaks:
    input:
        treatment=v2_peak_treatment,
        control=v2_peak_control,
        treatment_qc=lambda wc: v2_parent_qc(V2_RUN_BY_ID[wc.run_id]["study_id"], "treatment"),
        control_qc=lambda wc: v2_parent_qc(V2_RUN_BY_ID[wc.run_id]["study_id"], "control")
    output:
        peaks=f"{v2_root()}/peaks/{{run_id}}/peaks.bed",
        native=f"{v2_root()}/peaks/{{run_id}}/native_peaks.bed",
        status=f"{v2_root()}/peaks/{{run_id}}/status.tsv"
    params:
        mode=lambda wc: V2_RUN_BY_ID[wc.run_id]["macs2_mode"],
        mode_flag=lambda wc: "--broad" if V2_RUN_BY_ID[wc.run_id]["macs2_mode"] == "broad" else "",
        broad_cutoff=lambda wc: f"--broad-cutoff {V2_CONFIG['macs2']['broad_cutoff']}" if V2_RUN_BY_ID[wc.run_id]["macs2_mode"] == "broad" else "",
        native_name=lambda wc: f"{wc.run_id}_peaks.broadPeak" if V2_RUN_BY_ID[wc.run_id]["macs2_mode"] == "broad" else f"{wc.run_id}_peaks.narrowPeak",
        outdir=lambda wc: f"{v2_root()}/peaks/{wc.run_id}/macs2"
    shell:
        r"""
        mkdir -p {params.outdir}
        macs2 callpeak -t {input.treatment} -c {input.control} -f BAM \
          -g {V2_CONFIG[macs2][genome_size]} -n {wildcards.run_id} \
          --keep-dup {V2_CONFIG[macs2][keep_dup]} --qvalue {V2_CONFIG[macs2][qvalue]} \
          {params.mode_flag} {params.broad_cutoff} --outdir {params.outdir}
        cp {params.outdir}/{params.native_name} {output.native}
        cp {output.native} {output.peaks}
        printf 'status\tpeak_count\ncomplete\t%s\n' "$(wc -l < {output.peaks})" > {output.status}
        """
