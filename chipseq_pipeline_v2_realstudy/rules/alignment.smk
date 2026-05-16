def align_all():
    outs = []
    for row in RUNS:
        alg = row["aligner"]
        run_id = row["run_id"]
        for cond in ("con", "treat"):
            base = f"results/{run_id}/{alg}/{cond}/aligned.sorted.bam"
            outs.extend([base, f"{base}.bai"])
    return outs


def ingest_fastq_align_all():
    outs = []
    for row in selected_fastq_rows():
        base = f"data/aligned/{row['study_id']}/{row['role']}/aligned.sorted.bam"
        outs.extend([base, f"{base}.bai"])
    return outs


def ingest_fastq_input(wc):
    row = find_manifest_row(wc.study_id, wc.role)
    return row["local_path"]


def ingest_fastq_index(wc):
    row = find_manifest_row(wc.study_id, wc.role)
    assembly = row.get("assembly", "unknown") or "unknown"
    return f"indexes/{assembly}/bowtie2_index"


rule align_bowtie2_realstudy:
    input:
        r1="results/{run_id}/{cond}/reads_R1.fasta",
        r2="results/{run_id}/{cond}/reads_R2.fasta",
    output:
        bam="results/{run_id}/bowtie2/{cond}/aligned.sorted.bam",
        bai="results/{run_id}/bowtie2/{cond}/aligned.sorted.bam.bai",
    threads: 4
    params:
        index=lambda wc: f"indexes/{find_row(wc.run_id).get('assembly', 'unknown')}/bowtie2_index"
    shell:
        r"""
        bowtie2 -p {threads} -f -x {params.index} -1 {input.r1} -2 {input.r2} \
          | samtools view -b - \
          | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """


rule align_downloaded_fastq_realstudy:
    input:
        fastq=ingest_fastq_input
    output:
        bam="data/aligned/{study_id}/{role}/aligned.sorted.bam",
        bai="data/aligned/{study_id}/{role}/aligned.sorted.bam.bai",
    threads: 4
    params:
        index=ingest_fastq_index
    shell:
        r"""
        bowtie2 -p {threads} -x {params.index} -U {input.fastq} \
          | samtools view -b - \
          | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """


rule align_done:
    input: align_all()


rule ingest_align_done:
    input: ingest_fastq_align_all()
