# alignment.smk — produces sorted & indexed BAMs per (run_id, cond) for the chosen aligner

# Bowtie2
rule align_bowtie2:
    input:
        r1 = "results/{run_id}/{cond}/reads_R1.fasta",
        r2 = "results/{run_id}/{cond}/reads_R2.fasta",
    output:
        bam = "results/{run_id}/bowtie2/{cond}/aligned.sorted.bam",
        bai = "results/{run_id}/bowtie2/{cond}/aligned.sorted.bam.bai",
    threads: 4
    params:
        bt2 = lambda wc: bowtie2_index(find_row(wc.run_id))
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        bowtie2 -f -x {params.bt2} -1 {input.r1} -2 {input.r2} \
          | samtools view -b - \
          | samtools sort -o {output.bam}
        samtools index {output.bam}
        """

# BWA-MEM
rule align_bwa_mem:
    input:
        r1 = "results/{run_id}/{cond}/reads_R1.fasta",
        r2 = "results/{run_id}/{cond}/reads_R2.fasta",
    output:
        bam = "results/{run_id}/bwa-mem/{cond}/aligned.sorted.bam",
        bai = "results/{run_id}/bwa-mem/{cond}/aligned.sorted.bam.bai",
    threads: 4
    params:
        bwa = lambda wc: bwa_index(find_row(wc.run_id))
    shell:
        r"""
        mkdir -p $(dirname {output.bam})
        bwa mem {params.bwa} {input.r1} {input.r2} \
          | samtools view -b - \
          | samtools sort -o {output.bam}
        samtools index {output.bam}
        """

# Collector returns only the BAMs for the aligner chosen in SAMPLES
def align_all():
    outs = []
    for r in SAMPLES:
        rid = r["run_id"]; alg = r["aligner"]
        for cond in ("con", "treat"):
            base = f"results/{rid}/{alg}/{cond}/aligned.sorted.bam"
            outs += [base, f"{base}.bai"]
    return outs

rule align_done:
    input: align_all()

