# alignment.smk — produces sorted & indexed BAMs using Bowtie2

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
        bowtie2 -f -x {params.bt2} -1 {input.r1} -2 {input.r2} \
          | samtools view -b - \
          | samtools sort -o {output.bam}
        samtools index {output.bam}
        """

def align_all():
    outs = []
    for r in SAMPLES:
        rid = r["run_id"]
        for cond in ("con", "treat"):
            base = f"results/{rid}/bowtie2/{cond}/aligned.sorted.bam"
            outs += [base, f"{base}.bai"]
    return outs

rule align_done:
    input: align_all()
