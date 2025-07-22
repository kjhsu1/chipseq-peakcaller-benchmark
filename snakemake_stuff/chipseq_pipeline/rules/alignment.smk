"""Align simulated reads with Bowtie2 or BWA-MEM"""

rule align_reads:
    input:
        fa = "reads/{genome}/{id}/reads.fa"
    output:
        bam = "align/{genome}/{aligner}/{id}/aligned.bam",
        bai = "align/{genome}/{aligner}/{id}/aligned.bam.bai"
    params:
        idx = lambda wc: config['indexes'][wc.genome][wc.aligner]
    log:
        "logs/align_reads/{genome}_{aligner}_{id}.log"
    conda:
        "envs/align.yml"
    shell:
        """
        if [ "{wildcards.aligner}" = "bowtie2" ]; then
            bowtie2 -x {params.idx} -f {input.fa} | samtools sort -o {output.bam} -
        else
            bwa mem {params.idx} {input.fa} | samtools sort -o {output.bam} -
        fi
        samtools index {output.bam}
        """ > {log} 2>&1
