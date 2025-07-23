"""Simulate CHIP-seq reads and probability mass files"""

rule simulate_reads:
    input:
        lambda wildcards: f"data/{wildcards.genome}"
    output:
        fasta = "reads/{genome}/{id}/reads.fa",
        pmf   = "reads/{genome}/{id}/pmf.csv"
    params:
        coverage = lambda wc: wc.coverage,
        num_peaks = lambda wc: wc.num_peaks,
        peak_tallness = lambda wc: wc.peak_tallness,
        peak_broadness = lambda wc: wc.peak_broadness,
        fragment_length = lambda wc: wc.fragment_length
    log:
        "logs/simulate_reads/{genome}_{id}.log"
    conda:
        "envs/sim.yml"
    shell:
        """
        python chipsim.py \
            --genome {input} \
            --coverage {params.coverage} \
            --num-peaks {params.num_peaks} \
            --tallness {params.peak_tallness} \
            --broadness {params.peak_broadness} \
            --fragment-length {params.fragment_length} \
            --out-fasta {output.fasta} \
            --out-pmf {output.pmf} \
            > {log} 2>&1
        """
