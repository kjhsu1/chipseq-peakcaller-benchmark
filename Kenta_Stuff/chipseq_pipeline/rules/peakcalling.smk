# peakcalling.smk — MACS2 peak calling using Bowtie2 alignments

def macs2_inputs(wc):
    r = find_row(wc.run_id)
    return {
        "treat": f"results/{wc.run_id}/bowtie2/treat/aligned.sorted.bam",
        "ctrl":  f"results/{wc.run_id}/bowtie2/con/aligned.sorted.bam",
    }

# MACS2
rule call_peaks_macs2:
    input: macs2_inputs
    output:
        narrow = "results/{run_id}/peaks/macs2/{run_id}_peaks.narrowPeak"
    params:
        gsize = lambda wc: macs2_gsize(find_row(wc.run_id)),
        flags = lambda wc: macs2_flags()
    shell:
        r"""
        macs2 callpeak \
          -t {input.treat} \
          -c {input.ctrl} \
          -g {params.gsize} \
          -n {wildcards.run_id} \
          {params.flags}
        """


# Collector returns only the peak files for the caller chosen in SAMPLES

def peaks_all():
    outs = []
    for r in SAMPLES:
        rid = r["run_id"]
        outs.append(f"results/{rid}/peaks/macs2/{rid}_peaks.narrowPeak")
    return outs

rule peaks_done:
    input: peaks_all()
