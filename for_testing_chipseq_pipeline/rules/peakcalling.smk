# peakcalling.smk — pick the caller per run; inputs reference the chosen aligner’s BAMs

def macs2_inputs(wc):
    r = find_row(wc.run_id)
    return {
        "treat": f"results/{wc.run_id}/{r['aligner']}/treat/aligned.sorted.bam",
        "ctrl":  f"results/{wc.run_id}/{r['aligner']}/con/aligned.sorted.bam",
    }

def epic2_inputs(wc):
    r = find_row(wc.run_id)
    return {
        "treat": f"results/{wc.run_id}/{r['aligner']}/treat/aligned.sorted.bam",
        "ctrl":  f"results/{wc.run_id}/{r['aligner']}/con/aligned.sorted.bam",
    }

# MACS2
rule call_peaks_macs2:
    input: unpack(macs2_inputs)
    output:
        narrow = "results/{run_id}/peaks/macs2/{run_id}_peaks.narrowPeak"
    params:
        gsize = lambda wc: macs2_gsize(find_row(wc.run_id)),
        flags = lambda wc: macs2_flags(),
        outdir = lambda wc: f"results/{wc.run_id}/peaks/macs2"
    shell:
        r"""
        macs2 callpeak \
          -t {input.treat} \
          -c {input.ctrl} \
          -g {params.gsize} \
          -n {wildcards.run_id} \
          {params.flags} \
          --outdir {params.outdir}
        """

# EPIC2
rule call_peaks_epic2:
    # conda: "../envs/epic2.yml"
    input: unpack(epic2_inputs)
    output:
        bed = "results/{run_id}/peaks/epic2/{run_id}_domains.bed"
    params:
        chromsizes = lambda wc: config["peakcallers"]["epic2"]["chromsizes"][find_row(wc.run_id)["genome"]],
        egf        = lambda wc: config["peakcallers"]["epic2"]["effective_genome_fraction"][find_row(wc.run_id)["genome"]],
        flags      = lambda wc: epic2_flags()
    shell:
        r"""
        epic2 \
          --treatment {input.treat} \
          --control  {input.ctrl} \
          --chromsizes {params.chromsizes} \
          --effective-genome-fraction {params.egf} \
          --output {output.bed} \
          {params.flags}
        """

# Collector returns only the peak files for the caller chosen in SAMPLES
def peaks_all():
    outs = []
    for r in SAMPLES:
        rid = r["run_id"]
        if r["peakcaller"] == "macs2":
            outs.append(f"results/{rid}/peaks/macs2/{rid}_peaks.narrowPeak")
        else:  # epic2
            outs.append(f"results/{rid}/peaks/epic2/{rid}_domains.bed")
    return outs

rule peaks_done:
    input: peaks_all()
