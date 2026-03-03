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
        bed = "results/{run_id}/peaks/macs2/{run_id}_peaks.bed"
    params:
        gsize = lambda wc: macs2_gsize(find_row(wc.run_id)),
        flags = lambda wc: macs2_flags(),
        mode = lambda wc: find_row(wc.run_id).get("macs2_mode", "narrow"),
        mode_flag = lambda wc: ("--broad" if find_row(wc.run_id).get("macs2_mode", "narrow") == "broad" else ""),
        native_peak = lambda wc: (
            f"results/{wc.run_id}/peaks/macs2/{wc.run_id}_peaks.broadPeak"
            if find_row(wc.run_id).get("macs2_mode", "narrow") == "broad"
            else f"results/{wc.run_id}/peaks/macs2/{wc.run_id}_peaks.narrowPeak"
        ),
        ctrl_arg = lambda wc, input: (f"-c {input.ctrl}" if find_row(wc.run_id)["use_control"] else ""),
        outdir = lambda wc: f"results/{wc.run_id}/peaks/macs2"
    shell:
        r"""
        if [[ "{params.mode}" != "narrow" && "{params.mode}" != "broad" ]]; then
          echo "Unsupported macs2_mode: {params.mode}" >&2
          exit 1
        fi
        macs2 callpeak \
          -t {input.treat} \
          {params.ctrl_arg} \
          -g {params.gsize} \
          -n {wildcards.run_id} \
          {params.mode_flag} \
          {params.flags} \
          --outdir {params.outdir}
        cp {params.native_peak} {output.bed}
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
        ctrl_arg   = lambda wc, input: (f"--control {input.ctrl}" if find_row(wc.run_id)["use_control"] else ""),
        flags      = lambda wc: epic2_flags()
    shell:
        r"""
        epic2 \
          --treatment {input.treat} \
          {params.ctrl_arg} \
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
            outs.append(f"results/{rid}/peaks/macs2/{rid}_peaks.bed")
        else:  # epic2
            outs.append(f"results/{rid}/peaks/epic2/{rid}_domains.bed")
    return outs

rule peaks_done:
    input: peaks_all()
