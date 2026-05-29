# peakcalling.smk — pick the caller per run; inputs reference the chosen aligner’s BAMs

def macs2_inputs(wc):
    r = find_row(wc.run_id)
    return {
        "treat": result_path(wc.run_id, r["aligner"], "treat", "aligned.sorted.bam"),
        "ctrl":  result_path(wc.run_id, r["aligner"], "con", "aligned.sorted.bam"),
    }

def epic2_inputs(wc):
    r = find_row(wc.run_id)
    return {
        "treat": result_path(wc.run_id, r["aligner"], "treat", "aligned.sorted.bam"),
        "ctrl":  result_path(wc.run_id, r["aligner"], "con", "aligned.sorted.bam"),
    }

def homer_inputs(wc):
    r = find_row(wc.run_id)
    return {
        "treat": result_path(wc.run_id, r["aligner"], "treat", "aligned.sorted.bam"),
        "ctrl":  result_path(wc.run_id, r["aligner"], "con", "aligned.sorted.bam"),
    }

# MACS2
rule call_peaks_macs2:
    input: unpack(macs2_inputs)
    output:
        bed = result_path("{run_id}", "peaks", "macs2", "{run_id}_peaks.bed")
    params:
        gsize = lambda wc: macs2_gsize(find_row(wc.run_id)),
        flags = lambda wc: macs2_flags(),
        mode = lambda wc: find_row(wc.run_id).get("macs2_mode", "narrow"),
        mode_flag = lambda wc: ("--broad" if find_row(wc.run_id).get("macs2_mode", "narrow") == "broad" else ""),
        native_peak = lambda wc: (
            result_path(wc.run_id, "peaks", "macs2", f"{wc.run_id}_peaks.broadPeak")
            if find_row(wc.run_id).get("macs2_mode", "narrow") == "broad"
            else result_path(wc.run_id, "peaks", "macs2", f"{wc.run_id}_peaks.narrowPeak")
        ),
        ctrl_arg = lambda wc, input: (f"-c {input.ctrl}" if find_row(wc.run_id)["use_control"] else ""),
        outdir = lambda wc: result_path(wc.run_id, "peaks", "macs2")
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
        rm -f {params.outdir}/{wildcards.run_id}_peaks.xls
        rm -f {params.outdir}/{wildcards.run_id}_peaks.gappedPeak
        """

# EPIC2
rule call_peaks_epic2:
    # conda: "../envs/epic2.yml"
    input: unpack(epic2_inputs)
    output:
        bed = result_path("{run_id}", "peaks", "epic2", "{run_id}_domains.bed")
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

# HOMER
rule make_homer_tagdirs:
    input: unpack(homer_inputs)
    output:
        treat_done = result_path("{run_id}", "peaks", "homer", "{run_id}_treat_tagdir", "tagInfo.txt"),
        ctrl_done = result_path("{run_id}", "peaks", "homer", "{run_id}_ctrl_tagdir", "tagInfo.txt"),
    params:
        treat_dir = lambda wc: result_path(wc.run_id, "peaks", "homer", f"{wc.run_id}_treat_tagdir"),
        ctrl_dir = lambda wc: result_path(wc.run_id, "peaks", "homer", f"{wc.run_id}_ctrl_tagdir"),
    shell:
        r"""
        makeTagDirectory {params.treat_dir} {input.treat} -format sam
        makeTagDirectory {params.ctrl_dir} {input.ctrl} -format sam
        """

rule call_peaks_homer:
    input:
        treat_done = result_path("{run_id}", "peaks", "homer", "{run_id}_treat_tagdir", "tagInfo.txt"),
        ctrl_done = result_path("{run_id}", "peaks", "homer", "{run_id}_ctrl_tagdir", "tagInfo.txt"),
    output:
        txt = result_path("{run_id}", "peaks", "homer", "{run_id}_peaks.txt"),
        bed = result_path("{run_id}", "peaks", "homer", "{run_id}_peaks.bed"),
    params:
        treat_dir = lambda wc: result_path(wc.run_id, "peaks", "homer", f"{wc.run_id}_treat_tagdir"),
        ctrl_dir = lambda wc: result_path(wc.run_id, "peaks", "homer", f"{wc.run_id}_ctrl_tagdir"),
        style = lambda wc: homer_flags(find_row(wc.run_id))["style"],
        flags = lambda wc: homer_flags(find_row(wc.run_id))["flags"],
        ctrl_arg = lambda wc: (f"-i {result_path(wc.run_id, 'peaks', 'homer', f'{wc.run_id}_ctrl_tagdir')}" if find_row(wc.run_id)["use_control"] else ""),
    shell:
        r"""
        findPeaks {params.treat_dir} -style {params.style} {params.ctrl_arg} {params.flags} -o {output.txt}
        pos2bed.pl {output.txt} > {output.bed}
        """

# Collector returns only the peak files for the caller chosen in SAMPLES
def peaks_all():
    outs = []
    for r in SAMPLES:
        rid = r["run_id"]
        if r["peakcaller"] == "macs2":
            outs.append(result_path(rid, "peaks", "macs2", f"{rid}_peaks.bed"))
        elif r["peakcaller"] == "epic2":
            outs.append(result_path(rid, "peaks", "epic2", f"{rid}_domains.bed"))
        else:  # homer
            outs.append(result_path(rid, "peaks", "homer", f"{rid}_peaks.bed"))
    return outs

rule peaks_done:
    input: peaks_all()
