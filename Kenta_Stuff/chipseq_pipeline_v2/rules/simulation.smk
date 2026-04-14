# simulation.smk — emits FASTA R1/R2, pmf.csv, and planted peak centers per (run_id, cond)

EMIT_PMF_CSV = bool(config.get("emit_pmf_csv", True))
PMF_OUTPUT_PATTERN = (
    "results/{run_id}/{cond}/pmf.csv"
    if EMIT_PMF_CSV else
    "results/{run_id}/{cond}/pmf.disabled"
)

def sim_outputs_for(run_id, cond):
    base = f"results/{run_id}/{cond}"
    outputs = [
        f"{base}/reads_R1.fasta",
        f"{base}/reads_R2.fasta",
        f"{base}/planted_peaks.bed",
    ]
    outputs.append(f"{base}/pmf.csv" if EMIT_PMF_CSV else f"{base}/pmf.disabled")
    return outputs

rule simulate_reads:
    output:
        r1  = "results/{run_id}/{cond}/reads_R1.fasta",
        r2  = "results/{run_id}/{cond}/reads_R2.fasta",
        peaks = "results/{run_id}/{cond}/planted_peaks.bed",
        pmf = PMF_OUTPUT_PATTERN,
    params:
        cov       = lambda wc: (find_row(wc.run_id)["coverage_ctrl"]
                                if wc.cond == "con" else find_row(wc.run_id)["coverage_treat"]),
        tfcount   = lambda wc: (find_row(wc.run_id)["tf_peak_count_ctrl"]
                                if wc.cond == "con" else find_row(wc.run_id)["tf_peak_count_treat"]),
        tf_sigma  = lambda wc: find_row(wc.run_id)["tf_sigma"],
        tf_enrich = lambda wc: find_row(wc.run_id)["tf_enrich"],
        frag_len  = lambda wc: find_row(wc.run_id)["fragment_length"],
        read_len  = lambda wc: find_row(wc.run_id)["read_length"],
        nb_k      = lambda wc: find_row(wc.run_id)["nb_k"],
        seed      = lambda wc: find_row(wc.run_id)["seed"],
        fasta     = lambda wc: fasta_path(find_row(wc.run_id)),
        acc_bed   = lambda wc: acc_bed_path(find_row(wc.run_id)),
        gc_bias   = lambda wc: gc_bias_path(find_row(wc.run_id)),
        acc_weight = lambda wc: config.get("acc_weight", 1.0),
        tf_exp    = lambda wc: find_row(wc.run_id)["tf_exp"],
        gc_exp    = lambda wc: find_row(wc.run_id)["gc_exp"],
        acc_exp   = lambda wc: find_row(wc.run_id)["acc_exp"],
        map_coverage_pct = lambda wc: find_row(wc.run_id)["map_coverage_pct"],
        map_sigma = lambda wc: find_row(wc.run_id)["map_sigma"],
        map_enrich = lambda wc: find_row(wc.run_id)["map_enrich"],
        map_exp = lambda wc: find_row(wc.run_id)["map_exp"],
        pmf_arg = lambda wc, output: (
            f"--pmf_csv {output.pmf}" if EMIT_PMF_CSV else "--skip_pmf_csv"
        ),
        skip_pmf = lambda wc: "true" if not EMIT_PMF_CSV else "false",
    shell:
        r"""
        python -m scripts.updated_chip_seq \
          --fasta {params.fasta} \
          --coverage {params.cov} \
          --tf_peak_count {params.tfcount} \
          --fragment_length {params.frag_len} \
          --read_length {params.read_len} \
          --tf_sigma {params.tf_sigma} \
          --tf_enrich {params.tf_enrich} \
          --accessibility_bed {params.acc_bed} \
          --acc_weight {params.acc_weight} \
          --gc_bias_params {params.gc_bias} \
          --tf_exp {params.tf_exp} \
          --gc_exp {params.gc_exp} \
          --acc_exp {params.acc_exp} \
          --map_coverage_pct {params.map_coverage_pct} \
          --map_sigma {params.map_sigma} \
          --map_enrich {params.map_enrich} \
          --map_exp {params.map_exp} \
          --seed {params.seed} \
          --nb_k {params.nb_k} \
          --output_fasta1 {output.r1} \
          --output_fasta2 {output.r2} \
          {params.pmf_arg} \
          --planted_peaks_bed {output.peaks}
        if [ "{params.skip_pmf}" = "true" ]; then
          touch {output.pmf}
        fi
        """

def sim_all():
    outs = []
    for r in SAMPLES:
        outs += sim_outputs_for(r["run_id"], "con")
        outs += sim_outputs_for(r["run_id"], "treat")
    return outs

rule sim_done:
    input: sim_all()
