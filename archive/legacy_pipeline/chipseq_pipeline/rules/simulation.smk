# simulation.smk — emits FASTA R1/R2 and pmf.csv per (run_id, cond)

def sim_outputs_for(run_id, cond):
    base = f"results/{run_id}/{cond}"
    return [
        f"{base}/reads_R1.fasta",
        f"{base}/reads_R2.fasta",
        f"{base}/pmf.csv",
    ]

rule simulate_reads:
    output:
        r1  = "results/{run_id}/{cond}/reads_R1.fasta",
        r2  = "results/{run_id}/{cond}/reads_R2.fasta",
        pmf = "results/{run_id}/{cond}/pmf.csv",
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
        fasta     = lambda wc: fasta_path(find_row(wc.run_id)),
        acc_bed   = lambda wc: acc_bed_path(find_row(wc.run_id)),
        gc_bias   = lambda wc: gc_bias_path(find_row(wc.run_id)),
        tf_exp    = lambda wc: find_row(wc.run_id)["tf_exp"],
        gc_exp    = lambda wc: find_row(wc.run_id)["gc_exp"],
        acc_exp   = lambda wc: find_row(wc.run_id)["acc_exp"],
    shell:
        r"""
        python scripts/updated_chip_seq.py \
          --fasta {params.fasta} \
          --coverage {params.cov} \
          --tf_peak_count {params.tfcount} \
          --fragment_length {params.frag_len} \
          --read_length {params.read_len} \
          --tf_sigma {params.tf_sigma} \
          --tf_enrich {params.tf_enrich} \
          --accessibility_bed {params.acc_bed} \
          --gc_bias_params {params.gc_bias} \
          --tf_exp {params.tf_exp} \
          --gc_exp {params.gc_exp} \
          --acc_exp {params.acc_exp} \
          --nb_k {params.nb_k} \
          --output_fasta1 {output.r1} \
          --output_fasta2 {output.r2} \
          --pmf_csv {output.pmf}
        """

def sim_all():
    outs = []
    for r in SAMPLES:
        outs += sim_outputs_for(r["run_id"], "con")
        outs += sim_outputs_for(r["run_id"], "treat")
    return outs

rule sim_done:
    input: sim_all()

