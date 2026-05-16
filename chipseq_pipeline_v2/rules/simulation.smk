# simulation.smk — emits FASTA R1/R2, optional PMF artifact, and planted peak centers

EMIT_PMF_CSV = config_bool("emit_pmf_csv", True)
EMIT_PMF_DISABLED_MARKER = config_bool("emit_pmf_disabled_marker", True)


def pmf_output_path(cond):
    base = f"results/{{run_id}}/{cond}"
    if EMIT_PMF_CSV:
        return f"{base}/pmf.csv"
    if EMIT_PMF_DISABLED_MARKER:
        return f"{base}/pmf.disabled"
    return f"{base}/pmf.omitted"


def pmf_output(cond):
    path = pmf_output_path(cond)
    if EMIT_PMF_CSV or EMIT_PMF_DISABLED_MARKER:
        return path
    return temp(path)


def retained_sim_outputs():
    outputs = []
    for row in SAMPLES:
        run_id = row["run_id"]
        outputs.append(f"results/{run_id}/treat/planted_peaks.bed")
        outputs.append(f"results/{run_id}/treat/planted_peak_intervals.bed")
        if EMIT_PMF_CSV or EMIT_PMF_DISABLED_MARKER:
            outputs.append(pmf_output_path("con").format(run_id=run_id))
            outputs.append(pmf_output_path("treat").format(run_id=run_id))
    return outputs


rule simulate_reads_con:
    output:
        r1=temp("results/{run_id}/con/reads_R1.fasta"),
        r2=temp("results/{run_id}/con/reads_R2.fasta"),
        peaks=temp("results/{run_id}/con/planted_peaks.bed"),
        peak_intervals=temp("results/{run_id}/con/planted_peak_intervals.bed"),
        pmf=pmf_output("con"),
    params:
        cov=lambda wc: find_row(wc.run_id)["coverage_ctrl"],
        tfcount=lambda wc: find_row(wc.run_id)["tf_peak_count_ctrl"],
        tf_sigma=lambda wc: find_row(wc.run_id)["tf_sigma"],
        tf_enrich=lambda wc: find_row(wc.run_id)["tf_enrich"],
        frag_len=lambda wc: find_row(wc.run_id)["fragment_length"],
        read_len=lambda wc: find_row(wc.run_id)["read_length"],
        nb_k=lambda wc: find_row(wc.run_id)["nb_k"],
        seed=lambda wc: find_row(wc.run_id)["seed"],
        tf_seed=lambda wc: find_row(wc.run_id)["treat_tf_seed"],
        map_seed=lambda wc: find_row(wc.run_id)["control_map_seed"],
        read_seed=lambda wc: find_row(wc.run_id)["control_read_seed"],
        fasta=lambda wc: fasta_path(find_row(wc.run_id)),
        acc_bed=lambda wc: acc_bed_path(find_row(wc.run_id)),
        gc_bias=lambda wc: gc_bias_path(find_row(wc.run_id)),
        acc_weight=lambda wc: config.get("acc_weight", 1.0),
        tf_exp=lambda wc: find_row(wc.run_id)["tf_exp"],
        gc_exp=lambda wc: find_row(wc.run_id)["gc_exp"],
        acc_exp=lambda wc: find_row(wc.run_id)["acc_exp"],
        map_coverage_pct=lambda wc: find_row(wc.run_id)["map_coverage_pct"],
        map_sigma=lambda wc: find_row(wc.run_id)["map_sigma"],
        map_enrich=lambda wc: find_row(wc.run_id)["map_enrich"],
        map_exp=lambda wc: find_row(wc.run_id)["map_exp"],
        pmf_arg=lambda wc, output: (
            f"--pmf_csv {output.pmf}" if EMIT_PMF_CSV else "--skip_pmf_csv"
        ),
        skip_pmf=lambda wc: "true" if not EMIT_PMF_CSV else "false",
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
          --tf_seed {params.tf_seed} \
          --map_seed {params.map_seed} \
          --read_seed {params.read_seed} \
          --nb_k {params.nb_k} \
          --output_fasta1 {output.r1} \
          --output_fasta2 {output.r2} \
          {params.pmf_arg} \
          --planted_peaks_bed {output.peaks} \
          --planted_peak_intervals_bed {output.peak_intervals}
        if [ "{params.skip_pmf}" = "true" ]; then
          touch {output.pmf}
        fi
        """


rule simulate_reads_treat:
    output:
        r1=temp("results/{run_id}/treat/reads_R1.fasta"),
        r2=temp("results/{run_id}/treat/reads_R2.fasta"),
        peaks="results/{run_id}/treat/planted_peaks.bed",
        peak_intervals="results/{run_id}/treat/planted_peak_intervals.bed",
        pmf=pmf_output("treat"),
    params:
        cov=lambda wc: find_row(wc.run_id)["coverage_treat"],
        tfcount=lambda wc: find_row(wc.run_id)["tf_peak_count_treat"],
        tf_sigma=lambda wc: find_row(wc.run_id)["tf_sigma"],
        tf_enrich=lambda wc: find_row(wc.run_id)["tf_enrich"],
        frag_len=lambda wc: find_row(wc.run_id)["fragment_length"],
        read_len=lambda wc: find_row(wc.run_id)["read_length"],
        nb_k=lambda wc: find_row(wc.run_id)["nb_k"],
        seed=lambda wc: find_row(wc.run_id)["seed"],
        tf_seed=lambda wc: find_row(wc.run_id)["treat_tf_seed"],
        map_seed=lambda wc: find_row(wc.run_id)["treat_map_seed"],
        read_seed=lambda wc: find_row(wc.run_id)["treat_read_seed"],
        fasta=lambda wc: fasta_path(find_row(wc.run_id)),
        acc_bed=lambda wc: acc_bed_path(find_row(wc.run_id)),
        gc_bias=lambda wc: gc_bias_path(find_row(wc.run_id)),
        acc_weight=lambda wc: config.get("acc_weight", 1.0),
        tf_exp=lambda wc: find_row(wc.run_id)["tf_exp"],
        gc_exp=lambda wc: find_row(wc.run_id)["gc_exp"],
        acc_exp=lambda wc: find_row(wc.run_id)["acc_exp"],
        map_coverage_pct=lambda wc: find_row(wc.run_id)["map_coverage_pct"],
        map_sigma=lambda wc: find_row(wc.run_id)["map_sigma"],
        map_enrich=lambda wc: find_row(wc.run_id)["map_enrich"],
        map_exp=lambda wc: find_row(wc.run_id)["map_exp"],
        pmf_arg=lambda wc, output: (
            f"--pmf_csv {output.pmf}" if EMIT_PMF_CSV else "--skip_pmf_csv"
        ),
        skip_pmf=lambda wc: "true" if not EMIT_PMF_CSV else "false",
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
          --tf_seed {params.tf_seed} \
          --map_seed {params.map_seed} \
          --read_seed {params.read_seed} \
          --nb_k {params.nb_k} \
          --output_fasta1 {output.r1} \
          --output_fasta2 {output.r2} \
          {params.pmf_arg} \
          --planted_peaks_bed {output.peaks} \
          --planted_peak_intervals_bed {output.peak_intervals}
        if [ "{params.skip_pmf}" = "true" ]; then
          touch {output.pmf}
        fi
        """


rule sim_done:
    input: retained_sim_outputs()
