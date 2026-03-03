# Master Snakefile

configfile: "config.yaml"

from itertools import product

def get_peakcaller_list(cfg):
    if "peakcaller_list" in cfg:
        return cfg["peakcaller_list"]
    if isinstance(cfg.get("peakcallers"), list):
        return cfg["peakcallers"]
    raise ValueError(
        "Config must define 'peakcaller_list' (preferred) or legacy list-valued "
        "'peakcallers'. The 'peakcallers' parameter dictionary is not a sweep list."
    )

# ---------- sweep catalog ----------
def build_samples(cfg):
    S = []
    rid = 0
    peakcaller_list = get_peakcaller_list(cfg)
    use_control_values = cfg.get("use_control", [True])
    for (
        genome, acc_key, gc_key, frag, read, nbk, aligner, peakcaller,
        tf_exp, gc_exp, acc_exp, map_coverage_pct, map_sigma, map_enrich,
        map_exp, use_control
    ) in product(
        cfg["genomes"], cfg["acc_beds"], cfg["gc_bias_sets"],
        cfg["fragment_length"], cfg["read_length"], cfg["nb_k"],
        cfg["aligners"], peakcaller_list,
        cfg["tf_exp"], cfg["gc_exp"], cfg["acc_exp"],
        cfg.get("map_coverage_pct", [0.0]),
        cfg.get("map_sigma", [5.0]),
        cfg.get("map_enrich", [1.0]),
        cfg.get("map_exp", [1.0]),
        use_control_values
    ):
        for cov_t, cov_c in product(cfg["coverage_treat"], cfg["coverage_ctrl"]):
            for tpc, tsig, tenr in product(cfg["tf_peak_count_treat"],
                                           cfg["tf_sigma"], cfg["tf_enrich"]):
                rid += 1
                S.append({
                    "run_id": f"{rid:04d}",
                    "id_ctrl":  f"{rid:04d}_con",
                    "id_treat": f"{rid:04d}_treat",
                    # shared
                    "genome": genome,
                    "acc_key": acc_key,
                    "gc_key": gc_key,
                    "fragment_length": frag,
                    "read_length": read,
                    "nb_k": nbk,
                    "aligner": aligner,
                    "peakcaller": peakcaller,
                    "tf_exp": tf_exp,
                    "gc_exp": gc_exp,
                    "acc_exp": acc_exp,
                    "map_coverage_pct": map_coverage_pct,
                    "map_sigma": map_sigma,
                    "map_enrich": map_enrich,
                    "map_exp": map_exp,
                    "use_control": use_control,
                    # per-condition
                    "coverage_ctrl":  cov_c,
                    "coverage_treat": cov_t,
                    "tf_peak_count_ctrl": 0,      # control has no TF peaks
                    "tf_peak_count_treat": tpc,
                    "tf_sigma": tsig,
                    "tf_enrich": tenr,
                })
    return S

SAMPLES = build_samples(config)

# ---------- helpers shared by modules ----------
def find_row(run_id):
    return next(r for r in SAMPLES if r["run_id"] == run_id)

def fasta_path(row):    return config["genome_paths"][row["genome"]]
def acc_bed_path(row):  return config["accessibility_paths"][row["genome"]][row["acc_key"]]
def gc_bias_path(row):  return config["gc_bias_paths"][row["gc_key"]]

def bowtie2_index(row): return config["indexes"][row["genome"]]["bowtie2_index"]
def bwa_index(row):     return config["indexes"][row["genome"]]["bwa_index"]

def macs2_gsize(row):   return config["peakcallers"]["macs2"]["genome_size"][row["genome"]]
def macs2_flags():      return config["peakcallers"]["macs2"].get("flags", "")
def epic2_flags():      return config["peakcallers"].get("epic2", {}).get("flags", "")

# ---------- parameter manifest ----------
rule write_params_table:
    output: config["params_table"]
    run:
        import pandas as pd
        pd.DataFrame(SAMPLES).to_csv(output[0], index=False)

# ---------- include stage modules ----------
include: "rules/simulation.smk"
include: "rules/alignment.smk"
include: "rules/peakcalling.smk"

# ---------- global default target ----------
rule all:
    input:
        rules.sim_done.input,
        rules.align_done.input,
        rules.peaks_done.input,
        config["params_table"]
    default_target: True
