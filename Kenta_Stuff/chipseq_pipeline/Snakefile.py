# Master Snakefile

configfile: "config.yaml"

from itertools import product

# ---------- sweep catalog ----------
def build_samples(cfg):
    S = []
    rid = 0
    for genome, acc_key, gc_key, frag, read, nbk, aligner, peakcaller, tf_exp, gc_exp, acc_exp in product(
        cfg["genomes"], cfg["acc_beds"], cfg["gc_bias_sets"],
        cfg["fragment_length"], cfg["read_length"], cfg["nb_k"],
        cfg["aligners"], cfg["peakcallers"],
        cfg["tf_exp"], cfg["gc_exp"], cfg["acc_exp"]
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

# ---------- helpers shared by modules (no leading underscore) ----------
def find_row(run_id):
    return next(r for r in SAMPLES if r["run_id"] == run_id)

def fasta_path(row):    return config["genome_paths"][row["genome"]]
def acc_bed_path(row):  return config["accessibility_paths"][row["genome"]][row["acc_key"]]
def gc_bias_path(row):  return config["gc_bias_paths"][row["gc_key"]]

def bowtie2_index(row): return config["indexes"][row["genome"]]["bowtie2_index"]
def bwa_index(row):     return config["indexes"][row["genome"]]["bwa_index"]

def macs2_gsize(row):   return config["peakcallers"]["macs2"]["genome_size"][row["genome"]]
def macs2_flags():      return config["peakcallers"]["macs2"].get("flags", "")
def epic2_gsize(row):   return config["peakcallers"]["epic2"]["genome_size"][row["genome"]]
def epic2_flags():      return config["peakcallers"]["epic2"].get("flags", "")

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

# ---------- global default target (Pattern B) ----------
rule all:
    input:
        rules.sim_done.input,
        rules.align_done.input,
        rules.peaks_done.input,
        config["params_table"]
    default_target: True

