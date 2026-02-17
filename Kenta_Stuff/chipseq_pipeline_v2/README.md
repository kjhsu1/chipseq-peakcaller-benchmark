# Chipseq Pipeline

This directory contains a Snakemake workflow for simulating CHIP-seq reads,
aligning them and calling peaks.

## Layout
- `Snakefile.py` – entry point for the workflow.
- `config.yaml` – base configuration options. Bias exponents (`tf_exp`, `gc_exp`, `acc_exp`) allow reshaping TF, GC and accessibility PMFs.
- `all_comb_alg_peak.config.yaml` – example config exercising every aligner/peakcaller combination.
- `envs/` – Conda environment definitions.
- `rules/` – Snakemake rule files.
- `scripts/` – helper scripts used by the rules.
- `tests/` – unit tests for pipeline utilities.

## Running Tests
```bash
source ../snakemake_stuff/setup.sh
pytest tests
```

## Precision/Recall Stats
Compute precision/recall/F1 for called peaks versus planted centers in archived sweep results.

```bash
python scripts/peak_pr_stats.py \
  --results-dir archived_results/results_tf_gcacc_ctrl_sweep32_clean64
```

This writes a timestamped folder under `archived_results/` with:
- `per_run_stats.csv`
- `group_summary.csv`
- `run_filter_manifest.txt`

## Config Conventions
- `peakcaller_list` controls which callers are included in the sweep.
- `peakcallers` is the caller-parameter dictionary (flags/genome sizes/etc).
- Simulation emits `planted_peaks.bed` for each `{run_id}/{cond}`:
  - treatment contains 1-bp planted TF center intervals
  - control is empty when `tf_peak_count_ctrl = 0`
