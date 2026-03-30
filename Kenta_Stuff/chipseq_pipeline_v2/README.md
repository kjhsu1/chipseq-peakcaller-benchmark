# Chipseq Pipeline

This directory contains a Snakemake workflow for simulating CHIP-seq reads,
aligning them and calling peaks.

## Layout
- `Snakefile.py` – entry point for the workflow.
- `config.yaml` – base configuration options. Bias exponents (`tf_exp`, `gc_exp`, `acc_exp`, `map_exp`) allow reshaping TF, GC, accessibility and mappability PMFs.
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

## Control-Depth Sweep (128 runs)
This workflow supports a path-scoped calibration + 4-group sweep process entirely under
`Kenta_Stuff/chipseq_pipeline_v2`.

### 1) Bias calibration pilot
```bash
cd /Users/kentahsu/Code/KorfLab/Background_Forked/Kenta_Stuff/chipseq_pipeline_v2
rm -rf results
conda run -n chipseq_align snakemake -s Snakefile.py --configfile configs/calibration_bias_pilot.yaml --cores 8
conda run -n chipseq_align python scripts/control_depth_calibrate.py \
  --results-dir results \
  --params-csv results/params/run_params.csv \
  --output-dir archived_results/bias_calibration_YYYYMMDD_HHMMSS
```

### 2) Full sweep (4 groups x 32 runs)
Repeat these commands for each config:
- `configs/sweep128_proxy_uniform.yaml`
- `configs/sweep128_proxy_bumpy.yaml`
- `configs/sweep128_ctrltreat_uniform.yaml`
- `configs/sweep128_ctrltreat_bumpy.yaml`

```bash
cd /Users/kentahsu/Code/KorfLab/Background_Forked/Kenta_Stuff/chipseq_pipeline_v2
rm -rf results
conda run -n chipseq_align snakemake -s Snakefile.py --configfile configs/sweep128_proxy_uniform.yaml --cores 8
cp -a results archived_results/sweep128_proxy_uniform_YYYYMMDD_HHMMSS
```

### 3) Cumulative group evaluation and figures
```bash
cd /Users/kentahsu/Code/KorfLab/Background_Forked/Kenta_Stuff/chipseq_pipeline_v2
conda run -n chipseq_align python scripts/control_depth_eval.py \
  --input-dirs \
    archived_results/sweep128_proxy_uniform_YYYYMMDD_HHMMSS \
    archived_results/sweep128_proxy_bumpy_YYYYMMDD_HHMMSS \
    archived_results/sweep128_ctrltreat_uniform_YYYYMMDD_HHMMSS \
    archived_results/sweep128_ctrltreat_bumpy_YYYYMMDD_HHMMSS \
  --output-dir archived_results/eval_sweep128_YYYYMMDD_HHMMSS
```

### Outputs
- `tables/category_method_ratio_summary.csv`
- `tables/category_method_summary.csv`
- `tables/figure_table_manifest.csv`
- per category and method combo:
  - `<category>/<aligner+peakcaller>/figures/pr_f1_vs_ratio.png`
  - `<category>/<aligner+peakcaller>/figures/fdr_inflation_vs_ratio.png`
  - `<category>/<aligner+peakcaller>/figures/interaction_heatmap.png`

## Eight-Category Sweep
The current workflow supports eight peakcaller-based experiment categories:
- `shotgun_flat_narrow_none`
- `shotgun_flat_broad_none`
- `flatearth_peaks_broad_integrated`
- `flatearth_plateaus_broad_integrated`
- `realistic_peaks_wavy_narrow_integrated`
- `realistic_peaks_hilly_narrow_integrated`
- `realistic_plateaus_wavy_broad_integrated`
- `realistic_plateaus_hilly_broad_integrated`

Three configs extend those categories with higher-value enrichment settings:
- `configs/flatearth_broad_integrated_enrich128.yaml`
- `configs/realistic_peaks_integrated_enrich128.yaml`
- `configs/realistic_plateau_integrated_enrich128.yaml`

### TF-vs-Histone Boundary
- `tf_sigma` controls planted signal width.
- `tf_exp` controls contrast/sharpness of the TF PMF after normalization; it does not set the width.
- The simulator accepts any `tf_sigma`, but the benchmark's primary category mapping uses:
  - TF-like / narrow: `tf_sigma <= 5`
  - histone-like / broad: `tf_sigma >= 15`
  - intermediate values (`5 < tf_sigma < 15`) are treated as ambiguous and excluded from the headline 8-category summaries

## Config Conventions
- `peakcaller_list` controls which callers are included in the sweep.
- `peakcallers` is the caller-parameter dictionary (flags/genome sizes/etc).
- `macs2_mode` controls MACS2 decoding strategy (`narrow` or `broad`) and defaults to `["narrow"]`.
- `use_control` toggles whether peakcallers receive control/input BAM (`[true, false]` sweeps both modes).
- MACS2 outputs are normalized to `results/{run_id}/peaks/macs2/{run_id}_peaks.bed` for both modes.
- Mappability bias is optional. `map_coverage_pct` deterministically sets the number of mappability Gaussian centers:
  - `num_map_peaks = round(num_bins * map_coverage_pct / 100)`
  - `map_sigma`, `map_enrich`, and `map_exp` control the Gaussian shape.
- Simulation emits `planted_peaks.bed` for each `{run_id}/{cond}`:
  - treatment contains 1-bp planted TF center intervals
  - control is empty when `tf_peak_count_ctrl = 0`
