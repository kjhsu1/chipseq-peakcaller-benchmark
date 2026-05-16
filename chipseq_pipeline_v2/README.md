# Chipseq Pipeline

This directory contains a Snakemake workflow for simulating CHIP-seq reads,
aligning them and calling peaks.

## Current Baseline
- The current controlled baseline is the six TF-clean categories archived under
  `analysis_outputs/tfclean_balanced_288_config_report_20260505/`.
- The older eight-category wording is now historical context only; the active
  benchmark baseline for current summaries is the TF-clean six.
- Real-study-conditioned work is separated into the sibling pipeline
  `../chipseq_pipeline_v2_realstudy`.
- Before running tests or helper scripts in either pipeline, execute:

```bash
source ../snakemake_stuff/setup.sh
```

## Layout
- `Snakefile.py` – entry point for the workflow.
- `config.yaml` – base configuration options. Bias exponents (`tf_exp`, `gc_exp`, `acc_exp`, `map_exp`) allow reshaping TF, GC, accessibility and mappability PMFs.
- `all_comb_alg_peak.config.yaml` – example config exercising every aligner/peakcaller combination.
- `envs/` – Conda environment definitions.
- `rules/` – Snakemake rule files.
- `scripts/` – helper scripts used by the rules.
- `tests/` – unit tests for pipeline utilities.
- `docs/epic2_tfclean_sweep.md` – matched EPIC2 follow-up sweep plan and
  cluster launch notes.

## Running Tests
```bash
source ../snakemake_stuff/setup.sh
pytest tests
```

## Benchmark Tracks
The current controlled benchmark baseline is the six-category TF-clean report
under `analysis_outputs/tfclean_balanced_288_config_report_20260505/`.

This pipeline remains the truth-aware controlled simulator track. The
realism-oriented sibling track lives in
`../chipseq_pipeline_v2_realstudy` and should be used for real-study-conditioned
reference modeling, ingest bookkeeping, and ontology-driven follow-up work.

Older notes that assumed the realstudy prototype would live on a separate legacy
feature branch are obsolete. The current implementation path is the sibling
pipeline layout in this worktree.

## Real ENCODE SKN-1 BigWig Rebuild
The repo also carries a lightweight provenance bundle for a local rebuild of
real `C. elegans` SKN-1 ENCODE/modERN signal tracks from experiment
`ENCSR012BJM`.

Tracked files are limited to the rebuild recipe and metadata:
- `scripts/rebuild_real_encode_skn1_bigwigs.sh`
- `real_data_rebuilds/encode_skn1_ce11_pooled_rebuild_20260503_154535/README.md`
- `real_data_rebuilds/encode_skn1_ce11_pooled_rebuild_20260503_154535/logs/`
- `real_data_rebuilds/encode_skn1_ce11_pooled_rebuild_20260503_154535/meta/`

Heavy artifacts remain intentionally untracked:
- downloaded BAMs
- coordinate-sorted / pooled BAMs and `.bai`
- generated BigWigs

Example rebuild invocation:
```bash
cd chipseq_pipeline_v2
conda run -n chipseq_align bash scripts/rebuild_real_encode_skn1_bigwigs.sh \
  --treat1 /path/to/ENCFF775SIN.bam \
  --treat2 /path/to/ENCFF927CRC.bam \
  --control /path/to/ENCFF904DXM.bam \
  --outdir real_data_rebuilds/encode_skn1_ce11_pooled_rebuild_YYYYMMDD_HHMMSS
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
- `group_summary_counts_based.csv`
- `group_summary_mean_of_runs.csv`
- `metric_definition.md`
- `run_filter_manifest.txt`

## Control-Depth Sweep (128 runs)
This workflow supports a path-scoped calibration + 4-group sweep process entirely under
`chipseq_pipeline_v2`.

### 1) Bias calibration pilot
```bash
cd chipseq_pipeline_v2
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
cd chipseq_pipeline_v2
rm -rf results
conda run -n chipseq_align snakemake -s Snakefile.py --configfile configs/sweep128_proxy_uniform.yaml --cores 8
cp -a results archived_results/sweep128_proxy_uniform_YYYYMMDD_HHMMSS
```

### 3) Cumulative group evaluation and figures
```bash
cd chipseq_pipeline_v2
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

### Balanced Control-Coverage Sweep
The balanced sweep configs keep the same parameter-family composition at every
control-coverage x-point and use separate curves for fixed treatment coverage.

Production configs:
- `configs/balanced_flatearth_peaks_broad_integrated_288.yaml`
- `configs/balanced_flatearth_plateaus_broad_integrated_288.yaml`
- `configs/balanced_realistic_peaks_wavy_narrow_integrated_288.yaml`
- `configs/balanced_realistic_peaks_hilly_narrow_integrated_288.yaml`
- `configs/balanced_realistic_plateaus_wavy_broad_integrated_288.yaml`
- `configs/balanced_realistic_plateaus_hilly_broad_integrated_288.yaml`

TF parameter revision:
- The original `balanced_*_288.yaml` configs are preserved for traceability.
- Revised `balanced_tfclean_*_288.yaml` configs use `tf_exp: 1.0` so the
  planted TF PMF is not sharpened after the Gaussian/enrichment model is built.
- Peak-shaped `balanced_tfclean_*` configs use `tf_sigma: 5`; the previous
  `tf_sigma: 1.5` with `tf_exp: 4.0` concentrated most treatment fragments
  into only a few exact fragment-start bins around planted centers, which made
  visual pileups look like repeated duplicates rather than a realistic enriched
  region.
- Peak-shaped `balanced_tfclean_*` configs use `tf_enrich: [1500, 2500]`.
  Plateau `balanced_tfclean_*` configs keep `tf_sigma: 20` and use
  `tf_enrich: [1000, 2500]`.
- Revised production and pilot configs set `emit_pmf_csv: false`; they create
  `pmf.disabled` marker files instead of large PMF tables.
- Hilly configs keep `map_enrich: [10]`; mappability-only tuning put this at
  the mild end of the tested range, with p95 treatment depth about 2x the
  nonzero median background depth.
- Track reruns in `balanced_288_run_attempt_history.log`; archive resulting
  run directories under the Quobyte `archived_results/` area and regenerate
  analysis summaries under `analysis_outputs/`.

Revised production configs:
- `configs/balanced_tfclean_flatearth_peaks_broad_integrated_288.yaml`
- `configs/balanced_tfclean_flatearth_plateaus_broad_integrated_288.yaml`
- `configs/balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288.yaml`
- `configs/balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288.yaml`
- `configs/balanced_tfclean_realistic_plateaus_wavy_broad_integrated_288.yaml`
- `configs/balanced_tfclean_realistic_plateaus_hilly_broad_integrated_288.yaml`

Shape-check pilot configs:
- `configs/pilot_shapecheck_flatearth_peaks_broad_2.yaml`
- `configs/pilot_shapecheck_flatearth_plateaus_broad_2.yaml`
- `configs/pilot_shapecheck_realistic_peaks_wavy_2.yaml`
- `configs/pilot_shapecheck_realistic_peaks_hilly_2.yaml`
- `configs/pilot_shapecheck_realistic_plateaus_wavy_2.yaml`
- `configs/pilot_shapecheck_realistic_plateaus_hilly_2.yaml`

New sweep knobs:
- `seed` can now be swept to create true replicate runs with different planted peaks.
- `tf_seed` optionally fixes planted TF centers independently of other random simulation steps and defaults to `seed`.
- `map_seed` optionally fixes mappability Gaussian centers independently of TF placement and defaults to `seed`.
- production sweeps now derive condition-specific seeds internally so treatment
  and control no longer share the same read or mappability RNG stream.
- `emit_pmf_csv` controls whether simulation writes `pmf.csv` outputs:
  - set `true` for PMF/overlay sanity-check runs
  - set `false` for production sweeps that do not need PMF-derived diagnostics
- `emit_pmf_disabled_marker` controls whether runs with `emit_pmf_csv: false`
  retain `pmf.disabled` marker files:
  - set `true` to keep explicit no-PMF markers
  - set `false` to omit retained PMF artifacts entirely
- `retain_bams` controls whether aligned BAM/BAI files are kept after
  peakcalling:
  - set `true` for visual QC or curated run packs
  - set `false` for production sweeps that only need planted peaks and called peaks
- balanced TF-clean configs now use the MACS2 preset
  `benchmark_control_sensitive_default`; legacy permissive settings remain
  reproducible through the `benchmark_permissive_legacy` preset or explicit
  `peakcallers.macs2.flags`.


### TF-vs-Histone Boundary
- `tf_sigma` controls planted signal width.
- `tf_exp` controls contrast/sharpness of the TF PMF after normalization; it does not set the width.
- The simulator accepts any `tf_sigma`, but the benchmark's primary category mapping uses:
  - TF-like / narrow: `tf_sigma <= 5`
  - histone-like / broad: `tf_sigma >= 15`
  - intermediate values (`5 < tf_sigma < 15`) are treated as ambiguous and excluded from the headline 8-category summaries

## Sequential Slurm Sweeps
The `slurm/` directory contains sequential Slurm submission assets for the two
128-run enrich configs:
- `realistic_peaks_integrated_enrich128.sbatch`
- `realistic_plateau_integrated_enrich128.sbatch`
- `submit_two_sweeps.sh`

These jobs are designed to:
- run one after the other with `afterok`
- use `4` CPUs, `4G` RAM, and `3:00:00` walltime per job
- write all results under `/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2`
- use a Conda env at `/quobyte/ikorfgrp/home/kjhsu/envs/chipseq_pipeline_v2_peak_for_hpc`

Each job executes Snakemake from a job-specific Quobyte work directory so the
workflow's relative `results/` outputs never write into the repo checkout.

For the balanced 288-run configs, the sequential submitters are:
- `submit_balanced_288_series.sh` for the full six-config chain
- `submit_balanced_288_hilly_tail.sh` for rerunning only the remaining hilly pair
- `submit_balanced_288_hilly_tail_25cpu.sh` for rerunning the remaining hilly pair with `25` CPUs passed to both Slurm and Snakemake
- `submit_balanced_tfclean_288_series_20cpu.sh` for the revised TF-clean
  six-config chain, using `20` CPUs and ordering jobs from fastest to slowest
- `submit_balanced_tfclean_288_series_128cpu.sh` for the revised TF-clean
  six-config chain with direct-to-archive `results/` output and `128` CPUs
- `balanced_tfclean_288_publicgrp_8cpu_singlejob.sbatch` for running the six
  revised TF-clean configs sequentially inside one `publicgrp` allocation,
  which is capped at `8` CPUs per job on `publicgrp-high-qos`, so the queue
  wait only happens once
- `balanced_tfclean_288_publicgrp_low_128cpu_singlejob.sbatch` for the same
  six-config single-job run on `publicgrp` `low`, using per-config workdir
  results plus `mv` into `archived_results` on success
- `submit_tfclean_benchmark_publicgrp_low.sh` submits one representative
  TF-clean config at `8`, `16`, and `32` CPUs on `publicgrp` `low` so the
  useful scaling point can be measured before launching a full six-config rerun
- `submit_balanced_tfclean_288_singlejob_scratch.sh` submits the current
  preferred six-config TF-clean rerun: one `ikorfgrp` `high` allocation at
  `25` CPUs, each config runs sequentially on local `/scratch`, and completed
  config results are copied to Quobyte in the background while the next config
  computes
- `submit_epic2_tfclean_1728_series.sh` submits the matched six-config EPIC2
  comparison sweep. This is cluster-only; local validation is limited to config
  parsing and run-table checks unless Snakemake/EPIC2 are installed locally.

The balanced sbatch jobs now sanitize copied Python bytecode caches in the
job-specific pipeline workdir and set `PYTHONDONTWRITEBYTECODE=1` before
launching Snakemake, which avoids stale or partially written `.pyc` files from
breaking parallel `python -m scripts.updated_chip_seq` imports.

Balanced-report helpers:
- `scripts/peak_pr_stats.py` computes per-run precision/recall/F1 for one archived config directory
- `scripts/build_balanced_288_config_report.py` consumes those per-run stats and writes one folder per balanced config with:
  - `pr_recall_f1_vs_ctrl_coverage.png`
  - `plot_point_summary.csv`
  - `data_info.md`
  - plus a root `README.md` and optional copied `attempt_history.log`
- `scripts/investigate_peak_recovery.py` traces each run back to planted peak centers and called peak intervals
- `scripts/summarize_peak_recovery_patterns.py` summarizes false-negative and false-positive seed/locus patterns from those trace outputs
- `balanced_288_run_attempt_history.log` records the known Slurm/log attempt history for the six balanced configs and their retries


## Config Conventions
- `peakcaller_list` controls which callers are included in the sweep.
- `peakcallers` is the caller-parameter dictionary (flags/genome sizes/etc).
- `macs2_mode` controls MACS2 decoding strategy (`narrow` or `broad`) and defaults to `["narrow"]`.
- `use_control` toggles whether peakcallers receive control/input BAM (`[true, false]` sweeps both modes).
- `seed` is the shared fallback seed for the simulator.
- `tf_seed` and `map_seed` can be added to configs or sweeps; if omitted they inherit from `seed`.
- TF planting, mappability-center placement, and read-count noise now use separate RNG streams, so fixing `tf_seed` keeps planted TF centers aligned across categories with the same genome and `tf_peak_count_treat`.
- MACS2 outputs are normalized to `results/{run_id}/peaks/macs2/{run_id}_peaks.bed` for both modes.
- Mappability bias is optional. `map_coverage_pct` deterministically sets the number of mappability Gaussian centers:
  - `num_map_peaks = round(num_bins * map_coverage_pct / 100)`
  - `map_sigma`, `map_enrich`, and `map_exp` control the Gaussian shape.
- Simulation emits `planted_peaks.bed` for each `{run_id}/{cond}`:
  - treatment contains 1-bp planted TF center intervals
  - control is empty when `tf_peak_count_ctrl = 0` and is treated as a temporary intermediate
- Simulation FASTAs (`reads_R1.fasta`, `reads_R2.fasta`) are temporary
  intermediates and are removed after alignment consumes them.
- MACS2 retains the normalized `results/{run_id}/peaks/macs2/{run_id}_peaks.bed`
  plus the native broad/narrow peak file; the extra `.xls` and `.gappedPeak`
  side products are pruned after normalization.
