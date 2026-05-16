# Active Scripts Index

This folder contains scripts used by the active controlled benchmark. Historical
or low-utility helpers are archived under
`archive/legacy_scripts/chipseq_pipeline_v2_scripts/`.

## Simulation And Workflow Helpers

- `updated_chip_seq.py`: generates simulated treatment/control reads and truth
  files for the controlled benchmark.
- `seed_helpers.py`: derives deterministic condition-specific seeds.
- `lib.py`: shared small helpers retained from the original workflow.

## Evaluation And Reporting

- `peak_pr_stats.py`: computes per-run precision, recall, F1, and count-based
  summaries from called peaks and truth files.
- `eval_helpers.py`: resolves peak paths and writes metric definition sidecars.
- `control_depth_eval.py`: aggregates archived sweep outputs by category and
  caller.
- `build_balanced_288_config_report.py`: builds the six-category PR/recall/F1
  report from repaired per-run outputs.
- `audit_control_plateaus.py`: checks whether repeated high-control outcomes
  reflect cloned outputs or plausible caller saturation.

## Controlled Cleanup And Scientific Audits

- `summarize_current_tfclean_categories.py`: summarizes the current six
  TF-clean category outputs.
- `category_summary_lib.py`: helpers for category summaries.
- `score_sim_realism.py`: BAM-only realism scorecard entrypoint.
- `realism_metrics_lib.py`: realism metric helpers.
- `summarize_hilly_density_tuning.py`: summarizes hilly-density tuning pilots.
- `hilly_tuning_lib.py`: hilly tuning helper functions.

## Focused Investigation Helpers

- `investigate_peak_recovery.py`: inspects peak recovery behavior for selected
  runs.
- `summarize_peak_recovery_patterns.py`: summarizes peak recovery
  investigations.
- `control_depth_calibrate.py`: older calibration helper still referenced for
  control-depth exploration.
- `rebuild_real_encode_skn1_bigwigs.sh`: one-off rebuild helper for the local
  ENCODE SKN-1 reference assets.
