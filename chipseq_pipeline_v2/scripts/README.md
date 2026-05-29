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
- `summarize_balanced_tfclean_288_progress.py`: writes a compact milestone
  summary for selected local `results_balanced_tfclean_*_288` MACS2 roots so
  heartbeat-style progress checks can use one artifact instead of repeated
  manual directory inspection. It now also reports canonical selectors,
  `launch_state` (`not_started`, `in_progress`, `score_ready`, `scored`), and
  a coarse `recommended_action` so the full local six-config/`1728` wave can
  be managed from one status table. The default `all_six` call writes the
  canonical current artifact, while subset calls now write selector-scoped
  artifact roots so they do not overwrite the authoritative full-wave view.
- `build_balanced_tfclean_288_decision_summary.py`: reads the authoritative
  `all_six` progress CSV and writes a heartbeat-friendly decision summary with
  one batch-level next action plus a compact per-root snapshot.
- `build_homer_staged_report.py`: builds the current staged HOMER 128-run
  combined report from whichever category stats directories have finished so
  far, and can be made strict with `--require-all`.
- `summarize_homer_staged_categories.py`: writes a compact staged HOMER summary
  across the completed category stats roots and records which categories are
  still missing.
- `compare_peakcaller_staged_categories.py`: writes a paired HOMER-vs-MACS2
  staged category comparison so the caller decision can be refreshed as each
  matched stats directory finishes.
- `../score_hilly_macs2_controls.sh`: guarded local entrypoint for scoring the
  two remaining hilly MACS2 staged control categories once all `128` peak BEDs
  exist for each category.
- `../run_balanced_tfclean_288_local.sh`: local sequential entrypoint for the
  validated six `balanced_tfclean_*_288.yaml` MACS2 configs; defaults to
  `--dry-run` and can be switched to full execution with `--run`. It accepts
  either raw config paths, canonical selectors like `wavy_peak_narrow`, or
  batch selectors like `wavy_pair`, `hilly_pair`, `flatearth_pair`, and
  `all_six`.
  If a selected result root is already running locally, Snakemake may report
  the existing directory lock during a second dry-run against the same root.
- `../refresh_balanced_tfclean_288_local_report.sh`: guarded local scorer and
  report refresher for completed `results_balanced_tfclean_*_288` MACS2 runs;
  it rebuilds per-config stats and the combined local 288 report once each
  selected config has all `288` peak BEDs. It also accepts canonical selectors
  like `hilly_plateau_broad`, batch selectors like `wavy_pair`, and now emits a
  clean not-ready message if the corresponding local result root does not exist
  yet.
- `../refresh_balanced_tfclean_288_local_progress.sh`: one-command wrapper for
  `summarize_balanced_tfclean_288_progress.py`; it refreshes the compact local
  `288` milestone artifact for selectors like `wavy_pair`, `hilly_pair`, or
  `all_six` using the local `background_project` environment.
- `../refresh_balanced_tfclean_288_local_decision.sh`: one-command wrapper for
  the authoritative `all_six` progress refresh plus
  `build_balanced_tfclean_288_decision_summary.py`; it updates both the full
  local-wave checkpoint and the heartbeat-friendly batch recommendation in one
  step.
- `../launch_balanced_tfclean_288_remaining_local.sh`: reads the refreshed
  local `288` progress artifact, selects only configs whose `launch_state` is
  `not_started`, and forwards those canonical selectors into
  `run_balanced_tfclean_288_local.sh`. It defaults to `--dry-run` and supports
  `--run` for real execution. It can target the whole wave or a subset such as
  `flatearth_pair` or `hilly_pair`, and then intersects that subset with the
  currently `not_started` roots. `--limit N` can further cap the number of
  untouched configs launched from that ordered remaining set.
- `../refresh_homer_staged_artifacts.sh`: one-command refresh entrypoint for
  the staged HOMER combined report, staged HOMER summary, and staged
  HOMER-vs-MACS2 comparison artifacts.
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
