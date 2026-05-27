# Goal

Make the `chipseq_pipeline_v2_realstudy` production realstudy ChIPs ontology
workflow runnable, self-verifying, and audit-ready so that one Slurm
production run completes with exit code `0`, produces every expected run output
from `metadata/prototype_run_table.csv`, builds final ontology/performance
summaries and plots, captures runtime/resource provenance, and leaves a
reproducibility package proving exactly what produced the results.

## Context

- Project: `kjhsu1/chipseq-peakcaller-benchmark`
- Branch: `main`
- Working dir: `chipseq_pipeline_v2_realstudy`
- Pipeline: Snakemake workflow in `Snakefile.py`
- Main production launcher:
  - `slurm/submit_chips_realsim_ontology_128cpu_2tb.sh`
  - `slurm/chips_realsim_ontology_128cpu_2tb.sbatch`
- Main config files:
  - `config.yaml`
  - `configs/chips_cluster_full.yaml`
- Target workflow flags:
  - `enable_chips_targets=true`
  - `enable_chips_ontology_targets=true`
- Expected parameter table:
  - `metadata/prototype_run_table.csv`
- Current output families:
  - Simulated outputs: `results_chips/{run_id}/...`
  - Ontology outputs: `analysis_outputs/chips_ontology/{run_id}/...`
  - Summary outputs:
    `analysis_outputs/chips_ontology/summary/...`
- Important existing behavior:
  - The sbatch script already touches `RUN_COMPLETE` on success and
    `RUN_FAILED` on failure.
  - The ontology workflow already defines per-run `region_metrics.csv`,
    per-run `classified.csv`, `combined_region_metrics.csv`, summary CSVs,
    and `ontology_f1_heatmap.png`.
- Audience:
  - Me and future reviewers who need to trust that the production benchmark
    completed, is complete, and is reproducible.

## Success Criteria

All 11 must pass.

1. Successful end-to-end realstudy production run
   - `${DEST_ROOT}/RUN_COMPLETE` exists
   - `${DEST_ROOT}/RUN_FAILED` does not exist
   - Slurm log shows final workflow exit code `0`
2. Complete expected run set produced
   - unique expected `run_id` count from
     `metadata/prototype_run_table.csv` equals unique completed `run_id` count
3. All simulated alignment outputs present
   - for every expected `run_id`, all four Bowtie2 BAM/BAM index files exist
4. All peak-call outputs present
   - for every expected `run_id`,
     `results_chips/{run_id}/peaks/macs2/{run_id}_peaks.bed` exists
5. Ontology scoring complete for all runs
   - for every expected `run_id`,
     `analysis_outputs/chips_ontology/{run_id}/region_metrics.csv` and
     `classified.csv` exist and are readable
6. Combined ontology table built
   - `analysis_outputs/chips_ontology/combined_region_metrics.csv` exists,
     is non-empty, and covers every expected run
7. Final performance summary tables built
   - required summary CSVs exist, are non-empty, and are readable
8. Final performance plots built
   - required PNG list exists, is non-empty, and is readable
9. Coverage of ontology classes is acceptable
   - default threshold: at least `3` ontology classes with at least `10`
     completed runs each, unless documented otherwise
10. Runtime/resource report captured
   - wall time, max RAM, CPU allocation, disk I/O if available, completed job
     count, failed job count, Slurm job ID, output root, and log path are
     written to summary outputs
11. Reproducibility package complete
   - final output root contains exact configs, params table, Slurm scripts,
     commit hash, validation report, resource report, and top conclusions

## Required Files And Commands

- Inspect and update as needed:
  - `Snakefile.py`
  - `rules/alignment.smk`
  - `rules/chips_simulation.smk`
  - `rules/ontology_analysis.smk`
  - `scripts/evaluate_by_region_ontology.py`
  - `slurm/submit_chips_realsim_ontology_128cpu_2tb.sh`
  - `slurm/chips_realsim_ontology_128cpu_2tb.sbatch`
- Fast feedback:
  - `python -m py_compile scripts/*.py`
- Tests:
  - `pytest tests`
- Dry-run:
  - `snakemake -s Snakefile.py --configfile config.yaml configs/chips_cluster_full.yaml --config enable_chips_targets=true enable_chips_ontology_targets=true --dry-run`
- Production:
  - `bash slurm/submit_chips_realsim_ontology_128cpu_2tb.sh`
- Post-run validation:
  - `python scripts/validate_chips_ontology_production_run.py --run-table metadata/prototype_run_table.csv --output-root <DEST_ROOT> --write-report <DEST_ROOT>/reproducibility/validation_report.md`

## Operating Rules

1. Plan first and keep `PLAN.md` current.
2. Work autonomously unless genuinely blocked.
3. Self-verify after every meaningful edit.
4. Debug failures directly and log them in `LOG.md`.
5. Use real workflow paths and real output checks.
6. No placeholders or fake validation.
7. Preserve prior goal history and scientific outputs.
8. Preserve Snakemake exit codes in Slurm post-run logic.
9. Keep success checks measurable with exact files, counts, or pass/fail.
10. Re-read all 11 success criteria before stopping.

## Final Deliverable

At the end, report:

1. PASS/FAIL/NOT RUN for all 11 success criteria with evidence paths
2. files created or modified and their purpose
3. commands run and the important outputs
4. production run result: Slurm job ID, output root, completion markers, exit code
5. expected vs actual run counts across alignments, peaks, ontology, and combined tables
6. final summary CSV and PNG outputs with readability confirmation
7. runtime/resource report contents
8. reproducibility package path and contents
9. decisions made, especially coverage threshold and required plots
10. real remaining limitations or follow-ups
