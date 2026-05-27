# Plan

## Realstudy Production Validation Wave

1. Replace the broad historical goal framing with the new realstudy production
   objective while preserving prior history in `LOG.md`. In progress.
2. Inspect the current realstudy workflow files named in `GOAL.md` and confirm:
   - exact output paths for alignments, peak BEDs, ontology outputs, and summary
     outputs
   - current run-table shape from `metadata/prototype_run_table.csv`
   - current failure points from the latest Slurm attempts
3. Build or update a production validation script:
   - target path:
     `scripts/validate_chips_ontology_production_run.py`
   - inputs:
     `--run-table`, `--output-root`, `--summary-dir`, `--write-report`
   - outputs:
     Markdown report plus machine-readable JSON or CSV summary
   - checks:
     success criteria `1` through `9` and `11` where locally possible
4. Add runtime/resource capture to the production path:
   - prefer post-run logic in
     `slurm/chips_realsim_ontology_128cpu_2tb.sbatch`
   - use `sacct` when available
   - do not mask Snakemake failures
5. Add reproducibility packaging:
   - create `${DEST_ROOT}/reproducibility/`
   - copy exact configs, params table, Slurm scripts, commit hash, validation
     outputs, resource report, and top conclusions
6. Confirm and, if needed, extend ontology summary plotting:
   - required now:
     `ontology_f1_heatmap.png`
   - likely add:
     `control_response_by_ontology.png`
7. Add ontology class coverage validation:
   - infer the correct ontology/class column from generated outputs
   - implement default threshold:
     at least `3` ontology classes with at least `10` completed runs each
   - write coverage summary to summary outputs
8. Run the tight feedback loop:
   1. `python -m py_compile scripts/*.py`
   2. `pytest tests`
   3. production dry-run with `config.yaml` and
      `configs/chips_cluster_full.yaml`
   4. validation script against any existing production output root, if useful
9. Reassess cluster readiness:
   - if local/static checks are clean, identify the exact remaining cluster
     blocker
   - current known blocker to verify first:
     Bowtie2 `ce11` index prefix points to a directory instead of the `ce11`
     basename prefix
10. Finish with an 11-criterion pass/fail table and exact evidence paths.

## Current Strategy

Start by making the goal-tracking files match the new prompt, then inspect the
current realstudy workflow and existing output directories so the new validator
and Slurm provenance logic use the workflow's actual paths instead of guessed
ones.
