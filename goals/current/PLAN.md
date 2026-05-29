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
This production-validation framing coexists with the earlier controlled-wave
work; preserve both additively in the tracking files.

## Historical Controlled + Simulator Plan Context

1. Reproduce and document the current ChIPs realstudy failure points without
   changing code. Done.
2. Repair the ingest-to-alignment handoff so downloaded FASTQs are real
   workflow outputs, not out-of-band assumptions.
3. Split processed-BAM studies from FASTQ-needs-alignment studies so ENCODE BAM
   inputs do not get forced through `data/aligned/...` paths meant for FASTQ
   ingest.
4. Remove stale parse-time metadata assumptions in `Snakefile.py` so regenerated
   manifests and run tables do not silently disagree with the DAG.
5. Unify local path configuration for references and Bowtie2 indexes so local
   dry-runs and local execution use the same override mechanism across ingest
   alignment and ChIPs simulation.
6. Make peak-calling genome size assembly-aware instead of hard-coding one worm
   scale value for every study.
7. Align ChIPs model-learning inputs with the real-study peak BED generation
   logic so the BAM used for `chips learn` matches the replicate structure used
   to define the training peaks.
8. Tighten manifest/local-file bookkeeping so `downloaded` and `local_exists`
   reflect the current checkout accurately and do not overstate local readiness.
9. Future-proof run IDs against parameter-axis growth so later expansions do not
   create silent collisions.
10. Re-run the realstudy test suite and local Snakemake dry-runs, then reassess
   whether the ChIPs path is truly locally runnable or still cluster-only.

1. Record the current goal and initialize progress/log tracking. Done.
2. Document the controlled-sweep bug-fix state as prerequisite context. Done.
3. Write simulator compatibility criteria, compare ChIPs/ChIPulate/isChIP/ChIPsim,
   and choose the simulator for the realistic benchmark path.
   Done.
4. Revise the 11-step plan around ChIPs and the ontology handoff. Done.
5. Add ChIPs workflow config/rules/Slurm launch assets in
   `chipseq_pipeline_v2_realstudy`. Done.
6. Add matched EPIC2 sweep config and Slurm/shell launch assets in
   `chipseq_pipeline_v2`. Done.
7. Add `/goals` operating instructions to root `AGENTS.md`. Done.
8. Prune low-utility active scripts and add a plain-language script index. Done.
9. Run fast checks, tests, and dry-runs; record proof and limitations. Done in
   the local `background_project` Conda environment.

## Next Cluster Step

Run these from the cluster once the branch is available there:

```bash
cd chipseq_pipeline_v2
source ../snakemake_stuff/setup.sh
snakemake -s Snakefile.py --configfile configs/epic2_tfclean_realistic_peaks_wavy_narrow_integrated_288.yaml --dry-run
bash slurm/submit_epic2_tfclean_1728_series.sh
```

```bash
cd chipseq_pipeline_v2_realstudy
source ../snakemake_stuff/setup.sh
snakemake -s Snakefile.py --config enable_chips_targets=true --dry-run
bash slurm/submit_chips_realsim.sh
```

## Local Validation Command Pattern

On this Mac, use:

```bash
eval "$(/opt/anaconda3/bin/conda shell.zsh hook)"
conda activate background_project
HOME=/private/tmp/chipseq_snakemake_home snakemake -s Snakefile.py --dry-run
```

## ChIPs Environment and Tiny Smoke Wave

1. Confirm whether `background_project` exposes `chips`. Done: it did not at
   the start of this wave.
2. Install or expose ChIPs inside `background_project` so workflow commands do
   not depend on another Conda environment. Done.
3. Add or reuse a tiny local smoke config that resolves one ChIPs run only, with
   very small coverage settings and no large sweep targets. Done.
4. Validate in a tight loop:
   - `python -m py_compile scripts/*.py`
   - `pytest tests`
   - Snakemake dry-run for the tiny ChIPs config
   - one tiny ChIPs-controlled output target if the installed binary is usable
   Done.
5. Record proof, limitations, and whether the realstudy pipeline is locally
   runnable beyond dry-run. Done.

## Local and Slurm ChIPs Execution Wave

1. Inspect the current ChIPs configs, docs, and launch scripts. Done.
2. Add a local execution config that uses local staged asset paths instead of
   cluster production paths. Done.
3. Add a small local launcher for dry-run and full-run modes without Slurm.
   Done.
4. Update docs so local Snakemake, tiny smoke, and cluster Slurm commands are
   all explicit and separate. Done.
5. Validate compile/tests plus local launcher dry-runs and tiny smoke execution.
   Done.

## Broad/Narrow-Safe Peakcaller Retry Wave

1. Append this new wave across `GOAL.md`, `PLAN.md`, `PROGRESS.md`, and
   `LOG.md` without deleting earlier history. In progress.
2. Audit the six current TF-clean categories for planted shape, category name,
   current decoding mode, and evaluation mode; record any confounds.
3. Fix the planted-shape versus decoding-mode mismatch so category semantics
   and caller modes agree.
4. Research reputable peakcallers with both broad and narrow support that are
   practical on both the local Mac and the HPC cluster.
5. Choose the best-fit replacement for EPIC2 and document why it is a better
   benchmark fit.
6. Replace the EPIC2 path in `chipseq_pipeline_v2` with the selected caller,
   including configs, Snakemake rules, evaluation helpers, docs, and Slurm
   wrappers.
7. Run the smallest useful local validation loops for the replacement caller in
   both broad-like and narrow-like settings.
8. Tune the replacement caller parameters in tight loops and log retained
   settings plus rejected settings.
9. Scale into staged batch runs, generate six-category plots per batch, and
   inspect for low precision/recall anomalies before escalating further.
10. Recheck the MACS2 six-category rerun path for any remaining planted-shape
    versus decoding mismatch and patch/document it if found.
11. Keep the HPC realstudy continuation commands, status, and blockers current
    while local controlled-pipeline work progresses.

## Current Strategy Shift

The evidence from matched staged controls now changes the local benchmark
strategy:

1. Keep HOMER integration and staged reporting intact as a documented broad +
   narrow replacement candidate that was fully tested locally.
2. Treat MACS2 as the validated leading path for larger reruns because the
   full six-category matched staged control comparison is now complete and
   substantially stronger than HOMER for both broad and narrow categories.
3. Use the completed six-category MACS2-vs-HOMER staged comparison as the
   decision artifact for the next larger rerun stage while preserving the
   planted-shape/decoder guardrails.
4. The next concrete local execution stage is the validated six-config
   `balanced_tfclean_*_288` MACS2 rerun path, using
   `run_balanced_tfclean_288_local.sh` for dry-runs or full local sequential
   execution.
5. Avoid frequent runtime polling during long local runs; prefer milestone
   checks and heartbeat-style resumes when no decision can be made yet.
6. Treat heartbeat-first monitoring as the default for long-running jobs in
   this wave: do one early checkpoint to estimate the next meaningful gate,
   then defer until completion, score-ready output, a clear phase transition,
   or a meaningful failure changes what we can do next.
7. Use `scripts/summarize_balanced_tfclean_288_progress.py` as the preferred
   low-token checkpoint for active local `288` MACS2 roots so heartbeat resumes
   can inspect one compact milestone artifact before deciding whether a scoring
   refresh or a later defer is warranted.
8. Use `refresh_balanced_tfclean_288_local_progress.sh` as the normal operator
   entrypoint for that checkpoint refresh so the local runner, progress
   checkpoint, and score refresher all have consistent selector-driven shell
   wrappers.
9. Use the refreshed `all_six` local-progress artifact as the operator view for
   the full local `1728` wave: it should make not-started vs in-progress vs
   score-ready roots explicit before any decision to launch the remaining four
   configs or refresh scoring.
10. Use `launch_balanced_tfclean_288_remaining_local.sh` when local capacity is
    available for the next wave step so only the `not_started` configs are
    targeted and already-active roots are left alone.
11. Prefer subset-driven remaining launches like `flatearth_pair` or
    `hilly_pair` when we want to advance the local six-config/`1728` wave in
    scientifically coherent chunks rather than immediately launching every
    untouched root.
12. Use bounded remaining launches with `--limit N` when local capacity or
    risk tolerance suggests advancing only the next one or two untouched roots
    from a selected subset.
13. Treat `analysis_outputs/tfclean_balanced_288_local_progress_current/` as
    the authoritative full-wave checkpoint only for `all_six`; subset refreshes
    should use selector-scoped progress artifacts so the six-config operator
    view is not overwritten accidentally.
14. Use the decision artifact derived from the authoritative `all_six`
    progress table to make heartbeat resumes explicit about whether the batch
    should defer, launch more, or score, instead of re-deriving that judgment
    from raw counts each time.
15. Prefer `refresh_balanced_tfclean_288_local_decision.sh` as the normal
    heartbeat checkpoint command for the local six-config/`1728` wave because
    it refreshes both the authoritative progress table and the batch-level
    recommendation together.

## Current Completion Audit State

1. The local controlled peakcaller-retry wave is now evidence-complete through
   the staged `128` caller comparison and the full six-category local MACS2
   `288` rerun wave.
2. The local proof artifacts to treat as authoritative are:
   - `chipseq_pipeline_v2/analysis_outputs/peakcaller_staged_comparison_current/`
   - `chipseq_pipeline_v2/analysis_outputs/tfclean_balanced_288_local_current/`
   - `chipseq_pipeline_v2/analysis_outputs/tfclean_balanced_288_local_decision_current/`
3. The remaining non-local gap is the HPC realstudy continuation front, which
   is still open until cluster execution finishes and produces ontology metrics
   and plots.
4. The current best next move on that open front is to carry forward the
   repaired Bowtie2 basename-prefix config for cluster use, then treat the
   missing production reference assets and the absent post-resubmission cluster
   outputs as the remaining evidence gap rather than as a local workflow-shape
   ambiguity.
