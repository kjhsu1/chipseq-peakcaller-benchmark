## 2026-05-16 21:49 PDT - ChIPs Environment and Tiny Smoke Wave

- New `/goal` objective: install/expose ChIPs in `background_project` and run
  tight local validation/smoke checks for the realstudy ChIPs pipeline without
  launching large sweeps.
- Starting worktree: `/Users/kentahsu/Code/KorfLab/Background_Forked/.codex_worktrees/chips-realstudy-bugfix-wave`.
- Starting branch: `codex/chips-realstudy-bugfix-wave`.
- Confirmed `background_project` path: `/opt/anaconda3/envs/background_project`.
- Finding: `chips` was not on `PATH` after activating `background_project`.
- Finding: `/opt/anaconda3/envs/alignment/bin/chips` exists but is not a safe
  substitute as called directly; it failed with `_chips: command not found`.
- Decision: attempt a real install into `background_project` via Conda using
  `conda-forge` and `bioconda`.
- Conda install attempt failed because package `chips` is not available for the
  active `osx-arm64` channels. Web/package check indicates Bioconda packaging is
  not native for this Mac architecture.
- Important naming trap: `/opt/anaconda3/envs/alignment/bin/chips` is an EMBOSS
  wrapper, not the ChIPs ChIP-seq simulator. It should not be used for this
  pipeline.
- Cloned upstream ChIPs source from `https://github.com/gymreklab/chips` into
  `/private/tmp/chips-source`.
- Source build initially failed in vendored zlib `1.2.8` because modern macOS
  SDK headers conflicted with zlib's old `TARGET_OS_MAC`/`fdopen` fallback.
- Applied a temporary build-only patch in `/private/tmp/chips-source` to avoid
  that old zlib fallback on `__APPLE__`, then rebuilt successfully.
- Installed ChIPs into `/opt/anaconda3/envs/background_project/bin/chips`.
- Validation: `conda activate background_project && chips --version` prints
  `chips-2.4`.
- Added `configs/chips_tiny_smoke.yaml` for a one-run local smoke with
  `coverage_treat=0.001`, `coverage_ctrl=0.001`, `seed=11`, and ChIPs
  `numcopies_histone=1`.
- Kept local smoke data out of git via `.gitignore`:
  `chipseq_pipeline_v2_realstudy/data/local_smoke/` and
  `chipseq_pipeline_v2_realstudy/results_chips/`.
- Copied the toy `ce11_1pct` FASTA into the realstudy local-smoke data folder
  and indexed it there with `samtools faidx`, avoiding writes to the old
  controlled pipeline.
- Validation: `python -m py_compile scripts/*.py` passed in
  `chipseq_pipeline_v2_realstudy`.
- Validation: `pytest tests` passed in `chipseq_pipeline_v2_realstudy`:
  `17 passed`.
- First tiny Snakemake dry-run attempt failed only because this Snakemake
  version lets `--quiet` consume target-like positional values. Reran without
  `--quiet`.
- Validation: tiny ChIPs dry-run resolved exactly one `chips_simreads_control`
  job.
- Validation: tiny ChIPs actual target completed exactly one
  `chips_simreads_control` job and wrote paired control FASTQs.
- Output proof: each smoke FASTQ has `28` lines, equal to `7` read pairs under
  four-line FASTQ format.
- Validation: default realstudy dry-run still passes after the smoke config and
  ignore-rule additions: Snakemake reports nothing to do for the default
  ingest-prep targets.
- Updated local-runnability assessment: the ChIPs simulator itself is now
  locally usable in `background_project`, and the realstudy pipeline can execute
  a tiny ChIPs WCE/control rule locally. Full local execution is still not the
  default path because it requires full real-study downloads and Bowtie2 index
  assets; cluster execution remains the intended path for production sweeps.

## 2026-05-16 - Local and Slurm ChIPs Execution Wave

- New `/goal` objective: make the ChIPs realstudy Snakemake workflow explicitly
  runnable in both local mode and cluster Slurm mode.
- Current Slurm mode exists through `slurm/submit_chips_realsim.sh` and
  `slurm/chips_realsim_singlejob.sbatch`.
- Current local mode exists only as tests, dry-runs, and the tiny WCE/control
  smoke target. There is not yet a clean full local launcher for staged assets.
- Decision: add a minimal local config plus one shell launcher. This is worth
  adding despite helper-script restraint because it gives a clear public
  interface for local Snakemake execution without changing the workflow rules.
- Added `configs/chips_local_full.yaml`, which points full local runs to staged
  local paths under `data/local/` for references/indexes while preserving the
  existing manifest-driven raw-data paths.
- Added `run_chips_realsim_local.sh`.
  - Default: `MODE=dry-run`.
  - Full local run after staging assets: `MODE=run CORES=4 ./run_chips_realsim_local.sh`.
  - Tiny/local target mode can be selected through `CONFIG_FILES` and `TARGETS`.
- Added local preflight checks for `configs/chips_local_full.yaml` so missing
  full-run assets produce a clear checklist before Snakemake starts.
- Updated `.gitignore` so local full-run data under
  `chipseq_pipeline_v2_realstudy/data/local/` stays out of git.
- Updated `chipseq_pipeline_v2_realstudy/README.md` and
  `docs/chips_realsim_workflow.md` with separate local dry-run, local full-run,
  tiny smoke, and Slurm launch commands.
- Validation: `python -m py_compile scripts/*.py` passed after the local launch
  additions.
- Validation: `pytest tests` passed after the local launch additions:
  `17 passed`.
- Validation: `bash -n run_chips_realsim_local.sh slurm/submit_chips_realsim.sh
  slurm/chips_realsim_singlejob.sbatch` passed.
- Validation: default local launcher now exits with code `3` and a clear
  missing-asset checklist when full local assets are absent:
  `data/local/references/ce11/genome.fa`,
  `data/local/references/ce11/genome.fa.fai`, and
  `data/local/indexes/ce11/bowtie2/ce11*.bt2`.
- Validation: tiny smoke launcher dry-run using
  `CONFIG_FILES='config.yaml configs/chips_tiny_smoke.yaml'` and explicit
  control FASTQ targets passed; Snakemake reported the requested files are
  present and up to date.
- Validation: tiny smoke launcher `MODE=run` with the same explicit targets
  also passed; Snakemake reported nothing to do because the previous smoke
  outputs are already present.

# Log

## ChIPs Ontology Analysis and Cluster Sweep Wave

- Started new `/goal` wave additively per user request; previous `/goals`
  sections are preserved.
- Created isolated worktree and branch:
  `.codex_worktrees/chips-ontology-analysis` on
  `codex/chips-ontology-analysis`, based on current `main`.
- Confirmed current realstudy ChIPs rules include:
  `chips_learn_model`, `chips_simreads_treat`, `chips_simreads_control`,
  `align_bowtie2_chips`, and `call_peaks_macs2_chips`.
- Confirmed ChIPs treatment simulation uses `chips simreads -t bed` with a
  learned model and real peak BED, while control simulation uses
  `chips simreads -t wce`.
- Confirmed ontology scripts already exist as standalone tooling:
  `classify_regions.py`, `ontology_lib.py`, and
  `evaluate_by_region_ontology.py`, but they are not wired into the ChIPs
  Snakemake target graph.
- Confirmed current genome FASTA and Bowtie2 index staging/building are
  prerequisite assets, not Snakemake-built workflow steps.
- Inspected latest successful cluster script patterns, especially
  `chipseq_pipeline_v2/slurm/balanced_tfclean_288_publicgrp_low_128cpu_singlejob.sbatch`
  and `submit_balanced_tfclean_288_series_128cpu.sh`, plus current realstudy
  ChIPs Slurm assets.
- Scientific validity report must be written before implementation and should
  clarify that ChIPs provides realistic simulated reads and input peak
  templates/model parameters, not a closed-form intrinsic PMF.
- Wrote `chipseq_pipeline_v2_realstudy/docs/chips_ontology_scientific_validity.md`
  before workflow implementation. The report concludes the valid downstream
  analysis is empirical template-region recovery: compare final simulated
  MACS2 calls against the real-study peak BED used as the ChIPs treatment
  template, then summarize recovery by region ontology class and control depth.
- Implemented ontology DAG wiring in
  `chipseq_pipeline_v2_realstudy/rules/ontology_analysis.smk`. New stages:
  score template-region recovery, classify regions with the existing ontology
  helper, combine classified region tables, and evaluate aggregate summaries.
- Added `scripts/score_chips_ontology_regions.py` to compute region-level
  treatment/control read counts, binned background variability, enrichment,
  peak-shape proxy metrics, and truth-template recovery labels.
- Added `scripts/combine_csv_tables.py` for a small reusable table-combine
  step instead of embedding brittle shell one-liners in Snakemake.
- Added `enable_chips_ontology_targets` to `config.yaml` and included ontology
  targets in `Snakefile.py` only when that switch is enabled, preserving the
  old ChIPs workflow default behavior.
- Updated ChIPs simread rules so Snakemake reserves the configured ChIPs thread
  count for `chips simreads` instead of passing `--thread` as an untracked
  parameter.
- Added cluster launch assets:
  `chipseq_pipeline_v2_realstudy/slurm/chips_realsim_ontology_128cpu_2tb.sbatch`
  and `chipseq_pipeline_v2_realstudy/slurm/submit_chips_realsim_ontology_128cpu_2tb.sh`.
  The sbatch script requests `128` CPUs and `2000G` RAM, copies the pipeline to
  a per-job work directory, symlinks high-volume outputs into an archived
  results directory, and runs Snakemake in `background_project`.
- Updated `chipseq_pipeline_v2_realstudy/docs/chips_realsim_workflow.md` with
  the ontology target switch, empirical truth-template caveat, and cluster
  submit command.
- Validation proof: `conda run -n background_project python -m py_compile`
  passed for new and touched ontology scripts.
- Validation proof: `conda run -n background_project pytest tests/test_ontology_helpers.py`
  passed: `1 passed`.
- Validation proof: full realstudy test suite passed with
  `conda run -n background_project pytest tests`: `18 passed`.
- Validation proof: ontology-enabled local dry-run passed with
  `HOME=/private/tmp/chipseq_snakemake_home snakemake --dry-run --quiet -j 1
  -s Snakefile.py --configfile config.yaml configs/chips_local_dryrun.yaml
  --config enable_chips_ontology_targets=true`. The resulting DAG has `518`
  jobs total, including `72` score jobs, `72` classify jobs, one combine job,
  and one ontology evaluation job.

## ChIPs Realstudy Bug-Hunt Findings

- Re-ran `chipseq_pipeline_v2_realstudy` tests in `background_project`:
  `15 passed`.
- Re-ran the local ChIPs dry-run with
  `--configfile config.yaml --configfile configs/chips_local_dryrun.yaml --config enable_chips_targets=true --dry-run`.
- The dry-run fails with `MissingInputException` for
  `data/raw/geo_gse67028_celegans_h3k9me2_adult/SRR1917669.fastq.gz`, showing
  that the ingest alignment path expects the FASTQ file directly instead of a
  workflow-produced output.
- `rules/ingest_real_data.smk` currently downloads to a `.done` marker target,
  but `rules/alignment.smk` requires the FASTQ path itself as input. This is a
  real DAG bug, not just a missing local file.
- Selected ENCODE BAM rows (`needs_alignment=false`) are still incompatible with
  the current ChIPs learning and real-study peak-calling path shapes, because
  those rules assume `data/aligned/{study_id}/{role}/aligned.sorted.bam` even
  when the manifest source is a processed BAM download.
- `Snakefile.py` loads `RUNS` and `DATA_MANIFEST` at parse time before rules can
  rebuild `metadata/prototype_run_table.csv` and `metadata/data_manifest.csv`,
  creating stale in-memory workflow state risk.
- `configs/chips_local_dryrun.yaml` only overrides the ChIPs reference/index
  lookup. It does not override the ingest alignment Bowtie2 index path in
  `rules/alignment.smk`, so the local path story is inconsistent.
- `background_project` currently has `snakemake`, `macs2`, `bowtie2`, and
  `samtools`, but not `chips`, so even a logic-fixed workflow would still need
  a local tool install before true local execution.
- The local dry-run helper points to `../chipseq_pipeline_v2/data/genomes/...`,
  but the corresponding Bowtie2 index files are not present under the expected
  local path, reinforcing that the current local dry-run is only a partial
  workflow-shape check.
- The committed manifest says the selected GEO FASTQs are `downloaded` and
  `local_exists=True`, but those FASTQ files are absent in this checkout. That
  bookkeeping mismatch can mislead future runs.
- `sample_reads_from_intensity.py` currently builds `run_id` from only study ID,
  treatment depth, control depth, and seed. That is safe for the current narrow
  axis set, but it will collide if more axes are added later.

## Planned Repair Order

1. Fix the ingest/download/alignment DAG break.
2. Split FASTQ-ingest and processed-BAM study handling.
3. Remove stale parse-time metadata assumptions.
4. Unify local reference/index config wiring.
5. Add assembly-aware MACS2 genome-size handling.
6. Make ChIPs learning inputs match the replicate structure used to call the
   real-study peak BED.
7. Correct manifest bookkeeping and local-runnability documentation.
8. Harden run IDs and re-verify the local ChIPs path.

## ChIPs Realstudy Bug-Fix Implementation Notes

- Reworked `Snakefile.py` so runtime `RUNS` and selected-manifest rows are built
  from source config/manifests instead of trusting stale generated CSVs at parse
  time.
- Reworked the sampling helper layer so the run table builder lives in
  `scripts/realstudy_sampling_lib.py`, and run IDs now include the full active
  axis set (`fragment_length`, `read_length`, `aligner`, `peakcaller`,
  `macs2_mode`) instead of only study/depth/seed.
- Rewired real-study downloading so the download rule remains marker-based for
  Snakemake compatibility, but alignments now depend on download markers while
  consuming the explicit local FASTQ path written by the downloader.
- Reworked real-study BAM routing so processed BAM studies and FASTQ studies no
  longer share the same false `data/aligned/...` assumption.
- Added merged treatment/control BAM outputs for selected real studies, so
  MACS2 real-study peak calling and `chips learn` now use the same replicate
  structure.
- Unified local reference/index lookup through shared assembly helpers in
  `Snakefile.py`, and updated the local dry-run config so the override is not
  ChIPs-only anymore.
- Replaced hard-coded MACS2 genome size with assembly-aware lookup.

## New Validation Proof

- `cd chipseq_pipeline_v2_realstudy && conda activate background_project && python -m py_compile scripts/*.py`
  completed with exit code 0 after the bug-fix edits.
- `cd chipseq_pipeline_v2_realstudy && conda activate background_project && pytest tests`
  passed: `17 passed`.
- `cd chipseq_pipeline_v2_realstudy && conda activate background_project && HOME=/private/tmp/chipseq_snakemake_home snakemake -s Snakefile.py --configfile config.yaml --configfile configs/chips_local_dryrun.yaml --config enable_chips_targets=true --dry-run --quiet`
  passed. Updated dry-run DAG: `372` jobs total, including `3`
  `download_realstudy_file` jobs, `3` downloaded FASTQ alignments, `2`
  `merge_selected_realstudy_bams` jobs, `1` real-study MACS2 peak call, `1`
  `chips_learn_model` job, `72` ChIPs treatment simulations, `72` ChIPs
  control simulations, `144` simulated-read alignments, and `72` simulated
  MACS2 peak calls.
- `cd chipseq_pipeline_v2_realstudy && conda activate background_project && HOME=/private/tmp/chipseq_snakemake_home snakemake -s Snakefile.py --dry-run --quiet`
  still passes for the default ingest-prep path with `2` jobs.
- Local `chips` binary check still returns empty in `background_project`, so
  the workflow is now locally dry-run clean but not yet locally executable
  end-to-end without an additional tool install.

## Initial Context

- Active worktree: `feature-realstudy-benchmark-prototype`.
- Controlled pipeline: `chipseq_pipeline_v2`.
- Realstudy pipeline: `chipseq_pipeline_v2_realstudy`.
- The selected clean worm broad-mark realstudy inputs are downloaded locally and
  tracked through manifests, but raw FASTQs are intentionally ignored.
- Recent controlled-sweep bug fixes are uncommitted at goal start and must be
  treated as prerequisite context before interpreting old plots.

## Simulator Search Notes

- ChIPs is the strongest first-choice candidate because it provides `learn` and
  `simreads`, can infer model parameters from existing BAM+peak inputs, emits
  FASTQ, supports whole-cell extract control mode, and is packaged through
  Bioconda.
- ChIPulate is scientifically rich for TF binding biology but weaker for this
  benchmark because it does not naturally model genome-wide background fragments
  outside bound regions.
- isChIP is command-line and peak-template based, but the ChIPs paper reports
  that it does not infer model parameters from existing datasets.
- ChIPsim is mature Bioconductor infrastructure, but its current emphasis on
  nucleosome positioning makes it less directly aligned with the control-depth
  workflow than ChIPs.

## Repo Hygiene Notes

- Archived low-utility active scripts that were not referenced by current
  workflow rules: PMF plotting/overlay helpers, curated visual pack builder,
  old replot helper, and standalone BAM-to-BigWig helper.

## Implementation Notes

- `GOAL.md` is kept as the user's original `/goal` prompt with only minor
  markdown style cleanup.
- Added root `AGENTS.md` guidance for `/goals` progress updates, helper-script
  restraint, config-parameter restraint, and test placement.
- Added `chipseq_pipeline_v2/docs/chipseq_simulator_comparison.md` with a
  compatibility checklist, simulator comparison table, and ChIPs recommendation.
- Added `chipseq_pipeline_v2/docs/revised_11_step_plan_chips.md` to replace the
  naive real-study pileup-to-PMF plan with a ChIPs-centered realistic simulator
  plan while preserving ontology-window evaluation.
- Added ChIPs workflow rules under
  `chipseq_pipeline_v2_realstudy/rules/chips_simulation.smk`.
- ChIPs targets now default to `geo_gse67028_celegans_h3k9me2_adult`, the clean
  downloaded worm broad-mark study, instead of trying every selected manifest
  study by default.
- Added ChIPs Slurm assets:
  `chipseq_pipeline_v2_realstudy/slurm/chips_realsim_singlejob.sbatch` and
  `chipseq_pipeline_v2_realstudy/slurm/submit_chips_realsim.sh`.
- Added six EPIC2 TF-clean configs mirroring the repaired MACS2 TF-clean
  configs, plus EPIC2 Slurm assets:
  `chipseq_pipeline_v2/slurm/epic2_tfclean_288_singlejob.sbatch` and
  `chipseq_pipeline_v2/slurm/submit_epic2_tfclean_1728_series.sh`.
- EPIC2 local design count check: six configs resolve to 288 runs each, for
  `1728` total runs.
- Replaced the remaining active `dataclasses` import in
  `chipseq_pipeline_v2/scripts/control_depth_eval.py` with a plain
  `typing.NamedTuple`, honoring the no-dataclasses constraint.

## Validation Proof

- Local validation environment: `background_project`
  (`/opt/anaconda3/envs/background_project`), Snakemake `7.32.4`, Python
  `3.10.18`.
- `cd chipseq_pipeline_v2 && source ../snakemake_stuff/setup.sh && python -m py_compile scripts/*.py`
  completed with exit code 0. The setup script printed the expected local Mac
  warning about missing `/cvmfs/.../conda`.
- `cd chipseq_pipeline_v2 && source ../snakemake_stuff/setup.sh && pytest tests`
  passed: `23 passed`.
- `cd chipseq_pipeline_v2_realstudy && source ../snakemake_stuff/setup.sh && python -m py_compile scripts/*.py`
  completed with exit code 0. The setup script printed the same local Mac
  warning.
- `cd chipseq_pipeline_v2_realstudy && source ../snakemake_stuff/setup.sh && pytest tests`
  passed: `15 passed`.
- EPIC2 run-count validation printed six 288-run configs and `total: 1728`.
- ChIPs default target validation from `metadata/prototype_run_table.csv`:
  default ChIPs study `geo_gse67028_celegans_h3k9me2_adult`, `72` run rows
  out of `144` total realstudy prototype rows, with all three treatment and
  control sampling labels represented.
- `cd chipseq_pipeline_v2 && conda activate background_project && HOME=/private/tmp/chipseq_snakemake_home snakemake -s Snakefile.py --configfile configs/epic2_tfclean_realistic_peaks_wavy_narrow_integrated_288.yaml --dry-run --quiet`
  passed. Dry-run DAG: `1442` jobs total, including `576` Bowtie2 alignments,
  `288` EPIC2 peak calls, `288` control simulations, `288` treatment
  simulations, and one params table.
- `cd chipseq_pipeline_v2_realstudy && conda activate background_project && HOME=/private/tmp/chipseq_snakemake_home snakemake -s Snakefile.py --dry-run --quiet`
  passed. Dry-run DAG: `2` jobs total for the default ingest-prep path.
- `cd chipseq_pipeline_v2_realstudy && conda activate background_project && HOME=/private/tmp/chipseq_snakemake_home snakemake -s Snakefile.py --configfile config.yaml configs/chips_local_dryrun.yaml --config enable_chips_targets=true --dry-run --quiet`
  passed. Dry-run DAG: `367` jobs total, including `3` downloaded FASTQ
  alignments, `1` ingested real-study MACS2 peak call, `1` ChIPs learn-model
  job, `72` ChIPs treatment simulations, `72` ChIPs WCE control simulations,
  `144` simulated-read alignments, and `72` simulated MACS2 peak calls.
- Generated six EPIC2 expected run tables under `chipseq_pipeline_v2/results/params/`.
  Each table has `288` rows and `peakcaller=epic2`, for `1728` total expected
  EPIC2 runs.

## Blockers And Limits

- `snakemake_stuff/setup.sh` is written for the cluster and tries
  `/cvmfs/hpc.ucdavis.edu/sw/conda/root/bin/conda`, which is not present on this
  Mac. Local validation therefore uses `conda activate background_project`.
- Snakemake writes its source cache under `~/Library/Caches` by default, which
  is blocked by the sandbox. Local dry-runs pass when
  `HOME=/private/tmp/chipseq_snakemake_home` is set.
- Full EPIC2 and ChIPs execution remains cluster-side because it requires the
  complete bioinformatics toolchain and production references.
