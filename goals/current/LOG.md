## 2026-05-26 - Realstudy Production Validation Wave

- New `/goal` objective loaded from `~/Downloads/prompt.txt`. This objective
  replaces the older broad benchmark-wave focus for active work, but prior
  entries are preserved as historical context and prior implementation proof.
- Decision: keep `goals/current/` as the tracking location, since that is the
  existing repo convention named in `AGENTS.md`.
- Rewrote `GOAL.md`, `PLAN.md`, and `PROGRESS.md` to match the new production
  objective: one audit-ready realstudy ChIPs ontology run with explicit
  validation, resource capture, and reproducibility packaging.
- Carry-forward blocker from the most recent realstudy Slurm failure:
  Bowtie2 `ce11` index configuration points at a directory-like prefix
  (`indexes/ce11/bowtie2_index`) rather than the actual basename prefix that
  Bowtie2 expects. The likely correct resolved prefix is
  `/quobyte/ikorfgrp/home/kjhsu/results/reference_assets/ce11/bowtie2/ce11`.
- Carry-forward non-blocker: the prior `250G` realstudy failure was not a RAM
  issue. Accounting showed roughly `10.64G` max memory used before failure.
- Next validation order under the new objective:
  1. inspect the live realstudy workflow files named in the prompt
  2. confirm exact output paths from code
  3. add production validation and provenance capture
  4. rerun compile, tests, and dry-run
- Live workflow inspection confirmed:
  - expected simulated BAM outputs:
    `results_chips/{run_id}/bowtie2/{cond}/aligned.sorted.bam(.bai)`
  - expected final peak BEDs:
    `results_chips/{run_id}/peaks/macs2/{run_id}_peaks.bed`
  - expected per-run ontology outputs:
    `analysis_outputs/chips_ontology/{run_id}/region_metrics.csv` and
    `classified.csv`
  - expected summary outputs under:
    `analysis_outputs/chips_ontology/summary/`
- Fixed the production Bowtie2 prefix bug in tracked config:
  `chips.bowtie2_index_by_assembly.ce11` now points to
  `indexes/ce11/bowtie2_index/ce11` instead of the directory-only path.
- Added required production validation/provenance pieces:
  - `scripts/validate_chips_ontology_production_run.py`
  - `control_response_by_ontology.png`
  - `ontology_class_coverage.csv`
  - Slurm post-run resource reports and reproducibility package creation
- Validation proof in the new worktree
  (`codex/realstudy-production-validation`):
  - `source ../snakemake_stuff/setup.sh && python -m py_compile scripts/*.py`
    passed
  - `source ../snakemake_stuff/setup.sh && pytest tests`
    passed: `20 passed`
  - First production dry-run attempt failed for environment reasons only:
    Snakemake tried to write its source cache into a read-only OnDemand home
    path. Reran with writable `HOME=/tmp/chipseq_snakemake_home` and
    `XDG_CACHE_HOME=/tmp/chipseq_snakemake_cache`.
  - Second production dry-run attempt failed because the new worktree did not
    yet have the local untracked `references/` and `indexes/` symlink trees
    copied over from the older checkout.
  - Copied the local untracked asset-link trees from the older checkout into the
    new worktree so the production config could be exercised locally without
    changing tracked workflow code.
  - Production dry-run then passed:
    `snakemake -s Snakefile.py --configfile config.yaml configs/chips_cluster_full.yaml --config enable_chips_targets=true enable_chips_ontology_targets=true --dry-run`
    resolved a `1029`-job DAG.
- Production launch executed from the new worktree:
  `bash slurm/submit_chips_realsim_ontology_128cpu_2tb.sh`
  submitted Slurm job `14516460`.
- Live production output root:
  `/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2_realstudy/chips_runs/chips_realsim_ontology_128cpu_2tb_20260526_034056`
- Live run proof so far:
  - the run passed the old failure points and reached active
    `align_bowtie2_chips` execution
  - latest observed Snakemake progress: `206 / 1029` steps
  - latest observed artifact counts:
    - `76` aligned BAMs
    - `0` peak BEDs
    - `0` `region_metrics.csv`
    - `0` `classified.csv`
  - markers still show:
    - `RUN_COMPLETE=no`
    - `RUN_FAILED=no`
- Later production failure root cause: the cluster cannot resolve external
  download hosts for selected realstudy inputs (`www.encodeproject.org` and
  `ftp.sra.ebi.ac.uk`). The run failed in `download_realstudy_file`, not in
  ontology scoring or ChIPs simulation.
- Decision after user direction: stop spending time on live cluster fetches.
  Production now defaults to pre-staged realstudy inputs via
  `realstudy_require_local_inputs: true`; the ingest rule passes
  `--require-local` and fails fast if the manifest-listed `data/raw/` file is
  missing.
- Validation after the local-input change:
  - `pytest tests` passed: `23 passed`
  - `python -m py_compile scripts/*.py` passed with
    `PYTHONPYCACHEPREFIX=/tmp/chipseq_pycache`
  - production dry-run rebuilt the `1029`-job DAG and then stopped on an
    existing tracked/generated `metadata/data_manifest.csv` protected-output
    condition in this local worktree; the new local-input params hook itself
    was accepted after the wildcard signature fix.

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

## 2026-05-26 PDT - Broad/Narrow-Safe Peakcaller Retry Wave

- Started a new additive `/goals` wave from a fresh worktree on branch
  `codex/peakcaller-retry-sweep`, based on `main`, per the default branching
  rule.
- User objective for this wave:
  1. Continue `/goals` additively.
  2. Replace the failed EPIC2 comparison path with a better peakcaller that
     supports both broad and narrow decoding.
  3. Validate and tune locally in small loops and staged batches.
  4. Fix/document any category where planted peak shape and decoding mode do
     not match.
- Initial controlled-pipeline audit findings:
  - `chipseq_pipeline_v2/docs/category_semantics.md` explicitly admits that some
    category differences are confounded with MACS2 broad versus narrow mode.
  - `chipseq_pipeline_v2/docs/nomenclature_table.md` maps
    `balanced_tfclean_flatearth_peaks_broad_integrated_288` to canonical name
    `flatearth_peak_broad`, meaning a peak-like shape is currently paired with
    broad decoding by design.
  - `chipseq_pipeline_v2/configs/balanced_tfclean_flatearth_peaks_broad_integrated_288.yaml`
    confirms that this peak-like category is currently configured with
    `macs2_mode: ["broad"]`.
  - `chipseq_pipeline_v2/scripts/peak_pr_stats.py` switches truth evaluation to
    interval mode for any run with `macs2_mode == "broad"` or `peakcaller ==
    "epic2"`, so the decoding-mode choice also changes the truth/evaluation
    regime.
  - Working hypothesis: the user’s remembered confound is real; at least one
    peak-like category is being decoded with broad mode, and the workflow/docs
    currently encode that confound instead of preventing it.
- Next actions:
  1. Finish the category/mode audit.
  2. Research candidate peakcallers with both broad and narrow support.
  3. Pick the best-fit replacement and wire it into the controlled workflow.
- Replacement-caller research/result:
  - HOMER was chosen as the first real replacement candidate over EPIC2 because
    it provides explicit TF-style (`findPeaks -style factor`) and histone-style
    (`findPeaks -style histone`) modes, is Bioconda-installable on `osx-arm64`,
    and also has standard Linux/HPC installation paths.
  - Installed `homer-5.1` successfully into `background_project`.
- Initial HOMER integration work:
  - Added `homer` to `chipseq_pipeline_v2/config.yaml` caller options.
  - Added HOMER peakcalling rules using `makeTagDirectory`, `findPeaks`, and
    `pos2bed.pl`.
  - Added HOMER path resolution to `scripts/eval_helpers.py`.
  - Added narrow and broad HOMER pilot configs for local validation.
  - Flipped `balanced_tfclean_flatearth_peaks_broad_integrated_288` from
    `macs2_mode: ["broad"]` to `macs2_mode: ["narrow"]` as the first direct
    fix for the category/decoder confound.
- Validation proof:
  - `python -m py_compile chipseq_pipeline_v2/scripts/*.py` passed after the
    HOMER integration edits.
  - HOMER narrow pilot dry-run passed with an 8-job DAG:
    `simulate -> align -> make_homer_tagdirs -> call_peaks_homer`.
  - HOMER broad pilot dry-run also passed with the same DAG shape.
  - Actual HOMER narrow pilot executed successfully after building the missing
    local Bowtie2 index under
    `chipseq_pipeline_v2/data/indexes/ce11_1pct/bowtie2_index/`.
  - The run produced `results/0001/peaks/homer/0001_peaks.txt` and
    `results/0001/peaks/homer/0001_peaks.bed`.
- Important findings from the first actual HOMER run:
  - Using `makeTagDirectory ... -fragLength pe` caused HOMER to report a bogus
    fragment length (`-123456`) even though peak calling still completed.
  - Manual test without `-fragLength pe` gave a sensible fragment estimate
    (`148`), so the workflow was patched to drop that flag.
  - The current one-run pilot configs all emit to `results/0001/...`, so a
    narrow pilot and a broad pilot will collide unless results are cleaned or
    the workflow gains an explicit result-root/run-namespace mechanism. This is
    now a concrete validation-loop blocker for keeping multiple local pilot
    outputs side by side.

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

## 2026-05-18 Overnight Result Check For 64 CPU / 50G Jobs

- Live `squeue` no longer lists jobs `14193045`, `14193046`, or `14193047`.
- EPIC2 hilly peaks job `14193045` completed successfully:
  - Slurm state: `COMPLETED`, `ExitCode=0:0`
  - ran on `publicgrp/low` with `64` cores and `50G`
  - start/end: `2026-05-17T22:04:46` to `2026-05-17T23:10:05`
  - elapsed: `01:05:19`
  - max memory used: `8.97G`
  - result root:
    `/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/archived_results/epic2_tfclean_realistic_peaks_hilly_narrow_integrated_288_20260517_220447`
  - output counts: `288` EPIC2 `*_domains.bed`, `576` aligned BAMs
- EPIC2 hilly plateaus job `14193046` completed successfully:
  - Slurm state: `COMPLETED`, `ExitCode=0:0`
  - ran on `publicgrp/low` with `64` cores and `50G`
  - start/end: `2026-05-17T23:10:35` to `2026-05-18T00:18:51`
  - elapsed: `01:08:16`
  - max memory used: `8.83G`
  - result root:
    `/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/archived_results/epic2_tfclean_realistic_plateaus_hilly_broad_integrated_288_20260517_231037`
  - output counts: `288` EPIC2 `*_domains.bed`, `576` aligned BAMs
- Realstudy ontology job `14193047` failed:
  - Slurm state: `FAILED`, `ExitCode=1:0`
  - ran on `publicgrp/low` with `64` cores and `50G`
  - start/end: `2026-05-17T21:28:47` to `2026-05-17T23:42:25`
  - elapsed: `02:13:38`
  - max memory used: `50.00G`, exactly the reserved amount
  - result root:
    `/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2_realstudy/chips_runs/chips_realsim_ontology_128cpu_2tb_20260517_212849`
  - marker: `RUN_FAILED`
  - progress reached `278 of 1029 steps (27%) done`
  - failure traceback ended with
    `NotADirectoryError: [Errno 20] Not a directory: 'analysis_outputs/chips_ontology'`
    during Snakemake cleanup, after an upstream job failure.
  - resource signal: `Max Mem used = 50.00G`, so the 50G allocation was fully
    consumed and likely too small for this realstudy shape.
- Started EPIC2 metrics generation for the completed hilly-peaks config after
  confirming both remaining EPIC2 output blocks exist.

## 2026-05-18 Resume Follow-Up

- Rechecked the failed realstudy job `14193047` in the persisted work directory:
  `/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2_realstudy/work/14193047/pipeline`.
- The first hard failure was not the final Snakemake cleanup traceback. The
  upstream failed rule was `prepare_realstudy_ingest_plan`:
  `ModuleNotFoundError: No module named 'scripts'` from
  `python scripts/fetch_real_study_data.py`.
- The cleanup traceback
  `NotADirectoryError: [Errno 20] Not a directory: 'analysis_outputs/chips_ontology'`
  is still recorded, but it happened after the upstream rule failure.
- The resource signal remains important: job `14193047` used exactly `50.00G`
  max memory, so a second failure may still require a larger RAM allocation.
  For this resume, the requested `64` CPU / `50G` setting was preserved.
- Patched realstudy package-style script invocations:
  `python -m scripts.fetch_real_study_data`,
  `python -m scripts.download_real_study_file`,
  `python -m scripts.sample_reads_from_intensity`, and
  `python -m scripts.classify_regions`.
- Restored active cluster submission assets from the staged successful-submit
  pattern:
  `chipseq_pipeline_v2_realstudy/configs/chips_cluster_full.yaml`, staged-copy
  submission in `slurm/submit_chips_realsim_ontology_128cpu_2tb.sh`, and the
  active Slurm script defaults for `publicgrp/low`, `64` CPUs, `50G`, and
  `2-00:00:00`.
- Validation after patch:
  `cd chipseq_pipeline_v2_realstudy && source ../snakemake_stuff/setup.sh && python -m py_compile scripts/*.py && pytest tests`
  passed with `18 passed`.
- Validation after patch:
  `cd chipseq_pipeline_v2_realstudy && source ../snakemake_stuff/setup.sh && XDG_CACHE_HOME=/tmp/chipseq_snakemake_cache snakemake -s Snakefile.py --configfile config.yaml configs/chips_cluster_full.yaml --config enable_chips_targets=true enable_chips_ontology_targets=true --dry-run --quiet`
  passed. DAG size: `1029` jobs.
- Submitted patched realstudy ontology job `14200626`.
  Slurm verification: `publicgrp/low`, `64` CPUs, `50G`, `2-00:00:00`, state
  `PENDING` with reason `Priority` at submit check.
- EPIC2 hilly peaks metrics completed at:
  `/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/archived_results/epic2_tfclean_realistic_peaks_hilly_narrow_integrated_288_20260517_220447/analysis_outputs/peak_pr_stats`.
- EPIC2 hilly plateaus metrics completed at:
  `/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/archived_results/epic2_tfclean_realistic_plateaus_hilly_broad_integrated_288_20260517_231037/analysis_outputs/peak_pr_stats`.
- Hilly plateaus counts-based summary shows F1 increasing with control depth
  from `0.1512` at `0.5x` control to `0.3512` at `24x` control.
- Built the EPIC2 six-config report root:
  `/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/archived_results/epic2_tfclean_six_config_report_20260518_1053`.
  The report root contains `README.md` and `attempt_history.log`; per-config
  plots and `data_info.md` were written into each config's
  `analysis_outputs/peak_pr_stats` directory by the report builder.
- Validation after EPIC2 metrics/report work:
  `cd chipseq_pipeline_v2 && source ../snakemake_stuff/setup.sh && python -m py_compile scripts/*.py && pytest tests`
  passed with `23 passed`.
- At user direction, canceled pending realstudy job `14200626` and resubmitted
  with the same `64` CPUs, `publicgrp/low`, and `2-00:00:00` walltime but
  higher memory: `250G`.
- New realstudy job: `14201143`.
  Slurm verification shows `ReqTRES=cpu=64,mem=250G,node=1,billing=64`,
  state `PENDING`, reason `Priority`.
- New staged submission source:
  `/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2_realstudy/submission_sources/chips_realsim_ontology_submit_20260518_110844`.

## 2026-05-18 Realstudy 250G Failure Investigation

- Realstudy job `14201143` failed after `01:30:36`, from
  `2026-05-18T12:39:55` to `2026-05-18T14:10:31`.
- This was not a RAM exhaustion failure. Slurm summary reports `250G`
  reserved and `10.64G` max memory used.
- The first workflow failure was `align_bowtie2_chips`, not the later
  Snakemake cleanup traceback.
- Root error from the run log:
  `(ERR): "indexes/ce11/bowtie2_index" does not exist or is not a Bowtie 2 index`.
- Follow-on pipe errors were:
  `[main_samview] fail to read the header from "-"` and
  `samtools sort: failed to read header from "-"`.
- The active/staged path `indexes/ce11/bowtie2_index` is a symlink to the
  directory `/quobyte/ikorfgrp/home/kjhsu/results/reference_assets/ce11/bowtie2`,
  whose index files are named `ce11.1.bt2`, `ce11.2.bt2`, etc.
- Bowtie2 `-x` needs the basename prefix, so the correct prefix is likely
  `/quobyte/ikorfgrp/home/kjhsu/results/reference_assets/ce11/bowtie2/ce11`
  or an equivalent repo-local symlink/prefix, not the directory symlink itself.

## 2026-05-26 Broad/Narrow-Safe Peakcaller Retry Wave Checkpoint

- New wave is running on branch `codex/peakcaller-retry-sweep` in worktree
  `/Users/kentahsu/Code/KorfLab/Background_Forked/.codex_worktrees/peakcaller-retry-sweep`.
- Constraint preserved: prior `/goals` history remains intact; new work is
  appended only.
- Confirmed a real category/decoder confound in the controlled benchmark:
  canonical `flatearth_peak_broad` is still a peak-like category
  (`tf_sigma=5`) but was configured with `macs2_mode: broad`, and the
  evaluator also changes truth mode from the same setting.
- Replacement-caller scan favored HOMER over the previous EPIC2 path because
  HOMER supports both:
  - `findPeaks -style factor` for narrow/factor-like decoding
  - `findPeaks -style histone` for broad/histone-like decoding
- Installed `homer-5.1` into `background_project` and verified `findPeaks`,
  `makeTagDirectory`, and `pos2bed.pl` are available on `PATH`.
- Integrated the first HOMER path into `chipseq_pipeline_v2`:
  - added HOMER config entries under `peakcallers.homer`
  - added Snakemake rules for tag-directory creation, peak calling, and BED
    export
  - added HOMER peak-path support in evaluation helpers
  - switched the active alternative caller list entry from EPIC2 to HOMER
- First direct confound fix applied:
  `configs/balanced_tfclean_flatearth_peaks_broad_integrated_288.yaml` now
  uses `macs2_mode: ["narrow"]` instead of `["broad"]`.
- Validation proof:
  `python -m py_compile chipseq_pipeline_v2/scripts/*.py` passed after the
  HOMER integration edits.
- Validation proof: HOMER narrow pilot dry-run passed for
  `configs/pilot_tfclean_homer_realistic_peaks_wavy_narrow_integrated_1.yaml`.
- Validation proof: HOMER broad pilot dry-run passed for
  `configs/pilot_tfclean_homer_realistic_plateaus_wavy_broad_integrated_1.yaml`.
- Initial actual narrow pilot succeeded, but it exposed a workflow-shape
  problem: both one-run pilots wrote to `results/0001/...`, so broad and narrow
  validation outputs would collide.
- Decision: add config-driven `result_root` support to the controlled workflow
  and namespace the pilot outputs instead of inventing fake extra run IDs.
- Implemented namespaced pilot roots:
  - narrow pilot now writes to `results_homer_narrow/...`
  - broad pilot now writes to `results_homer_broad/...`
  - params tables are namespaced under those same roots
- Built the missing local `ce11_1pct` Bowtie2 index in this worktree so actual
  local HOMER pilots could run without waiting on cluster assets.
- Important HOMER tuning finding:
  `makeTagDirectory -fragLength pe` produced a bogus fragment length
  (`-123456`) on our paired-end SAM inputs. Removing that flag yielded a
  sensible estimate and stable downstream peak calling.
- Validation proof: actual broad-mode HOMER pilot completed successfully under
  `results_homer_broad/0001/...` and produced
  `results_homer_broad/0001/peaks/homer/0001_peaks.bed`.
- Broad pilot first metric summary from
  `analysis_outputs/homer_pilot_wavy_plateau_broad/per_run_stats.csv`:
  - `precision=0.9524`
  - `recall=0.6667`
  - `f1=0.7843`
  - `total_called=21`
  - `tp_called=20`
  - interpretation: first broad-mode smoke behavior is strong enough to keep
    HOMER in play
- Validation proof: actual narrow-mode HOMER pilot completed successfully under
  `results_homer_narrow/0001/...` and produced
  `results_homer_narrow/0001/peaks/homer/0001_peaks.bed`.
- Narrow pilot first metric summary from
  `analysis_outputs/homer_pilot_wavy_peak_narrow/per_run_stats.csv`:
  - `precision=0.5789`
  - `recall=0.3667`
  - `f1=0.4490`
  - `total_called=19`
  - `tp_called=11`
  - interpretation: narrow-mode wiring works end-to-end, but recall is too low
    and tuning is the next local loop
- Current local conclusion:
  the replacement-caller path is now structurally integrated and validated in
  both broad and narrow smoke runs; the next bottleneck is caller tuning plus
  broader category-audit cleanup, not basic workflow wiring.

## 2026-05-26 Category Audit And Guardrail Pass

- Completed the explicit six-category audit against the active
  `balanced_tfclean_*_integrated_288` production configs.
- Audit result:
  - `balanced_tfclean_flatearth_peaks_broad_integrated_288` is peak-like
    (`tf_sigma=5`) and now uses `macs2_mode: narrow`
  - both plateau-like families use `tf_sigma=20` and `macs2_mode: broad`
  - both realistic peak-like families use `tf_sigma=5` and
    `macs2_mode: narrow`
  - the only remaining confound was not the config value anymore but the stale
    legacy naming/docs that still called the flatearth peak category `broad`
- Canonical naming cleanup applied:
  - `configs/category_name_map.yaml` now maps the legacy flatearth peak config
    to `flatearth_peak_narrow`
  - `docs/nomenclature_table.md` now records that the filename is legacy but
    the scientific label and decoding mode are narrow
- Documentation cleanup applied:
  - `docs/category_semantics.md` now states the decode contract explicitly:
    peak-like shapes must use narrow decoding and plateau-like shapes must use
    broad decoding
  - removed the stale limitation line claiming MACS2 broad/narrow mode is still
    confounded with the current category labels
- Added a code-level guardrail so this mismatch cannot silently recur:
  - `scripts/eval_helpers.py` now exposes `expected_decode_mode()`,
    `validate_decode_modes()`, and a centralized `truth_mode_for_row()`
  - `Snakefile.py` now validates every generated run row immediately after
    sample construction
  - `peak_pr_stats.py` now derives truth mode from the planted shape helper
    instead of trusting the legacy `macs2_mode` label alone
- Validation proof:
  - `python -m py_compile chipseq_pipeline_v2/scripts/*.py` passed
  - `cd chipseq_pipeline_v2 && HOME=/private/tmp/chipseq_snakemake_home pytest tests/test_eval_helpers.py`
    passed: `6 passed`
  - `cd chipseq_pipeline_v2 && HOME=/private/tmp/chipseq_snakemake_home snakemake -s Snakefile.py --configfile configs/pilot_tfclean_homer_realistic_peaks_wavy_narrow_integrated_1.yaml --dry-run`
    passed with the new validation active
- Current conclusion:
  the six-category semantics are now internally aligned across config,
  canonical naming, documentation, and evaluation code; the next local step is
  HOMER parameter tuning and staged batch growth.

## 2026-05-26 HOMER Narrow Tuning Loop

- Goal of this loop: improve the weak narrow-mode HOMER pilot without touching
  the simulator or the broad-mode path.
- Baseline narrow pilot remained:
  - result root: `results_homer_narrow`
  - metrics:
    - `precision=0.5789`
    - `recall=0.3667`
    - `f1=0.4490`
    - `total_called=19`
- Hypothesis 1: the default HOMER factor settings were too conservative, so
  lower filtering thresholds might recover more planted peaks.
- Added and ran two one-run tuning configs:
  - `configs/tune_tfclean_homer_realistic_peaks_wavy_narrow_relaxed_1.yaml`
    with `-F 2 -L 2 -C 0 -fdr 0.01 -tagThreshold 10`
  - `configs/tune_tfclean_homer_realistic_peaks_wavy_narrow_lenient_1.yaml`
    with `-F 2 -L 1.5 -C 0 -fdr 0.05 -tagThreshold 8`
- Validation proof:
  both tuning configs passed local dry-run and completed actual Snakemake runs
  successfully in `background_project`.
- Both relaxed-filter variants behaved almost identically:
  - `total_called=28`
  - `tp_called=13`
  - `precision=0.4643`
  - `recall=0.4333`
  - `f1=0.4483`
- Interpretation:
  the relaxed settings increased recall modestly (`0.3667 -> 0.4333`) but paid
  for it with a larger precision drop (`0.5789 -> 0.4643`), leaving F1
  essentially unchanged and slightly worse than baseline.
- Hypothesis 2: the main failure might be peak merging, so reduce HOMER peak
  size and minimum peak distance rather than only loosening significance
  thresholds.
- Added and ran a split-oriented config:
  - `configs/tune_tfclean_homer_realistic_peaks_wavy_narrow_split_1.yaml`
    with `-size 100 -minDist 120 -F 3 -L 3 -C 0 -fdr 0.01 -tagThreshold 8`
- Split-oriented result:
  - `total_called=28`
  - `tp_called=9`
  - `precision=0.3214`
  - `recall=0.3000`
  - `f1=0.3103`
- Interpretation:
  smaller HOMER factor windows did not solve the narrow-mode issue; they held
  the total call count at `28` but degraded overlap substantially, so this
  branch should not be used for staged scaling.
- Practical local conclusion from this first tuning pass:
  - broad-mode HOMER remains the stronger initial fit
  - narrow-mode HOMER works operationally but is still underperforming
  - simple threshold relaxation does not improve F1 enough to justify scaling
  - the next step should be either a more targeted narrow-mode HOMER strategy
    or a staged batch run only after deciding whether HOMER is still the right
    narrow-mode replacement path

## 2026-05-26 First Staged HOMER Batch Setup

- Re-checked replacement-caller viability against primary-source docs after the
  tuning loop:
  - MACS3 official docs confirm explicit narrow and broad calling modes
  - HOMER official docs confirm `findPeaks -style factor` and
    `findPeaks -style histone`
  - no better already-integrated broad+narrow candidate was found in the repo,
    so HOMER remains the active replacement path for the next staged local
    validation batch
- Decision:
  proceed with a first staged six-category HOMER batch instead of doing more
  one-run narrow tuning in isolation. The purpose is to inspect multi-category
  plot behavior before investing in the full `1728`-run replacement sweep.
- Added six canonical-name 128-run HOMER configs:
  - `configs/homer_tfclean_flatearth_peak_narrow_integrated_128.yaml`
  - `configs/homer_tfclean_flatearth_plateau_broad_integrated_128.yaml`
  - `configs/homer_tfclean_realistic_peaks_wavy_narrow_integrated_128.yaml`
  - `configs/homer_tfclean_realistic_peaks_hilly_narrow_integrated_128.yaml`
  - `configs/homer_tfclean_realistic_plateaus_wavy_broad_integrated_128.yaml`
  - `configs/homer_tfclean_realistic_plateaus_hilly_broad_integrated_128.yaml`
- The staged 128-run design is:
  - `coverage_treat = [5, 10]`
  - `coverage_ctrl = [0.5, 1, 2, 4, 8, 12, 16, 24]`
  - `tf_enrich = 2 values`
  - `seed = 4 values`
  - total per category = `2 x 8 x 2 x 4 = 128` runs
- Validation proof:
  all six configs dry-run cleanly and each resolves to exactly `128` runs.
- To get real runtime and anomaly signal before launching all six categories at
  once, started two representative local staged runs in `background_project`:
  - narrow representative:
    `configs/homer_tfclean_realistic_peaks_wavy_narrow_integrated_128.yaml`
  - broad representative:
    `configs/homer_tfclean_realistic_plateaus_wavy_broad_integrated_128.yaml`
- Runtime proof at launch:
  each representative staged run expands to `770` Snakemake jobs:
  - `128` treatment simulations
  - `128` control simulations
  - `256` Bowtie2 alignments
  - `128` HOMER tag-directory builds
  - `128` HOMER peak calls
  - `1` params-table write
  - `1` `all` target
- Early execution signal:
  both staged runs began successfully and advanced through the initial
  simulation phase without immediate decode-mode, path, or HOMER invocation
  errors.
- Current status:
  representative narrow and broad staged runs are in progress locally; next
  checkpoint is to capture completion, compute per-category performance tables,
  and inspect the first batch plots for low-precision or low-recall pathologies.

## 2026-05-26 Representative Batch Runtime Checkpoint

- Rechecked the representative 128-run staged sessions directly instead of
  inferring state from a stale file count.
- Important finding:
  the staged runs were not failed or empty; they were still active inside the
  existing exec sessions, and the previously inspected log files were simply
  partial snapshots.
- Narrow representative status at checkpoint:
  - config:
    `configs/homer_tfclean_realistic_peaks_wavy_narrow_integrated_128.yaml`
  - live session: `38941`
  - progress observed from session output: `172 / 770 steps (22%) done`
- Broad representative status at checkpoint:
  - config:
    `configs/homer_tfclean_realistic_plateaus_wavy_broad_integrated_128.yaml`
  - live session: `15911`
  - progress observed from session output: `177 / 770 steps (23%) done`
- Intermediate output counts at the same checkpoint still show only simulation
  artifacts and no alignments or called peaks yet:
  - narrow root: `75` treatment `pmf.disabled`, `79` control `pmf.disabled`,
    `0` BAMs, `0` HOMER peak BEDs
  - broad root: `71` treatment `pmf.disabled`, `83` control `pmf.disabled`,
    `0` BAMs, `0` HOMER peak BEDs
- Interpretation:
  these representative staged runs are still in the simulator-heavy front half,
  so no PR/recall/F1 evaluation can happen yet. The next meaningful checkpoint
  is either:
  - first appearance of aligned BAM outputs, or
  - completion of the staged categories

## 2026-05-26 Representative Batch Runtime Checkpoint 2

- Polled the live representative sessions again and confirmed both are still
  actively progressing.
- Narrow representative updated status:
  - live session: `38941`
  - latest observed progress: `229 / 770 steps (30%) done`
  - current artifact counts:
    - `119` treatment `pmf.disabled`
    - `116` control `pmf.disabled`
    - `0` aligned BAMs
    - `0` HOMER peak BEDs
- Broad representative updated status:
  - live session: `15911`
  - latest observed progress: `231 / 770 steps (30%) done`
  - current artifact counts:
    - `116` treatment `pmf.disabled`
    - `115` control `pmf.disabled`
    - `0` aligned BAMs
    - `0` HOMER peak BEDs
- Interpretation:
  the representative staged runs are still consuming the simulation phase in a
  healthy way, but they have not yet crossed the threshold where BAMs or peak
  calls exist. Evaluation, plot generation, and anomaly assessment remain
  blocked on reaching the alignment/calling half of the DAG.

## 2026-05-26 Representative Batch Runtime Checkpoint 3

- Both representative staged runs have now crossed out of pure simulation and
  into the alignment half of the DAG.
- Narrow representative updated status:
  - live session: `38941`
  - latest observed progress: `259 / 770 steps (34%) done`
  - current artifact counts:
    - `128` treatment `pmf.disabled`
    - `128` control `pmf.disabled`
    - `0` treatment BAMs
    - `3` control BAMs
    - `0` HOMER peak BEDs
- Broad representative updated status:
  - live session: `15911`
  - latest observed progress: `259 / 770 steps (34%) done`
  - current artifact counts:
    - `128` treatment `pmf.disabled`
    - `128` control `pmf.disabled`
    - `0` treatment BAMs
    - `3` control BAMs
    - `0` HOMER peak BEDs
- Interpretation:
  the representative batch design is now fully through the simulation stage for
  both categories and has started Bowtie2 alignments. This is the first hard
  proof that the staged HOMER 128-run configs are progressing through the real
  downstream half of the pipeline instead of only the front-end simulator
  stage. Metric generation is still waiting on treatment BAMs and peak BEDs.

## 2026-05-26 Representative Batch Runtime Checkpoint 4

- Polled both representative sessions again after the first alignment
  transition.
- Narrow representative updated status:
  - live session: `38941`
  - latest observed progress: `266 / 770 steps (35%) done`
  - current artifact counts:
    - `128` treatment `pmf.disabled`
    - `128` control `pmf.disabled`
    - `0` treatment BAMs
    - `10` control BAMs
    - `0` HOMER peak BEDs
- Broad representative updated status:
  - live session: `15911`
  - latest observed progress: `266 / 770 steps (35%) done`
  - current artifact counts:
    - `128` treatment `pmf.disabled`
    - `128` control `pmf.disabled`
    - `0` treatment BAMs
    - `9` control BAMs
    - `0` HOMER peak BEDs
- Interpretation:
  the runs are steadily consuming control-side alignments first. The next
  decisive checkpoint is the first treatment BAM, because that unlocks HOMER
  tag-directory generation and peak calling shortly afterward.

## 2026-05-26 Representative Batch Runtime Checkpoint 5

- Polled both representative sessions again specifically for the first
  treatment-side downstream artifacts.
- Narrow representative updated status:
  - live session: `38941`
  - latest observed progress: `273 / 770 steps (35%) done`
  - current downstream artifact counts:
    - `0` treatment BAMs
    - `18` control BAMs
    - `0` HOMER treat tag directories
    - `0` HOMER peak BEDs
- Broad representative updated status:
  - live session: `15911`
  - latest observed progress: `273 / 770 steps (35%) done`
  - current downstream artifact counts:
    - `0` treatment BAMs
    - `16` control BAMs
    - `0` HOMER treat tag directories
    - `0` HOMER peak BEDs
- Interpretation:
  both representative staged runs are still processing the control-alignment
  queue first. The critical downstream unlock remains unchanged: the first
  treatment BAM will be the point where HOMER tag-directory creation and peak
  calling can begin.

## 2026-05-26 Representative Batch Runtime Checkpoint 6

- Continued polling specifically for the first treatment-side downstream
  artifacts.
- Narrow representative updated status:
  - live session: `38941`
  - latest observed progress: `282 / 770 steps (37%) done`
  - current downstream artifact counts:
    - `0` treatment BAMs
    - `27` control BAMs
    - `0` HOMER treat tag directories
    - `0` HOMER peak BEDs
- Broad representative updated status:
  - live session: `15911`
  - latest observed progress: `282 / 770 steps (37%) done`
  - current downstream artifact counts:
    - `0` treatment BAMs
    - `25` control BAMs
    - `0` HOMER treat tag directories
    - `0` HOMER peak BEDs
- Interpretation:
  the downstream phase is still healthy, but the queue ordering continues to
  favor control BAM production first. No HOMER-side evaluation artifacts are
  available yet; metric generation remains gated on the first treatment BAM.

## 2026-05-26 Representative Batch Runtime Checkpoint 7

- Polled both representative sessions again for the first treatment-side
  artifacts and HOMER outputs.
- Narrow representative updated status:
  - live session: `38941`
  - latest observed progress: `293 / 770 steps (38%) done`
  - current downstream artifact counts:
    - `0` treatment BAMs
    - `37` control BAMs
    - `0` HOMER treat tag directories
    - `0` HOMER peak BEDs
- Broad representative updated status:
  - live session: `15911`
  - latest observed progress: `292 / 770 steps (38%) done`
  - current downstream artifact counts:
    - `0` treatment BAMs
    - `35` control BAMs
    - `0` HOMER treat tag directories
    - `0` HOMER peak BEDs
- Interpretation:
  both representative staged runs remain healthy and continue to advance
  through the control-alignment queue. The first treatment BAM is still the
  gating downstream event for HOMER tag-directory creation, peak calling, and
  staged metric generation.

## 2026-05-26 Representative Batch Runtime Checkpoint 8

- Re-polled both representative sessions and hit the first treatment-side
  downstream milestone.
- Narrow representative updated status:
  - live session: `38941`
  - latest observed progress: `304 / 770 steps (39%) done`
  - current downstream artifact counts:
    - `1` treatment BAM
    - `48` control BAMs
    - `0` HOMER treat tag directories
    - `0` HOMER peak BEDs
- Broad representative updated status:
  - live session: `15911`
  - latest observed progress: `304 / 770 steps (39%) done`
  - current downstream artifact counts:
    - `0` treatment BAMs
    - `47` control BAMs
    - `0` HOMER treat tag directories
    - `0` HOMER peak BEDs
- Interpretation:
  the narrow representative has reached the first treatment-alignment milestone,
  which is the earliest downstream proof that HOMER-side tag-directory creation
  can start soon. The next decisive checkpoint is the first HOMER treat
  tag-directory or the first HOMER peak BED.

## 2026-05-26 First Representative Batch Metrics And Plots

- Both representative staged sessions completed end-to-end successfully:
  - narrow session `38941`
  - broad session `15911`
- Completion proof on disk for each representative root:
  - `128` treatment BAMs
  - `128` control BAMs
  - `128` HOMER treat tag directories
  - `128` HOMER control tag directories
  - `128` HOMER peak BEDs
- Computed staged metrics with `scripts/peak_pr_stats.py` for:
  - `analysis_outputs/homer_tfclean_realistic_peaks_wavy_narrow_128_stats`
  - `analysis_outputs/homer_tfclean_realistic_plateaus_wavy_broad_128_stats`
- Built the first staged representative plot/report assets:
  - `analysis_outputs/homer_tfclean_realistic_peaks_wavy_narrow_128_stats/pr_recall_f1_vs_ctrl_coverage.png`
  - `analysis_outputs/homer_tfclean_realistic_plateaus_wavy_broad_128_stats/pr_recall_f1_vs_ctrl_coverage.png`
  - `analysis_outputs/homer_representative_128_report/README.md`
- Narrow representative headline behavior:
  - best F1: `0.52` at `coverage_ctrl=12`
  - worst recall in summary table: `0.20`
  - precision stays moderate (`~0.55` range), but recall is the limiting factor
- Broad representative headline behavior:
  - best F1: `0.6584` at `coverage_ctrl=1`
  - minimum precision across the shown summary points: `0.9144`
  - minimum recall across the shown summary points: `0.3208`
  - interpretation: broad HOMER behavior is much stronger and cleaner than the
    narrow HOMER path
- Example summary rows:
  - narrow at `coverage_ctrl=1`: precision `0.5588`, recall `0.4354`,
    F1 `0.4895`
  - narrow at `coverage_ctrl=12`: precision `0.5525`, recall `0.4917`,
    F1 `0.5203`
  - broad at `coverage_ctrl=1`: precision `0.9484`, recall `0.5042`,
    F1 `0.6584`
  - broad at `coverage_ctrl=4`: precision `0.9776`, recall `0.4583`,
    F1 `0.6241`
- First staged anomaly conclusion:
  the concern from pilot runs persists in the larger representative batch:
  HOMER is a plausible replacement on broad-like categories, but narrow-like
  categories remain substantially weaker because recall lags while precision
  stays only middling instead of excellent.
- Decision signal for the next phase:
  do not blindly scale HOMER to all six categories yet as if the replacement
  problem is solved. The representative batch evidence says we either need:
  - a more targeted narrow-mode HOMER strategy, or
  - a different replacement-caller strategy for the narrow side before treating
    the whole 1700-run replacement sweep as validated.

## 2026-05-26 Remaining Four Staged HOMER Categories Launched

- To satisfy the actual staged-validation goal shape, launched the remaining
  four `128`-run HOMER categories locally instead of stopping at the two
  representative categories:
  - `configs/homer_tfclean_flatearth_peak_narrow_integrated_128.yaml`
  - `configs/homer_tfclean_flatearth_plateau_broad_integrated_128.yaml`
  - `configs/homer_tfclean_realistic_peaks_hilly_narrow_integrated_128.yaml`
  - `configs/homer_tfclean_realistic_plateaus_hilly_broad_integrated_128.yaml`
- Live sessions:
  - flatearth peak narrow: `20626`
  - flatearth plateau broad: `81638`
  - hilly peak narrow: `23790`
  - hilly plateau broad: `91775`
- Initial validation proof:
  each launched config expanded to the same staged shape as the representative
  runs:
  - `770` Snakemake jobs total
  - `128` treatment simulations
  - `128` control simulations
  - `256` Bowtie2 alignments
  - `128` HOMER tag-directory jobs
  - `128` HOMER peak calls
  - `1` params-table write
  - `1` `all` target
- Early execution status:
  all four newly launched categories started successfully and entered the
  simulation phase without immediate decode-mode, path-resolution, or HOMER
  rule-shape failures.
- Current practical meaning:
  the staged validation wave is now truly moving toward a six-category HOMER
  report rather than stopping at a two-category probe. The next durable
  milestone is to harvest these four completions and build the first full
  six-category staged plot set.

## 2026-05-26 Early Runtime Checkpoint for Remaining Four HOMER Categories

- Polled the four newly launched staged HOMER sessions to confirm they are
  progressing past launch and into real work:
  - flatearth peak narrow: session `20626`
  - flatearth plateau broad: session `81638`
  - hilly peak narrow: session `23790`
  - hilly plateau broad: session `91775`
- Two sessions emitted clear buffered runtime progress:
  - `20626`: `97 / 770` jobs complete (`13%`)
  - `81638`: `99 / 770` jobs complete (`13%`)
- The other two sessions remained active but did not emit useful buffered text
  in that poll window:
  - `23790`
  - `91775`
- On-disk artifact counts for all four roots were still zero at this checkpoint
  for:
  - treatment BAMs
  - control BAMs
  - HOMER peak BEDs
- Interpretation:
  this is still an early simulation-phase checkpoint, not a scientific
  evaluation checkpoint. There is no anomaly evidence yet because the runs have
  not reached alignment or peak-calling outputs.

## 2026-05-26 Refined Runtime Checkpoint for Remaining Four HOMER Categories

- Follow-up polling clarified that the hilly pair were not stuck; they had
  simply started later and were still too early to emit on-disk artifacts in
  the previous poll window.
- Refined session status:
  - `20626` flatearth peak narrow: active and already through about `13%`
  - `81638` flatearth plateau broad: active and already through about `13%`
  - `23790` realistic peaks hilly narrow: active, early simulation phase
    (`4 / 770`, about `1%`)
  - `91775` realistic plateaus hilly broad: active, early simulation phase
    (`4 / 770`, about `1%`)
- Artifact spot-check on the two earlier-launched roots now shows real
  simulation outputs accumulating:
  - `results_homer_tfclean_flatearth_peak_narrow_128`: `57` treatment FASTA
    pairs, `61` control FASTA pairs
  - `results_homer_tfclean_flatearth_plateau_broad_128`: `56` treatment FASTA
    pairs, `63` control FASTA pairs
- Updated interpretation:
  all four remaining staged categories are alive and progressing. The first two
  are materially ahead, but none currently show evidence of decode-mode,
  HOMER-rule, or path-resolution failure.

## 2026-05-26 Runtime Checkpoint: First Two Remaining Categories Reach ~19%

- Re-polled the two earlier-launched remaining categories and confirmed they
  continue to advance through the simulation-heavy front half:
  - `20626` flatearth peak narrow: `147 / 770` jobs complete (`19%`)
  - `81638` flatearth plateau broad: `146 / 770` jobs complete (`19%`)
- Re-polled the two later-launched hilly categories; they remained active but
  did not emit additional buffered text in that poll window:
  - `23790`
  - `91775`
- On-disk scoring-readiness spot-check still showed zero for all four roots at
  this checkpoint for:
  - treatment BAMs
  - control BAMs
  - HOMER peak BEDs
- Interpretation:
  this staged wave is still not at the first scientific scoring checkpoint.
  The workflow remains healthy, but the current evidence is still runtime
  progress rather than evaluation output.

## 2026-05-26 Runtime Checkpoint: First Two Remaining Categories Reach ~22%

- Another poll confirms the two earlier-launched remaining categories continue
  to advance cleanly through the simulation-heavy front half:
  - `20626` flatearth peak narrow: `169 / 770` jobs complete (`22%`)
  - `81638` flatearth plateau broad: `168 / 770` jobs complete (`22%`)
- The two hilly sessions remained active but again did not emit useful buffered
  text in this poll window:
  - `23790`
  - `91775`
- Scoring-readiness spot-check still shows zero across all four remaining roots
  for:
  - treatment BAMs
  - control BAMs
  - HOMER peak BEDs
- Interpretation:
  the staged validation wave is still healthy but still pre-alignment. We are
  not yet at the first point where scientific performance metrics or six-
  category plots can be computed from these four categories.

## 2026-05-26 Runtime Checkpoint: First Two Remaining Categories Reach ~25%

- Another poll confirms continued healthy forward progress for the two
  earlier-launched remaining categories:
  - `20626` flatearth peak narrow: `190 / 770` jobs complete (`25%`)
  - `81638` flatearth plateau broad: `190 / 770` jobs complete (`25%`)
- The two later-launched hilly categories again remained active without
  producing additional buffered text in this poll window:
  - `23790`
  - `91775`
- Scoring-readiness spot-check is still zero across all four remaining roots
  for:
  - treatment BAMs
  - control BAMs
  - HOMER peak BEDs
- Interpretation:
  the staged wave remains healthy and is still in the simulation-heavy front
  half. The first six-category scoring milestone still has not been reached.

## 2026-05-26 MACS2 Six-Category Rerun Path Audit

- Used the current guarded workflow and active balanced TF-clean configs as the
  authoritative rerun path for the planned MACS2 retry wave.
- Verified the active six-category TF-clean production configs are decoder-safe
  on the shape/mode contract:
  - `balanced_tfclean_flatearth_peaks_broad_integrated_288.yaml`:
    `tf_sigma: [5]`, `macs2_mode: ["narrow"]`
  - `balanced_tfclean_flatearth_plateaus_broad_integrated_288.yaml`:
    `tf_sigma: [20]`, `macs2_mode: ["broad"]`
  - `balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288.yaml`:
    peak-like / narrow
  - `balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288.yaml`:
    peak-like / narrow
  - `balanced_tfclean_realistic_plateaus_wavy_broad_integrated_288.yaml`:
    plateau-like / broad
  - `balanced_tfclean_realistic_plateaus_hilly_broad_integrated_288.yaml`:
    plateau-like / broad
- Verified the workflow-level guardrail remains active:
  - `Snakefile.py` calls `validate_decode_modes(SAMPLES)` immediately after
    sample construction
  - `scripts/eval_helpers.py` defines the invariant:
    - peak-like / narrow for `tf_sigma <= 5`
    - plateau-like / broad for `tf_sigma >= 15`
    - intermediate values raise rather than silently mapping into headline
      categories
- Verified the canonical labeling/docs layer also now matches the guarded
  interpretation:
  - `configs/category_name_map.yaml`
  - `docs/category_semantics.md`
  - `docs/nomenclature_table.md`
- Added an explicit decoder-safety note to `chipseq_pipeline_v2/README.md`
  explaining that the flatearth peak config retains a stale filename token
  containing `broad`, but its planted shape and decode mode are narrow and the
  canonical label is `flatearth_peak_narrow`.
- Validation:
  - `python -m py_compile chipseq_pipeline_v2/scripts/*.py` passed after the
    documentation cleanup
- Conclusion:
  the active MACS2 six-category rerun path is now aligned with the planted
  shape contract and documented clearly enough that the old broad/narrow
  confound should not silently recur on the guarded path.

## 2026-05-26 First Downstream Handoff in Remaining HOMER Staged Batch

- The two earlier-launched remaining HOMER categories have now crossed out of
  pure simulation and into the downstream alignment phase:
  - `20626` flatearth peak narrow: `258 / 770` jobs complete (`34%`)
  - `81638` flatearth plateau broad: `257 / 770` jobs complete (`33%`)
- First alignment-side proof:
  - flatearth peak narrow:
    - `128` treatment FASTA pairs present
    - `126` control FASTA pairs present
    - `2` control BAMs present
    - `0` treatment BAMs present
    - `0` HOMER peak BEDs present
  - flatearth plateau broad:
    - `128` treatment FASTA pairs present
    - `127` control FASTA pairs present
    - `1` control BAM present
    - `0` treatment BAMs present
    - `0` HOMER peak BEDs present
- Session output confirms the first downstream gate directly:
  - `align_bowtie2` started on run `0096` control for both flatearth staged
    roots
  - one completed example showed a clean Bowtie2 alignment with `100% overall
    alignment rate` on the control sample before BAM sort/index completion
- The later-launched hilly staged categories are still much earlier:
  - `23790` realistic peaks hilly narrow: `8 / 770` jobs complete (`1%`)
  - `91775` realistic plateaus hilly broad: `8 / 770` jobs complete (`1%`)
  - each currently shows `4` treatment FASTA pairs and `4` control FASTA pairs
    on disk, but still no BAMs or HOMER peak BEDs
- Interpretation:
  the remaining staged six-category HOMER wave is now genuinely in the
  downstream half for the two flatearth categories, but still not yet at the
  first scoring point because no HOMER peak BEDs exist yet and treatment BAM
  generation has not started.

## 2026-05-26 Downstream Checkpoint: Control-Side Alignments Accumulating

- Follow-up polling confirms the two flatearth staged categories are now
  accumulating aligned control BAMs rather than just crossing the first
  alignment boundary:
  - `20626` flatearth peak narrow: `261 / 770` jobs complete (`34%`)
  - `81638` flatearth plateau broad: `261 / 770` jobs complete (`34%`)
- Updated artifact counts:
  - flatearth peak narrow:
    - `6` control BAMs
    - `0` treatment BAMs
    - `0` HOMER peak BEDs
  - flatearth plateau broad:
    - `5` control BAMs
    - `0` treatment BAMs
    - `0` HOMER peak BEDs
- Session output again showed clean control-side Bowtie2 completions for both
  flatearth roots, including additional completed runs such as:
  - peak-narrow control runs `0095` and `0080`
  - plateau-broad control runs `0095` and `0080`
- The hilly staged categories are still early simulation-phase jobs:
  - `23790` realistic peaks hilly narrow: `11 / 770` jobs complete (`1%`)
  - `91775` realistic plateaus hilly broad: `12 / 770` jobs complete (`2%`)
- Interpretation:
  the flatearth pair are now firmly in downstream control-side alignment, but
  the first scoring gate still has not opened because treatment BAMs and HOMER
  peak BEDs remain absent.

## 2026-05-26 Downstream Checkpoint: Control BAMs Reach Double Digits

- Another downstream poll confirms both flatearth staged categories are still
  progressing cleanly through control-side alignments:
  - `20626` flatearth peak narrow: `266 / 770` jobs complete (`35%`)
  - `81638` flatearth plateau broad: `266 / 770` jobs complete (`35%`)
- Updated downstream artifact counts:
  - flatearth peak narrow:
    - `10` control BAMs
    - `0` treatment BAMs
    - `0` HOMER treat tag directories
    - `0` HOMER peak BEDs
  - flatearth plateau broad:
    - `10` control BAMs
    - `0` treatment BAMs
    - `0` HOMER treat tag directories
    - `0` HOMER peak BEDs
- The hilly categories remain early:
  - `23790` realistic peaks hilly narrow: `12 / 770` jobs complete (`2%`)
  - `91775` realistic plateaus hilly broad: still active but too early to show
    downstream artifacts
- Interpretation:
  the flatearth pair continue to move deeper into the downstream half, but the
  treatment-side unlock has still not started, so HOMER peak-calling and staged
  scoring remain gated.

## 2026-05-26 Downstream Checkpoint: Control BAMs Continue Rising

- Another poll confirms the flatearth staged pair are still progressing cleanly
  through control-side alignment and have not stalled:
  - `20626` flatearth peak narrow: `271 / 770` jobs complete (`35%`)
  - `81638` flatearth plateau broad: `270 / 770` jobs complete (`35%`)
- Updated downstream artifact counts:
  - flatearth peak narrow:
    - `15` control BAMs
    - `0` treatment BAMs
    - `0` HOMER treat tag directories
    - `0` HOMER peak BEDs
  - flatearth plateau broad:
    - `14` control BAMs
    - `0` treatment BAMs
    - `0` HOMER treat tag directories
    - `0` HOMER peak BEDs
- The hilly categories are still too early to contribute downstream artifacts:
  - `23790` realistic peaks hilly narrow: still active, no BAMs yet
  - `91775` realistic plateaus hilly broad: still active, no BAMs yet
- Interpretation:
  the staged wave is healthy, but the same gating fact still holds: HOMER peak
  calling and therefore staged scoring remain blocked on the first treatment
  BAMs, which have not appeared yet.

## 2026-05-26 Downstream Checkpoint: Control BAMs Reach ~20

- Another downstream gate check confirms the flatearth pair continue to make
  healthy progress, but still exclusively on the control side:
  - `20626` flatearth peak narrow: `275 / 770` jobs complete (`36%`)
  - `81638` flatearth plateau broad: `275 / 770` jobs complete (`36%`)
- Updated downstream artifact counts:
  - flatearth peak narrow:
    - `20` control BAMs
    - `0` treatment BAMs
    - `0` HOMER treat tag directories
    - `0` HOMER peak BEDs
  - flatearth plateau broad:
    - `19` control BAMs
    - `0` treatment BAMs
    - `0` HOMER treat tag directories
    - `0` HOMER peak BEDs
- The hilly staged categories remain early simulation jobs:
  - `23790` realistic peaks hilly narrow: `15 / 770` jobs complete (`2%`)
  - `91775` realistic plateaus hilly broad: `16 / 770` jobs complete (`2%`)
- Interpretation:
  the batch remains healthy, but the core gating fact is unchanged: treatment
  BAM generation has still not started, so HOMER peak calling and the first
  full six-category staged scoring pass are still blocked on runtime progress.

## 2026-05-26 Control-First Scheduling Audit

- Inspected the active workflow rules to understand why the staged HOMER batch
  is moving through many control-side alignments before any treatment BAMs
  appear.
- Key workflow facts:
  - `simulation.smk` emits both control and treatment FASTA pairs independently
  - `alignment.smk` defines the same `align_bowtie2` rule shape for both
    conditions with no condition-specific priority
  - `peakcalling.smk` requires both treatment and control BAMs before HOMER can
    build tag directories or call peaks
- Strongest runtime evidence from the two flatearth staged roots:
  - flatearth peak narrow:
    - `128` treatment FASTA pairs currently present
    - `104` control FASTA pairs currently present
    - `24` control BAMs already produced
    - `0` treatment BAMs produced so far
  - flatearth plateau broad:
    - `128` treatment FASTA pairs currently present
    - `104` control FASTA pairs currently present
    - `24` control BAMs already produced
    - `0` treatment BAMs produced so far
- Interpretation:
  this pattern points to Snakemake scheduling/consumption order rather than a
  broken treatment path. Control FASTAs are being consumed and deleted by the
  align rule first, while treatment FASTAs are already present and waiting for
  their turn in the same align queue.
- Practical conclusion:
  the current state does not justify a workflow fix yet. The staged batch looks
  healthy, and the right next move is still to keep harvesting until the first
  treatment BAMs and then HOMER outputs appear.

## 2026-05-26 Downstream Checkpoint: Control BAMs Reach ~30

- Another gate-focused poll confirms the treatment-side unlock still has not
  started, but the flatearth staged pair continue to advance steadily through
  control-side alignments:
  - `20626` flatearth peak narrow: `286 / 770` jobs complete (`37%`)
  - `81638` flatearth plateau broad: `286 / 770` jobs complete (`37%`)
- Updated downstream artifact counts:
  - flatearth peak narrow:
    - `31` control BAMs
    - `0` treatment BAMs
    - `0` HOMER tag directories
    - `0` HOMER peak BEDs
  - flatearth plateau broad:
    - `30` control BAMs
    - `0` treatment BAMs
    - `0` HOMER tag directories
    - `0` HOMER peak BEDs
- Hilly staged categories remain early simulation jobs:
  - `23790` realistic peaks hilly narrow: `17 / 770` jobs complete (`2%`)
  - `91775` realistic plateaus hilly broad: still active, no downstream
    artifacts yet
- Interpretation:
  the batch remains healthy and the scheduling interpretation still holds, but
  the first full six-category scoring pass is still blocked on treatment BAM
  generation and then HOMER peak calling.

## 2026-05-26 Six-Category HOMER Report Builder Audit And Patch

- Audited the report-building path that will consume the remaining staged HOMER
  stats once all six category directories are available.
- Found a real labeling gap:
  `scripts/build_balanced_288_config_report.py` used each input directory's raw
  name as the category label, which would have carried stale names like
  `balanced_tfclean_flatearth_peaks_broad_integrated_288` into the staged HOMER
  report instead of the cleaned canonical label `flatearth_peak_narrow`.
- Patched the report builder to:
  - accept `--category-map` with default
    `configs/category_name_map.yaml`
  - load canonical labels via `scripts.category_summary_lib`
  - use canonical labels in plot titles and `data_info.md`
  - write the root report README with canonical names while still preserving
    the raw input directory name for traceability
  - generalize the root README title from `Balanced 288 Config Report` to
    `Balanced Config Report`
- Validation:
  - `python -m py_compile chipseq_pipeline_v2/scripts/build_balanced_288_config_report.py chipseq_pipeline_v2/scripts/category_summary_lib.py` passed
- Practical meaning:
  once the remaining HOMER stats exist, the first full staged six-category
  report can be built without reintroducing the stale flatearth broad/narrow
  naming confusion in its figures or summaries.

## 2026-05-26 Downstream Checkpoint: Control BAMs Reach ~40

- Latest gate-focused runtime snapshot:
  - `20626` flatearth peak narrow: still active, control-side alignment queue
    continuing
  - `81638` flatearth plateau broad: still active, control-side alignment queue
    continuing
- Updated downstream artifact counts:
  - flatearth peak narrow:
    - `42` control BAMs
    - `0` treatment BAMs
    - `0` HOMER peak BEDs
  - flatearth plateau broad:
    - `41` control BAMs
    - `0` treatment BAMs
    - `0` HOMER peak BEDs
  - hilly pair:
    - still no BAMs or HOMER peak BEDs
- Interpretation:
  the live staged batch remains healthy, the report path is now canonical-label
  ready, and the remaining blocker is still purely runtime: treatment BAMs must
  appear before HOMER calling and the first full six-category staged scoring
  pass can start.

## 2026-05-26 Treatment-Side Unlock Reached For Flatearth HOMER Categories

- The first treatment-side BAMs have now appeared in the remaining staged HOMER
  batch, which is the key downstream gate we had been waiting for.
- Current flatearth staged artifact counts:
  - flatearth peak narrow:
    - `48` control BAMs
    - `4` treatment BAMs
    - `0` HOMER tag directories
    - `0` HOMER peak BEDs
  - flatearth plateau broad:
    - `48` control BAMs
    - `2` treatment BAMs
    - `0` HOMER tag directories
    - `0` HOMER peak BEDs
- Session output confirms the queue has started launching treatment alignments,
  for example:
  - flatearth peak narrow:
    - completed treat alignment on run `0054`
    - launched treat alignment on run `0052`
  - flatearth plateau broad:
    - completed treat alignment on run `0054`
    - launched treat alignment on run `0052`
- Current session progress:
  - `20626` flatearth peak narrow: `306 / 770` jobs complete (`40%`)
  - `81638` flatearth plateau broad: `306 / 770` jobs complete (`40%`)
  - `23790` realistic peaks hilly narrow: `21 / 770` jobs complete (`3%`)
  - `91775` realistic plateaus hilly broad: `21 / 770` jobs complete (`3%`)
- Interpretation:
  the runtime blocker has materially changed. HOMER tag-directory creation is
  now on the critical path for the flatearth pair, so the next meaningful
  checkpoint is the first HOMER tagdir and then the first HOMER peak BED.

## 2026-05-26 Treatment BAMs Expanding, HOMER Still Pending

- Follow-up gate check shows the flatearth pair are now clearly progressing
  through multiple treatment-side alignments rather than only the first few:
  - flatearth peak narrow:
    - `48` control BAMs
    - `13` treatment BAMs
    - `0` HOMER tag directories
    - `0` HOMER peak BEDs
  - flatearth plateau broad:
    - `48` control BAMs
    - `12` treatment BAMs
    - `0` HOMER tag directories
    - `0` HOMER peak BEDs
- Session evidence confirms continued treatment-side alignment launches and
  completions, for example:
  - flatearth peak narrow:
    - completed treat alignments through runs such as `0052`, `0054`, `0128`
    - launched additional treat alignments such as `0050`, `0126`
  - flatearth plateau broad:
    - completed treat alignments through runs such as `0051`, `0052`, `0054`
    - launched additional treat alignments such as `0050`, `0063`
- Current session progress:
  - `20626` flatearth peak narrow: `315 / 770` jobs complete (`41%`)
  - `81638` flatearth plateau broad: `315 / 770` jobs complete (`41%`)
- Interpretation:
  the next gate is now very specific: once both BAM sides are sufficiently
  available for a run, HOMER tag-directory creation should begin. The first
  HOMER tagdir is now the most important next checkpoint.

## 2026-05-26 Treatment BAMs Reach ~20, HOMER Still Not Started

- Another gate-focused poll shows the flatearth pair have moved materially
  further through treatment-side alignment, but HOMER tag-directory creation
  has still not started yet.
- Updated flatearth staged counts:
  - flatearth peak narrow:
    - `48` control BAMs
    - `21` treatment BAMs
    - `0` HOMER tag directories
    - `0` HOMER peak BEDs
  - flatearth plateau broad:
    - `48` control BAMs
    - `20` treatment BAMs
    - `0` HOMER tag directories
    - `0` HOMER peak BEDs
- Current progress snapshots:
  - `20626` flatearth peak narrow: `324 / 770` jobs complete (`42%`)
  - `81638` flatearth plateau broad: `323 / 770` jobs complete (`42%`)
  - `23790` realistic peaks hilly narrow: `25 / 770` jobs complete (`3%`)
  - `91775` realistic plateaus hilly broad: `25 / 770` jobs complete (`3%`)
- Interpretation:
  the batch remains healthy and the next meaningful milestone is still the same:
  the first HOMER tagdir, followed by the first HOMER peak BED and then staged
  scoring.

## 2026-05-26 Treatment BAMs Reach ~30, HOMER Still Not Started

- Another gate-focused poll shows the flatearth pair continue to advance
  through treatment-side alignment, but HOMER tag-directory creation still has
  not started.
- Updated flatearth staged counts:
  - flatearth peak narrow:
    - `48` control BAMs
    - `30` treatment BAMs
    - `0` HOMER tag directories
    - `0` HOMER peak BEDs
  - flatearth plateau broad:
    - `48` control BAMs
    - `29` treatment BAMs
    - `0` HOMER tag directories
    - `0` HOMER peak BEDs
- Current progress snapshots:
  - `20626` flatearth peak narrow: `333 / 770` jobs complete (`43%`)
  - `81638` flatearth plateau broad: `332 / 770` jobs complete (`43%`)
- Interpretation:
  the staged batch remains healthy and is moving deeper into the treatment-side
  align queue, but the first HOMER tagdir is still the next real milestone
  before any staged scoring can begin.

## 2026-05-26 Treatment BAMs Reach ~40, HOMER Still Not Started

- Another gate check shows the flatearth pair are now well into treatment-side
  alignment, but HOMER tag-directory creation still has not started.
- Updated flatearth staged counts:
  - flatearth peak narrow:
    - `48` control BAMs
    - `40` treatment BAMs
    - `0` HOMER tag directories
    - `0` HOMER peak BEDs
  - flatearth plateau broad:
    - `48` control BAMs
    - `39` treatment BAMs
    - `0` HOMER tag directories
    - `0` HOMER peak BEDs
- Current progress snapshots:
  - `20626` flatearth peak narrow: `343 / 770` jobs complete (`45%`)
  - `81638` flatearth plateau broad: `342 / 770` jobs complete (`44%`)
- Interpretation:
  the batch is still healthy, but the first HOMER tagdir remains the true next
  milestone. Once that appears, peak calling and staged scoring should finally
  be able to begin.

## 2026-05-26 Treatment BAMs Reach ~50, HOMER Still Not Started

- Latest gate-focused poll:
  - flatearth peak narrow:
    - `48` control BAMs
    - `52` treatment BAMs
    - `0` HOMER tag directories
    - `0` HOMER peak BEDs
  - flatearth plateau broad:
    - `48` control BAMs
    - `50` treatment BAMs
    - `0` HOMER tag directories
    - `0` HOMER peak BEDs
- Current progress snapshots:
  - `20626` flatearth peak narrow: `354 / 770` jobs complete (`46%`)
  - `81638` flatearth plateau broad: `353 / 770` jobs complete (`46%`)
  - `23790` realistic peaks hilly narrow: `29 / 770` jobs complete (`4%`)
  - `91775` realistic plateaus hilly broad: `29 / 770` jobs complete (`4%`)
- Interpretation:
  the staged batch remains healthy and the flatearth pair are now deep into
  treatment-side alignment, but the first HOMER tagdir is still the next real
  milestone before staged scoring can begin.

## 2026-05-26 Treatment BAMs Reach ~60, HOMER Still Not Started

- Another gate-focused poll shows the flatearth pair continuing to progress
  through treatment-side alignment with no HOMER tag-directory creation yet.
- Updated flatearth staged counts:
  - flatearth peak narrow:
    - `48` control BAMs
    - `63` treatment BAMs
    - `0` HOMER tag directories
    - `0` HOMER peak BEDs
  - flatearth plateau broad:
    - `48` control BAMs
    - `61` treatment BAMs
    - `0` HOMER tag directories
    - `0` HOMER peak BEDs
- Current progress snapshots:
  - `20626` flatearth peak narrow: `365 / 770` jobs complete (`47%`)
  - `81638` flatearth plateau broad: `364 / 770` jobs complete (`47%`)
- Interpretation:
  the batch remains healthy and the next true milestone is still the first
  HOMER tagdir, followed by the first HOMER peak BED and the start of staged
  scoring.

## 2026-05-26 Both BAM Sides Substantial, HOMER Still Not Started

- Another gate-focused poll shows the flatearth pair now have substantial BAM
  coverage on both treatment and control sides, but HOMER tag-directory
  creation still has not started.
- Updated flatearth staged counts:
  - flatearth peak narrow:
    - `64` treatment BAMs
    - `59` control BAMs
    - `0` HOMER tag directories
    - `0` HOMER peak BEDs
  - flatearth plateau broad:
    - `64` treatment BAMs
    - `57` control BAMs
    - `0` HOMER tag directories
    - `0` HOMER peak BEDs
- Current progress snapshots from session output:
  - `20626` flatearth peak narrow: `376 / 770` jobs complete (`49%`)
  - `81638` flatearth plateau broad: at least `365 / 770` jobs complete in the
    captured poll window and still actively aligning
- Interpretation:
  the downstream prerequisites are now much more populated on both sides, so
  the absence of HOMER tagdirs is the only meaningful remaining gate before the
  first staged scoring pass can begin.

## 2026-05-26 Both BAM Sides Continue Growing, HOMER Still Not Started

- Latest gate-focused poll still shows no HOMER tag-directory creation, even
  though the flatearth pair now have large BAM counts on both treatment and
  control sides.
- Updated flatearth staged counts:
  - flatearth peak narrow:
    - `75` treatment BAMs
    - `64` control BAMs
    - `0` HOMER tag directories
    - `0` HOMER peak BEDs
  - flatearth plateau broad:
    - `72` treatment BAMs
    - `64` control BAMs
    - `0` HOMER tag directories
    - `0` HOMER peak BEDs
- Current progress snapshots from session output:
  - `20626` flatearth peak narrow: `391 / 770` jobs complete (`51%`)
  - `81638` flatearth plateau broad: at least `376 / 770` jobs complete in the
    captured poll window and still actively aligning
  - hilly pair remain early simulation-phase jobs (`33 / 770`, `4%`)
- Interpretation:
  the staged batch remains healthy, and the only meaningful remaining gate
  before scoring is still the first HOMER tagdir and then the first HOMER peak
  BED.

## 2026-05-26 HOMER Prerequisite Intersection Check

- To distinguish "HOMER has no eligible runs yet" from "Snakemake is still
  prioritizing alignments", checked the intersection of run IDs that already
  have both treatment and control BAMs on disk.
- Results:
  - flatearth peak narrow:
    - `95` run IDs with treatment BAMs
    - `64` run IDs with control BAMs
    - `48` run IDs with both BAM sides already present
  - flatearth plateau broad:
    - `92` run IDs with treatment BAMs
    - `64` run IDs with control BAMs
    - `46` run IDs with both BAM sides already present
- Example shared run IDs already eligible for HOMER:
  - `0025, 0026, 0027, 0028, 0029, 0030, 0031, 0032, 0041, 0042, 0043, 0044, 0045, 0046, 0047`
- Interpretation:
  the absence of HOMER tagdirs is no longer explainable by missing per-run BAM
  prerequisites. The evidence now points to scheduler ordering rather than a
  broken treatment path or missing inputs for specific flatearth runs.

## 2026-05-26 Direct Dry-Run Probe For HOMER Eligibility

- Ran direct Snakemake dry-run probes inside the activated `background_project`
  environment against a specific flatearth run that already has both BAMs:
  - target:
    `results_homer_tfclean_flatearth_peak_narrow_128/0025/peaks/homer/0025_treat_tagdir/tagInfo.txt`
  - target:
    `results_homer_tfclean_flatearth_peak_narrow_128/0025/peaks/homer/0025_peaks.bed`
- Important outcome:
  the dry-run immediately scheduled the expected HOMER rules:
  - `make_homer_tagdirs`
  - `call_peaks_homer`
- The only warning was the expected live-run lock message:
  - `Directory cannot be locked`
  which happened because the real staged Snakemake run is already active on the
  same working directory.
- Dry-run evidence for run `0025`:
  - `make_homer_tagdirs` input:
    - `results_homer_tfclean_flatearth_peak_narrow_128/0025/bowtie2/treat/aligned.sorted.bam`
    - `results_homer_tfclean_flatearth_peak_narrow_128/0025/bowtie2/con/aligned.sorted.bam`
  - `make_homer_tagdirs` output:
    - `results_homer_tfclean_flatearth_peak_narrow_128/0025/peaks/homer/0025_treat_tagdir/tagInfo.txt`
    - `results_homer_tfclean_flatearth_peak_narrow_128/0025/peaks/homer/0025_ctrl_tagdir/tagInfo.txt`
  - `call_peaks_homer` output:
    - `results_homer_tfclean_flatearth_peak_narrow_128/0025/peaks/homer/0025_peaks.txt`
    - `results_homer_tfclean_flatearth_peak_narrow_128/0025/peaks/homer/0025_peaks.bed`
- Interpretation:
  this is strong evidence that the staged workflow is not blocked by a hidden
  HOMER dependency or wrong output path. Eligible runs already exist, and
  Snakemake agrees that HOMER tagdir and peak-calling jobs are the next valid
  downstream work for those runs. The live delay is therefore best explained by
  scheduler ordering inside the already-running batch, not by a broken DAG.

## 2026-05-26 Flatearth HOMER Staged Categories Completed And Scored

- Important runtime milestone:
  the scheduler-ordering delay broke open naturally, and both flatearth staged
  HOMER categories completed end-to-end locally.
- Completion proof on disk:
  - `results_homer_tfclean_flatearth_peak_narrow_128`
    - `128` treatment BAMs
    - `128` control BAMs
    - `256` HOMER tag-info files
    - `128` HOMER peak BEDs
  - `results_homer_tfclean_flatearth_plateau_broad_128`
    - `128` treatment BAMs
    - `128` control BAMs
    - `256` HOMER tag-info files
    - `128` HOMER peak BEDs
- The narrow flatearth session `20626` reached `770 / 770` jobs complete in the
  live session log before closing.
- Local environment note:
  `snakemake_stuff/setup.sh` is cluster-specific here and failed on the Mac
  because it points to `/cvmfs/.../conda` and `$HOME/.conda/envs/background_project`.
  Per repo instructions, the correct local fallback was:
  `eval "$(/opt/anaconda3/bin/conda shell.zsh hook)" && conda activate background_project`.
- Computed staged metrics with `scripts/peak_pr_stats.py` for:
  - `analysis_outputs/homer_tfclean_flatearth_peak_narrow_128_stats`
  - `analysis_outputs/homer_tfclean_flatearth_plateau_broad_128_stats`
- Filter manifest for both stats roots confirms `included_runs: 128`.
- Flatearth narrow headline behavior from
  `analysis_outputs/homer_tfclean_flatearth_peak_narrow_128_stats/group_summary_counts_based.csv`:
  - best F1: `0.5371` at `coverage_ctrl=16`
  - nearby plateaus:
    - `coverage_ctrl=8`: precision `0.5815`, recall `0.4979`, F1 `0.5365`
    - `coverage_ctrl=12`: precision `0.5815`, recall `0.4979`, F1 `0.5365`
    - `coverage_ctrl=24`: precision `0.5815`, recall `0.4979`, F1 `0.5365`
  - low-control failure point:
    - `coverage_ctrl=0.5`: precision `0.5838`, recall `0.2250`, F1 `0.3248`
- Interpretation for flatearth narrow:
  narrow-mode HOMER is clearly functional and stronger than the tiny pilot
  looked, but it is still recall-limited and tops out around mid-`0.53` F1
  rather than looking decisively strong.
- Flatearth broad headline behavior from
  `analysis_outputs/homer_tfclean_flatearth_plateau_broad_128_stats/group_summary_counts_based.csv`:
  - best F1: `0.6556` at `coverage_ctrl=1`
  - other strong points:
    - `coverage_ctrl=4`: precision `0.9787`, recall `0.4833`, F1 `0.6471`
    - `coverage_ctrl=2`: precision `0.9437`, recall `0.4542`, F1 `0.6132`
  - high-control decline remains mostly recall-driven:
    - `coverage_ctrl=16`: precision `0.9818`, recall `0.3396`, F1 `0.5046`
- Interpretation for flatearth broad:
  broad-mode HOMER again looks like the cleaner replacement path, with
  excellent precision and materially stronger F1 than the narrow path.
- Current six-category staged status after scoring the flatearth pair:
  - completed and scored:
    - representative wavy narrow
    - representative wavy broad
    - flatearth narrow
    - flatearth broad
  - still early and not yet scoring-ready:
    - `results_homer_tfclean_realistic_peaks_hilly_narrow_128`
    - `results_homer_tfclean_realistic_plateaus_hilly_broad_128`
    - current artifact counts for both hilly roots remain `0` treatment BAMs,
      `0` control BAMs, and `0` HOMER peak BEDs at this checkpoint
- Decision signal:
  the staged evidence is now internally consistent across two broad-like and
  two narrow-like categories. Broad HOMER continues to look plausibly
  replacement-worthy, while narrow HOMER remains usable but notably weaker and
  still recall-limited.

## 2026-05-26 Hilly HOMER Categories Still In Early Simulation Phase

- Re-checked the two remaining staged `128`-run HOMER categories immediately
  after finishing flatearth scoring to see whether either had become
  scoring-ready.
- Live session snapshots:
  - hilly narrow session `23790`: `45 / 770` jobs complete (`6%`)
  - hilly broad session `91775`: `45 / 770` jobs complete (`6%`)
- Direct artifact counts confirm both are still pre-alignment at this
  checkpoint:
  - `results_homer_tfclean_realistic_peaks_hilly_narrow_128`
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` HOMER tag-info files
    - `0` HOMER peak BEDs
  - `results_homer_tfclean_realistic_plateaus_hilly_broad_128`
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` HOMER tag-info files
    - `0` HOMER peak BEDs
- Interpretation:
  the staged six-category HOMER report is now cleanly split into:
  - four categories already completed and scored locally
  - two hilly categories still in the front simulation phase
- Next gate remains simple:
  wait for the hilly pair to produce downstream BAM/HOMER artifacts, then score
  them and fold them into the first full six-category staged HOMER report.

## 2026-05-26 Hilly Recheck Confirms No Hidden Downstream Progress Yet

- Performed another live/runtime plus filesystem recheck on the two remaining
  hilly staged HOMER categories.
- Updated live session snapshots:
  - hilly narrow session `23790`: `48 / 770` jobs complete (`6%`)
  - hilly broad session `91775`: `47 / 770` jobs complete (`6%`)
- Both live sessions are still executing `simulate_reads_treat` and
  `simulate_reads_con` jobs only; no alignment or HOMER rules have appeared in
  their current logs.
- Direct artifact counts remain unchanged for both hilly roots:
  - `results_homer_tfclean_realistic_peaks_hilly_narrow_128`
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` HOMER tag-info files
    - `0` HOMER peak BEDs
  - `results_homer_tfclean_realistic_plateaus_hilly_broad_128`
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` HOMER tag-info files
    - `0` HOMER peak BEDs
- Interpretation:
  this rules out the possibility that the hilly pair silently became
  scoring-ready between turns. The first full six-category staged HOMER report
  is still waiting entirely on hilly-run runtime progress rather than on any
  new analysis or workflow change.

## 2026-05-26 Staged HOMER Combined-Report Path Wired And Validated

- Since the hilly pair are still runtime-gated, used this turn to remove the
  remaining glue work for the first full six-category staged HOMER report.
- Added `chipseq_pipeline_v2/scripts/build_homer_staged_report.py`.
  - It knows the six expected staged HOMER `128`-run stats roots.
  - It wraps `scripts/build_balanced_288_config_report.py`.
  - It can run in the current partial state or fail strictly with
    `--require-all`.
  - It writes `stage_status.txt` so the combined report root states explicitly
    whether it is partial or complete.
- Added staged HOMER stats-directory mappings to
  `chipseq_pipeline_v2/configs/category_name_map.yaml`, so combined-report
  README output uses canonical category names instead of raw stats-dir names.
- Documented the new wrapper in:
  - `chipseq_pipeline_v2/scripts/README.md`
  - `chipseq_pipeline_v2/README.md`
- Validation:
  - `python -m py_compile chipseq_pipeline_v2/scripts/build_homer_staged_report.py chipseq_pipeline_v2/scripts/build_balanced_288_config_report.py chipseq_pipeline_v2/scripts/category_summary_lib.py`
    passed.
  - Ran:
    `eval "$(/opt/anaconda3/bin/conda shell.zsh hook)" && conda activate background_project && cd chipseq_pipeline_v2 && MPLCONFIGDIR=/private/tmp/mplconfig_homer_stage python scripts/build_homer_staged_report.py --output-dir analysis_outputs/homer_tfclean_128_staged_report_current`
  - Output root created:
    `analysis_outputs/homer_tfclean_128_staged_report_current`
- Validation details for the current partial staged report:
  - `README.md` now lists canonical included categories:
    - `flatearth_peak_narrow`
    - `flatearth_plateau_broad`
    - `wavy_peak_narrow`
    - `wavy_plateau_broad`
  - `stage_status.txt` records:
    - `completed_category_count: 4`
    - `missing_category_count: 2`
    - missing categories:
      - `hilly_peak_narrow`
      - `hilly_plateau_broad`
- Local runtime note:
  setting `MPLCONFIGDIR=/private/tmp/mplconfig_homer_stage` avoids the
  unwritable Matplotlib config path issue during report generation on this Mac.
- Interpretation:
  the first full six-category staged HOMER report is now operationally
  one command away once the two hilly stats roots exist. The remaining blocker
  is purely runtime completion of those two categories, not any report-assembly
  or naming cleanup work.

## 2026-05-26 Staged HOMER Headline Summary Added

- Re-checked the hilly pair first:
  - hilly narrow session `23790`: reached `57 / 770` jobs complete (`7%`)
  - hilly broad session `91775`: reached `57 / 770` jobs complete (`7%`)
  - both are still running only `simulate_reads_treat` / `simulate_reads_con`
    jobs
  - both still have `0` treatment BAMs, `0` control BAMs, `0` HOMER tag-info
    files, and `0` HOMER peak BEDs
- Since scoring is still runtime-blocked, added
  `chipseq_pipeline_v2/scripts/summarize_homer_staged_categories.py` to give a
  compact current readout across the completed staged HOMER categories.
- Important implementation note:
  the first version imported `EXPECTED_HOMER_STATS` from
  `build_homer_staged_report.py`, which accidentally triggered report-building
  side effects because that script executes at import time. Fixed by making the
  summary script own its expected staged-category list instead of importing the
  executable wrapper.
- Documented the new summary helper in:
  - `chipseq_pipeline_v2/scripts/README.md`
  - `chipseq_pipeline_v2/README.md`
- Validation:
  - `python -m py_compile chipseq_pipeline_v2/scripts/summarize_homer_staged_categories.py`
    passed.
  - Ran:
    `eval "$(/opt/anaconda3/bin/conda shell.zsh hook)" && conda activate background_project && cd chipseq_pipeline_v2 && python scripts/summarize_homer_staged_categories.py --output-dir analysis_outputs/homer_staged_summary_current`
  - Output root created:
    `analysis_outputs/homer_staged_summary_current`
- The summary root now contains:
  - `README.md`
  - `category_headlines.csv`
  - `missing_categories.csv`
  - `mode_headlines.csv`
  - `stage_status.txt`
- Current staged HOMER headline signal from
  `analysis_outputs/homer_staged_summary_current/README.md`:
  - completed categories:
    - `flatearth_peak_narrow`: best F1 `0.5371`
    - `flatearth_plateau_broad`: best F1 `0.6556`
    - `wavy_peak_narrow`: best F1 `0.5200`
    - `wavy_plateau_broad`: best F1 `0.6584`
  - missing categories:
    - `hilly_peak_narrow`
    - `hilly_plateau_broad`
  - mode headlines:
    - `broad`: mean best F1 `0.6570` across `2` completed categories
    - `narrow`: mean best F1 `0.5285` across `2` completed categories
- Interpretation:
  the staged summary now makes the current decision signal explicit without
  digging through multiple stats roots: broad HOMER remains materially stronger
  than narrow HOMER even after the flatearth pair completed. The remaining
  blocker for the first full six-category staged verdict is still just the two
  hilly runs finishing.

## 2026-05-26 Staged HOMER Decision Artifact Added

- Re-checked the hilly pair again before making more artifact changes:
  - hilly narrow session `23790`: reached `65 / 770` jobs complete (`8%`)
  - hilly broad session `91775`: reached `65 / 770` jobs complete (`8%`)
  - both still show only simulation-phase jobs in their live logs
  - both still have `0` treatment BAMs, `0` control BAMs, and `0` HOMER peak
    BEDs on disk
- Added an explicit staged decision note to
  `scripts/summarize_homer_staged_categories.py`:
  - new output file:
    `analysis_outputs/homer_staged_summary_current/decision_summary.md`
  - README now also lists `decision_summary.md` in the staged summary outputs
- Validation:
  - `python -m py_compile chipseq_pipeline_v2/scripts/summarize_homer_staged_categories.py`
    passed.
  - Re-ran:
    `eval "$(/opt/anaconda3/bin/conda shell.zsh hook)" && conda activate background_project && cd chipseq_pipeline_v2 && python scripts/summarize_homer_staged_categories.py --output-dir analysis_outputs/homer_staged_summary_current`
- Current explicit staged decision from `decision_summary.md`:
  - broad mean best F1: `0.6570` across `2` completed categories
  - narrow mean best F1: `0.5285` across `2` completed categories
  - broad-minus-narrow mean best F1 gap: `0.1284`
  - decision read:
    - broad HOMER currently looks plausibly replacement-worthy
    - narrow HOMER is functional but still clearly weaker and recall-limited
  - remaining uncertainty:
    - this is still a partial staged verdict because the two hilly categories
      are missing
    - the full six-category decision should be refreshed after those staged
      stats finish
- Interpretation:
  the staged HOMER outputs now provide three layers of ready-to-refresh
  analysis:
  - combined report root
  - headline summary tables
  - explicit staged decision note
  The remaining blocker is still purely runtime progress on the two hilly
  staged runs.

## 2026-05-26 One-Command Staged HOMER Refresh Validated

- Re-checked the hilly pair first:
  - hilly narrow session `23790`: reached `77 / 770` jobs complete (`10%`)
  - hilly broad session `91775`: reached `77 / 770` jobs complete (`10%`)
  - both still show only simulation-phase jobs in the captured live logs
  - both still have `0` treatment BAMs, `0` control BAMs, and `0` HOMER peak
    BEDs on disk
- Added `chipseq_pipeline_v2/refresh_homer_staged_artifacts.sh` as the
  one-command local refresh entrypoint for the staged HOMER artifacts.
- The script:
  - activates the local Mac `background_project` environment through
    `/opt/anaconda3/bin/conda`
  - sets `MPLCONFIGDIR` to `/private/tmp/mplconfig_homer_stage` by default
  - refreshes `analysis_outputs/homer_tfclean_128_staged_report_current`
  - refreshes `analysis_outputs/homer_staged_summary_current`
- Fixed a portability issue found during validation:
  the first version used `BASH_SOURCE[0]`, which failed when run with `zsh`.
  The script now uses `$0` for its root-directory resolution.
- Documented the refresh entrypoint in:
  - `chipseq_pipeline_v2/README.md`
  - `chipseq_pipeline_v2/scripts/README.md`
- Validation:
  - `bash -n chipseq_pipeline_v2/refresh_homer_staged_artifacts.sh` passed.
  - `zsh chipseq_pipeline_v2/refresh_homer_staged_artifacts.sh` completed
    successfully.
  - Refreshed combined-report status still records:
    - `completed_category_count: 4`
    - `missing_category_count: 2`
    - missing categories:
      - `hilly_peak_narrow`
      - `hilly_plateau_broad`
  - Refreshed decision summary still records:
    - broad mean best F1: `0.6570` across `2` completed categories
    - narrow mean best F1: `0.5285` across `2` completed categories
    - broad-minus-narrow mean best F1 gap: `0.1284`
- Interpretation:
  once the hilly stats roots exist, the staged HOMER combined report, headline
  summary, and explicit decision note can now be refreshed together with one
  local command. The remaining gate is still runtime progress on the two hilly
  categories, not report assembly.

## 2026-05-26 Staged HOMER Headline Metric Wording Clarified

- Re-checked the hilly pair first:
  - hilly narrow session `23790`: reached `80 / 770` jobs complete (`10%`)
  - hilly broad session `91775`: reached `79 / 770` jobs complete (`10%`)
  - both still show only simulation-phase jobs in the captured live logs
  - both still have `0` treatment BAMs, `0` control BAMs, `0` HOMER tag-info
    files, and `0` HOMER peak BEDs on disk
- User question clarified that the report phrase "average F1" was ambiguous.
- Updated `scripts/summarize_homer_staged_categories.py` so refreshed staged
  artifacts define the headline metric explicitly as:
  mean of each completed category's best F1 across control-depth settings, not
  the mean F1 across all individual runs.
- Validation:
  - `python -m py_compile chipseq_pipeline_v2/scripts/summarize_homer_staged_categories.py`
    passed.
  - `zsh chipseq_pipeline_v2/refresh_homer_staged_artifacts.sh` completed
    successfully.
- Refreshed `analysis_outputs/homer_staged_summary_current/README.md` now says:
  - `mean of per-category best F1` means to find each category's best F1 across
    control-depth settings, then average those category-best values.
  - It is not the mean F1 across all individual runs.
- Refreshed `analysis_outputs/homer_staged_summary_current/decision_summary.md`
  now carries the same clarification under the current signal.
- The numeric current signal is unchanged:
  - broad mean of per-category best F1: `0.6570`
  - narrow mean of per-category best F1: `0.5285`
  - broad-minus-narrow gap: `0.1284`
- Interpretation:
  the staged decision artifact is now less likely to be misread. The remaining
  scientific gate is unchanged: wait for the two hilly categories to finish,
  score them, then refresh the six-category summary and decision.

## 2026-05-26 Matched MACS2 Control Group Started

- User correctly pointed out that HOMER scoring around the mid-`0.6` range is
  hard to interpret without a matched MACS2 control group on the same staged
  parameter sets.
- Added six matched MACS2 `128`-run control configs mirroring the HOMER staged
  sweep parameters exactly:
  - `configs/macs2_control_tfclean_flatearth_peak_narrow_integrated_128.yaml`
  - `configs/macs2_control_tfclean_flatearth_plateau_broad_integrated_128.yaml`
  - `configs/macs2_control_tfclean_realistic_peaks_wavy_narrow_integrated_128.yaml`
  - `configs/macs2_control_tfclean_realistic_peaks_hilly_narrow_integrated_128.yaml`
  - `configs/macs2_control_tfclean_realistic_plateaus_wavy_broad_integrated_128.yaml`
  - `configs/macs2_control_tfclean_realistic_plateaus_hilly_broad_integrated_128.yaml`
- Each config keeps the same staged HOMER axes:
  - `coverage_treat: [5, 10]`
  - `coverage_ctrl: [0.5, 1, 2, 4, 8, 12, 16, 24]`
  - `seed: [11, 23, 37, 53]`
  - the same category-specific `tf_sigma`, `tf_enrich`, bias, and
    mappability settings
  - the same narrow/broad decoding mode as the planted shape
- Validation:
  - all six matched MACS2 control configs expand to exactly `128` run rows.
  - dry-run target for
    `results_macs2_control_tfclean_flatearth_peak_narrow_128/0001/peaks/macs2/0001_peaks.bed`
    scheduled the expected five-job path:
    `simulate_reads_treat`, `simulate_reads_con`, two `align_bowtie2` jobs,
    and `call_peaks_macs2`.
- Tiny actual MACS2 control smoke:
  - ran one wavy narrow target:
    `results_macs2_control_tfclean_realistic_peaks_wavy_narrow_128/0001/peaks/macs2/0001_peaks.bed`
  - ran one wavy broad target:
    `results_macs2_control_tfclean_realistic_plateaus_wavy_broad_128/0001/peaks/macs2/0001_peaks.bed`
  - both completed end-to-end locally
  - narrow smoke used MACS2 narrow calling
  - broad smoke used MACS2 `--broad`
- Launched the first representative matched MACS2 control batches locally:
  - wavy narrow session `74595`
  - wavy broad session `32239`
- Both representative control batches expanded to `637` Snakemake jobs. This
  is lower than the HOMER `770`-job DAG because MACS2 does not need HOMER
  tag-directory stages.
- Initial runtime checkpoint:
  - wavy narrow MACS2 control was active and past the first simulation jobs
  - wavy broad MACS2 control was active and past the first simulation jobs
- Interpretation:
  this establishes an early MACS2 control group matched to the exact staged
  HOMER parameter sets. The next comparison should score the two representative
  MACS2 controls after completion and compare them against the already-scored
  HOMER wavy narrow and wavy broad categories before interpreting HOMER's
  mid-`0.6` broad result as a caller limit or a sweep/benchmark ceiling.

## 2026-05-26 Representative MACS2 Control Runtime Checkpoint

- Rechecked the two representative MACS2 control sessions:
  - wavy narrow MACS2 control session `74595`
  - wavy broad MACS2 control session `32239`
- Current live progress:
  - wavy narrow: `117 / 637` jobs complete (`18%`)
  - wavy broad: `154 / 637` jobs complete (`24%`)
- Current artifact counts:
  - `results_macs2_control_tfclean_realistic_peaks_wavy_narrow_128`
    - `6` treatment BAMs
    - `12` control BAMs
    - `1` MACS2 normalized peak BED
  - `results_macs2_control_tfclean_realistic_plateaus_wavy_broad_128`
    - `1` treatment BAM
    - `4` control BAMs
    - `1` MACS2 normalized peak BED
- The `params/` manifests are not present yet for either representative control
  root, so full `peak_pr_stats.py` scoring is not ready yet.
- Interpretation:
  the matched MACS2 controls are running and have already produced first peak
  outputs, so the path is not just smoke-clean. Scoring should wait until the
  full `128` peak BEDs and params tables exist.

## 2026-05-26 Representative MACS2 Control Runtime Checkpoint 2

- Rechecked the two representative MACS2 control sessions again:
  - wavy narrow MACS2 control session `74595`
  - wavy broad MACS2 control session `32239`
- Current live progress:
  - wavy narrow: `185 / 637` jobs complete (`29%`)
  - wavy broad: `214 / 637` jobs complete (`34%`)
- Current artifact counts:
  - `results_macs2_control_tfclean_realistic_peaks_wavy_narrow_128`
    - `6` treatment BAMs
    - `20` control BAMs
    - `1` MACS2 normalized peak BED
    - `1` params manifest
  - `results_macs2_control_tfclean_realistic_plateaus_wavy_broad_128`
    - `1` treatment BAM
    - `11` control BAMs
    - `1` MACS2 normalized peak BED
    - `1` params manifest
- Interpretation:
  both representative MACS2 controls are still running and now have params
  manifests, but each still has only the smoke/run-0001 peak BED. Full scoring
  remains gated on all `128` normalized MACS2 peak BEDs.

## 2026-05-26 Representative MACS2 Control Runtime Checkpoint 3

- Rechecked the two representative MACS2 control sessions again:
  - wavy narrow MACS2 control session `74595`
  - wavy broad MACS2 control session `32239`
- Current live progress:
  - wavy narrow: `249 / 637` jobs complete (`39%`)
  - wavy broad: `264 / 637` jobs complete (`41%`)
- Current artifact counts:
  - `results_macs2_control_tfclean_realistic_peaks_wavy_narrow_128`
    - `6` treatment BAMs
    - `28` control BAMs
    - `1` MACS2 normalized peak BED
    - `1` params manifest
  - `results_macs2_control_tfclean_realistic_plateaus_wavy_broad_128`
    - `1` treatment BAM
    - `17` control BAMs
    - `1` MACS2 normalized peak BED
    - `1` params manifest
- Interpretation:
  both matched MACS2 controls are still running through the simulation and
  alignment-heavy front half. Scoring remains gated on all `128` normalized
  MACS2 peak BEDs.

## 2026-05-26 Full Six-Category HOMER Staged Report Completed

- The two previously missing hilly HOMER staged categories now each have all
  `128` HOMER normalized peak BEDs:
  - `results_homer_tfclean_realistic_peaks_hilly_narrow_128`
  - `results_homer_tfclean_realistic_plateaus_hilly_broad_128`
- Scored both hilly categories with `scripts/peak_pr_stats.py`:
  - `analysis_outputs/homer_tfclean_realistic_peaks_hilly_narrow_128_stats`
  - `analysis_outputs/homer_tfclean_realistic_plateaus_hilly_broad_128_stats`
- The first refresh of the combined staged artifacts exposed a report-script
  edge case: `summarize_homer_staged_categories.py` assumed at least one
  missing category and failed with `KeyError: 'canonical_name'` when the report
  became complete.
- Fixed the summary script so the `0`-missing-category state writes an empty
  `missing_categories.csv` with stable columns and so the README only writes
  `- none` under completed categories when no completed categories exist.
- Validation:
  - `python -m py_compile scripts/summarize_homer_staged_categories.py` passed.
  - `zsh refresh_homer_staged_artifacts.sh` completed successfully.
- Refreshed evidence:
  - `analysis_outputs/homer_tfclean_128_staged_report_current/stage_status.txt`
    records `completed_category_count: 6` and `missing_category_count: 0`.
  - `analysis_outputs/homer_staged_summary_current/README.md` records:
    - `flatearth_peak_narrow`: best F1 `0.5371`
    - `flatearth_plateau_broad`: best F1 `0.6556`
    - `hilly_peak_narrow`: best F1 `0.5228`
    - `hilly_plateau_broad`: best F1 `0.6267`
    - `wavy_peak_narrow`: best F1 `0.5200`
    - `wavy_plateau_broad`: best F1 `0.6584`
  - `analysis_outputs/homer_staged_summary_current/mode_headlines.csv`
    records:
    - broad mean of per-category best F1 `0.6469` across `3` categories
    - narrow mean of per-category best F1 `0.5266` across `3` categories
    - broad mean best precision `0.9509`, recall `0.4903`
    - narrow mean best precision `0.5728`, recall `0.4875`
- Interpretation:
  HOMER is now validated through the first full local six-category staged
  report. Broad behavior remains plausibly replacement-worthy, but narrow
  remains weaker and recall-limited. The matched MACS2 control is still needed
  before interpreting the HOMER ceiling as caller-specific versus
  benchmark/simulation-limited.

## 2026-05-26 Representative MACS2 Control Runtime Checkpoint 4

- Rechecked the two representative MACS2 control sessions:
  - wavy narrow MACS2 control session `74595`
  - wavy broad MACS2 control session `32239`
- Current live progress:
  - wavy narrow: reached about `407 / 637` jobs complete (`64%`) in the live
    Snakemake output
  - wavy broad: reached about `392 / 637` jobs complete (`62%`) in the live
    Snakemake output
- Current artifact counts:
  - `results_macs2_control_tfclean_realistic_peaks_wavy_narrow_128`
    - `95` treatment BAMs
    - `65` control BAMs
    - `1` MACS2 normalized peak BED
  - `results_macs2_control_tfclean_realistic_plateaus_wavy_broad_128`
    - `74` treatment BAMs
    - `65` control BAMs
    - `1` MACS2 normalized peak BED
- Interpretation:
  both matched MACS2 controls continue to progress locally. They are still not
  score-ready because MACS2 peak calling has not yet fanned out to all `128`
  runs; scoring remains gated on all `128` normalized MACS2 peak BEDs.

## 2026-05-26 Representative MACS2 Controls Completed and Scored

- Both representative matched MACS2 control batches completed locally:
  - wavy narrow session `74595`: `637 / 637` Snakemake jobs complete
  - wavy broad session `32239`: `637 / 637` Snakemake jobs complete
- Both result roots now contain all `128` normalized MACS2 peak BEDs:
  - `results_macs2_control_tfclean_realistic_peaks_wavy_narrow_128`
  - `results_macs2_control_tfclean_realistic_plateaus_wavy_broad_128`
- Scored both MACS2 controls with `scripts/peak_pr_stats.py`:
  - `analysis_outputs/macs2_control_tfclean_realistic_peaks_wavy_narrow_128_stats`
  - `analysis_outputs/macs2_control_tfclean_realistic_plateaus_wavy_broad_128_stats`
- Matched wavy comparison against the already-scored HOMER categories:
  - MACS2 wavy narrow: best F1 `0.9181` at control `24`
    - precision `0.9258`, recall `0.9104`
  - HOMER wavy narrow: best F1 `0.5200` at control `12`
    - precision `0.5571`, recall `0.4875`
  - MACS2 wavy broad: best F1 `0.9059` at control `4`
    - precision `0.8973`, recall `0.9146`
  - HOMER wavy broad: best F1 `0.6584` at control `1`
    - precision `0.9484`, recall `0.5042`
- Interpretation:
  the HOMER mid-`0.6` broad ceiling and low narrow scores are not simply a
  staged-parameter or benchmark ceiling. On the exact same wavy parameter
  sets, MACS2 recovers both narrow and broad categories much better. This
  changes the next decision: HOMER is integrated and validated, but it does
  not currently look like the best replacement caller. The MACS2 path should
  be prioritized for the larger retry once the six matched MACS2 control
  categories are expanded/scored and the planted-shape/decoder safeguards are
  kept in place.

## 2026-05-26 Remaining Matched MACS2 Control Categories Launched

- Launched the remaining four matched MACS2 `128`-run control categories
  locally at `-j 4` each:
  - flatearth peak narrow session `25320`
  - flatearth plateau broad session `46562`
  - hilly peak narrow session `80065`
  - hilly plateau broad session `75384`
- Each category expanded to a `642`-job DAG:
  - `128` treatment simulation jobs
  - `128` control simulation jobs
  - `256` Bowtie2 alignment jobs
  - `128` MACS2 peak-calling jobs
  - `1` params table job
  - `1` all target
- Initial checkpoint:
  - all four DAGs started cleanly under the local `background_project`
    environment
  - flatearth narrow and flatearth broad had each reached `8 / 642` jobs in
    the captured live logs
  - hilly narrow and hilly broad had started simulation jobs but had not yet
    emitted a numeric progress checkpoint in the captured snippet
- Interpretation:
  this extends the MACS2 control from the two representative wavy categories
  to the full matched six-category control set. Once these four complete,
  score them with `scripts/peak_pr_stats.py` and compare the six-category MACS2
  control against the full six-category HOMER staged report before deciding
  whether HOMER should remain the replacement caller or whether the corrected
  MACS2 path should become the main large-sweep path.

## 2026-05-26 Remaining MACS2 Control Runtime Checkpoint

- Rechecked the four remaining matched MACS2 control sessions:
  - flatearth peak narrow session `25320`
  - flatearth plateau broad session `46562`
  - hilly peak narrow session `80065`
  - hilly plateau broad session `75384`
- Live Snakemake output:
  - flatearth peak narrow reached about `60 / 642` jobs complete (`9%`)
  - flatearth plateau broad reached about `63 / 642` jobs complete (`10%`)
  - both flatearth categories are still in the simulation-heavy front of the
    DAG
  - the two hilly session handles remain active, but the captured output did
    not yet include a numeric progress checkpoint after launch
- Current artifact counts:
  - flatearth peak narrow:
    - `38` treatment read-file sets
    - `40` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - flatearth plateau broad:
    - `29` treatment read-file sets
    - `48` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - hilly peak narrow:
    - `0` treatment read-file sets
    - `0` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - hilly plateau broad:
    - `0` treatment read-file sets
    - `0` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
- Note:
  a direct process-list check was blocked by the local sandbox, but the
  session handles are still active. This is not a scientific or workflow
  blocker; the reliable scoring gate remains the presence of all `128`
  normalized MACS2 peak BEDs per category.

## 2026-05-26 Remaining MACS2 Control Runtime Checkpoint 2

- Rechecked the same four matched MACS2 control sessions again.
- Live Snakemake output:
  - flatearth peak narrow reached about `112 / 642` jobs complete (`17%`)
  - flatearth plateau broad reached about `115 / 642` jobs complete (`18%`)
  - both flatearth categories remain in the simulation-heavy front of the DAG
  - hilly sessions are also producing files now; the earlier zero-output
    checkpoint was a timing/counting artifact, not a stuck-run signal
- Current artifact counts:
  - flatearth peak narrow:
    - `73` treatment read-file sets
    - `74` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - flatearth plateau broad:
    - `71` treatment read-file sets
    - `75` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - hilly peak narrow:
    - `3` treatment read-file sets
    - `1` control read-file set
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - hilly plateau broad:
    - `3` treatment read-file sets
    - `1` control read-file set
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
- Interpretation:
  all four remaining MACS2 control categories are alive and moving, but none is
  score-ready. The next useful action remains waiting for at least one category
  to produce all `128` normalized MACS2 peak BEDs, then scoring it immediately.

## 2026-05-26 Remaining MACS2 Control Runtime Checkpoint 3

- Rechecked the same four matched MACS2 control sessions.
- Live Snakemake output:
  - flatearth peak narrow reached about `184 / 642` jobs complete (`29%`)
  - flatearth plateau broad reached about `186 / 642` jobs complete (`29%`)
  - hilly peak narrow reached at least `4 / 642` jobs complete (`1%`)
  - hilly plateau broad reached at least `4 / 642` jobs complete (`1%`)
- Current artifact counts:
  - flatearth peak narrow:
    - `93` treatment read-file sets
    - `99` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - flatearth plateau broad:
    - `94` treatment read-file sets
    - `98` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - hilly peak narrow:
    - `3` treatment read-file sets
    - `1` control read-file set
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - hilly plateau broad:
    - `3` treatment read-file sets
    - `1` control read-file set
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
- Interpretation:
  the flatearth controls are approaching the transition from read generation
  into alignment, but no remaining category is score-ready yet.

## 2026-05-26 Remaining MACS2 Control Runtime Checkpoint 4

- Rechecked the four active matched MACS2 control sessions.
- Live Snakemake output:
  - flatearth peak narrow reached about `216 / 642` jobs complete (`34%`)
  - flatearth plateau broad reached about `217 / 642` jobs complete (`34%`)
  - hilly peak narrow and hilly plateau broad session handles remain active,
    but did not emit additional output in this polling window
- Current artifact counts:
  - flatearth peak narrow:
    - `113` treatment read-file sets
    - `112` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - flatearth plateau broad:
    - `111` treatment read-file sets
    - `113` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - hilly peak narrow:
    - `3` treatment read-file sets
    - `4` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - hilly plateau broad:
    - `4` treatment read-file sets
    - `4` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
- Interpretation:
  the flatearth pair is close to finishing read generation and should enter
  alignment soon. No remaining category is score-ready yet.

## 2026-05-26 Remaining MACS2 Control Runtime Checkpoint 5

- Rechecked the four active matched MACS2 control sessions.
- Live Snakemake output:
  - flatearth peak narrow reached about `242 / 642` jobs complete (`38%`)
  - flatearth plateau broad reached about `243 / 642` jobs complete (`38%`)
  - hilly peak narrow reached at least `8 / 642` jobs complete (`1%`)
  - hilly plateau broad reached at least `8 / 642` jobs complete (`1%`)
- Current artifact counts:
  - flatearth peak narrow:
    - `126` treatment read-file sets
    - `127` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - flatearth plateau broad:
    - `125` treatment read-file sets
    - `127` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - hilly peak narrow:
    - `4` treatment read-file sets
    - `4` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - hilly plateau broad:
    - `4` treatment read-file sets
    - `4` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
- Interpretation:
  flatearth read generation is almost complete and should start alignment next.
  No remaining category is score-ready yet.

## 2026-05-26 Remaining MACS2 Control Runtime Checkpoint 6

- Rechecked the four active matched MACS2 control sessions.
- Live Snakemake output:
  - flatearth peak narrow reached about `260 / 642` jobs complete (`40%`)
  - flatearth plateau broad reached about `260 / 642` jobs complete (`40%`)
  - both flatearth categories have entered Bowtie2 alignment
  - hilly peak narrow and hilly plateau broad session handles remain active,
    but did not emit additional output in this polling window
- Current artifact counts:
  - flatearth peak narrow:
    - `128` treatment read-file sets
    - `124` control read-file sets
    - `0` treatment BAMs
    - `4` control BAMs
    - `0` MACS2 peak BEDs
  - flatearth plateau broad:
    - `128` treatment read-file sets
    - `124` control read-file sets
    - `0` treatment BAMs
    - `4` control BAMs
    - `0` MACS2 peak BEDs
  - hilly peak narrow:
    - `4` treatment read-file sets
    - `4` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - hilly plateau broad:
    - `4` treatment read-file sets
    - `4` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
- Interpretation:
  flatearth is now in the alignment phase, but no remaining category is
  score-ready. Scoring still requires all `128` normalized MACS2 peak BEDs.

## 2026-05-26 Remaining MACS2 Control Runtime Checkpoint 7

- Rechecked the four active matched MACS2 control sessions.
- Live Snakemake output:
  - flatearth peak narrow reached about `264 / 642` jobs complete (`41%`)
  - flatearth plateau broad reached about `265 / 642` jobs complete (`41%`)
  - hilly peak narrow reached about `12 / 642` jobs complete (`2%`)
  - hilly plateau broad reached about `12 / 642` jobs complete (`2%`)
- Current artifact counts:
  - flatearth peak narrow:
    - `128` treatment read-file sets
    - `119` control read-file sets
    - `0` treatment BAMs
    - `9` control BAMs
    - `0` MACS2 peak BEDs
  - flatearth plateau broad:
    - `128` treatment read-file sets
    - `120` control read-file sets
    - `0` treatment BAMs
    - `8` control BAMs
    - `0` MACS2 peak BEDs
  - hilly peak narrow:
    - `8` treatment read-file sets
    - `4` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - hilly plateau broad:
    - `8` treatment read-file sets
    - `4` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
- Interpretation:
  the flatearth pair is actively aligning control read sets. Treatment BAMs
  and peak BEDs have not started yet, so no category is score-ready.

## 2026-05-26 Remaining MACS2 Control Runtime Checkpoint 8

- Rechecked the four active matched MACS2 control sessions.
- Live Snakemake output:
  - flatearth peak narrow reached about `269 / 642` jobs complete (`42%`)
  - flatearth plateau broad reached about `269 / 642` jobs complete (`42%`)
  - hilly peak narrow and hilly plateau broad session handles remain active,
    but did not emit additional output in this polling window
- Current artifact counts:
  - flatearth peak narrow:
    - `128` treatment read-file sets
    - `115` control read-file sets
    - `0` treatment BAMs
    - `13` control BAMs
    - `0` MACS2 peak BEDs
  - flatearth plateau broad:
    - `128` treatment read-file sets
    - `115` control read-file sets
    - `0` treatment BAMs
    - `13` control BAMs
    - `0` MACS2 peak BEDs
  - hilly peak narrow:
    - `8` treatment read-file sets
    - `4` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - hilly plateau broad:
    - `8` treatment read-file sets
    - `4` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
- Interpretation:
  the flatearth categories continue through control alignment. No remaining
  category is score-ready.

## 2026-05-26 Remaining MACS2 Control Runtime Checkpoint 9

- Rechecked the four active matched MACS2 control sessions.
- Live Snakemake output:
  - flatearth peak narrow reached about `274 / 642` jobs complete (`43%`)
  - flatearth plateau broad reached about `274 / 642` jobs complete (`43%`)
  - hilly peak narrow and hilly plateau broad session handles remain active,
    but did not emit additional output in this polling window
- Current artifact counts:
  - flatearth peak narrow:
    - `128` treatment read-file sets
    - `110` control read-file sets
    - `0` treatment BAMs
    - `18` control BAMs
    - `0` MACS2 peak BEDs
  - flatearth plateau broad:
    - `128` treatment read-file sets
    - `110` control read-file sets
    - `0` treatment BAMs
    - `18` control BAMs
    - `0` MACS2 peak BEDs
  - hilly peak narrow:
    - `8` treatment read-file sets
    - `4` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - hilly plateau broad:
    - `8` treatment read-file sets
    - `4` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
- Interpretation:
  flatearth remains in the control-alignment wave. Treatment BAMs and peak
  BEDs have not started yet, so scoring is still not valid.

## 2026-05-26 Remaining MACS2 Control Runtime Checkpoint 10

- Rechecked the four active matched MACS2 control sessions.
- Live Snakemake output:
  - flatearth peak narrow reached about `279 / 642` jobs complete (`43%`)
  - flatearth plateau broad reached about `280 / 642` jobs complete (`44%`)
  - hilly peak narrow reached about `16 / 642` jobs complete (`2%`)
  - hilly plateau broad reached about `16 / 642` jobs complete (`2%`)
- Current artifact counts:
  - flatearth peak narrow:
    - `128` treatment read-file sets
    - `104` control read-file sets
    - `0` treatment BAMs
    - `24` control BAMs
    - `0` MACS2 peak BEDs
  - flatearth plateau broad:
    - `128` treatment read-file sets
    - `104` control read-file sets
    - `0` treatment BAMs
    - `24` control BAMs
    - `0` MACS2 peak BEDs
  - hilly peak narrow:
    - `8` treatment read-file sets
    - `8` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - hilly plateau broad:
    - `9` treatment read-file sets
    - `7` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
- Interpretation:
  flatearth control alignment is progressing steadily, but no remaining
  category is score-ready.

## 2026-05-26 Remaining MACS2 Control Runtime Checkpoint 11

- Rechecked the four active matched MACS2 control sessions.
- Live Snakemake output:
  - flatearth peak narrow reached about `286 / 642` jobs complete (`45%`)
  - flatearth plateau broad reached about `287 / 642` jobs complete (`45%`)
  - hilly peak narrow and hilly plateau broad session handles remain active,
    but did not emit additional output in this polling window
- Current artifact counts:
  - flatearth peak narrow:
    - `128` treatment read-file sets
    - `97` control read-file sets
    - `0` treatment BAMs
    - `31` control BAMs
    - `0` MACS2 peak BEDs
  - flatearth plateau broad:
    - `128` treatment read-file sets
    - `97` control read-file sets
    - `0` treatment BAMs
    - `31` control BAMs
    - `0` MACS2 peak BEDs
  - hilly peak narrow:
    - `8` treatment read-file sets
    - `8` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - hilly plateau broad:
    - `9` treatment read-file sets
    - `7` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
- Interpretation:
  flatearth remains in control alignment. No remaining category is score-ready.

## 2026-05-26 Remaining MACS2 Control Runtime Checkpoint 12

- Rechecked the four active matched MACS2 control sessions.
- Live Snakemake output:
  - flatearth peak narrow reached about `293 / 642` jobs complete (`46%`)
  - flatearth plateau broad reached about `294 / 642` jobs complete (`46%`)
  - hilly peak narrow and hilly plateau broad session handles remain active,
    but did not emit additional output in this polling window
- Current artifact counts:
  - flatearth peak narrow:
    - `128` treatment read-file sets
    - `90` control read-file sets
    - `0` treatment BAMs
    - `38` control BAMs
    - `0` MACS2 peak BEDs
  - flatearth plateau broad:
    - `128` treatment read-file sets
    - `90` control read-file sets
    - `0` treatment BAMs
    - `38` control BAMs
    - `0` MACS2 peak BEDs
  - hilly peak narrow:
    - `8` treatment read-file sets
    - `8` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - hilly plateau broad:
    - `9` treatment read-file sets
    - `7` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
- Interpretation:
  flatearth is still in the control-alignment wave. No remaining category is
  score-ready.

## 2026-05-26 Remaining MACS2 Control Runtime Checkpoint 13

- Rechecked the four active matched MACS2 control sessions.
- Live Snakemake output:
  - flatearth peak narrow reached about `301 / 642` jobs complete (`47%`)
  - flatearth plateau broad reached about `301 / 642` jobs complete (`47%`)
  - hilly peak narrow reached about `20 / 642` jobs complete (`3%`)
  - hilly plateau broad reached about `20 / 642` jobs complete (`3%`)
- Current artifact counts:
  - flatearth peak narrow:
    - `128` treatment read-file sets
    - `82` control read-file sets
    - `0` treatment BAMs
    - `46` control BAMs
    - `0` MACS2 peak BEDs
  - flatearth plateau broad:
    - `128` treatment read-file sets
    - `82` control read-file sets
    - `0` treatment BAMs
    - `46` control BAMs
    - `0` MACS2 peak BEDs
  - hilly peak narrow:
    - `10` treatment read-file sets
    - `10` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - hilly plateau broad:
    - `11` treatment read-file sets
    - `9` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
- Interpretation:
  flatearth control alignment continues. Hilly simulation is also progressing,
  but no remaining category is score-ready.

## 2026-05-26 Remaining MACS2 Control Runtime Checkpoint 14

- Rechecked the four active matched MACS2 control sessions.
- Live Snakemake output:
  - flatearth peak narrow reached about `308 / 642` jobs complete (`48%`)
  - flatearth plateau broad reached about `309 / 642` jobs complete (`48%`)
  - both flatearth categories have now entered treatment alignment
  - hilly peak narrow and hilly plateau broad session handles remain active,
    but did not emit additional output in this polling window
- Current artifact counts:
  - flatearth peak narrow:
    - `122` treatment read-file sets
    - `80` control read-file sets
    - `6` treatment BAMs
    - `48` control BAMs
    - `0` MACS2 peak BEDs
  - flatearth plateau broad:
    - `122` treatment read-file sets
    - `80` control read-file sets
    - `6` treatment BAMs
    - `48` control BAMs
    - `0` MACS2 peak BEDs
  - hilly peak narrow:
    - `10` treatment read-file sets
    - `10` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - hilly plateau broad:
    - `11` treatment read-file sets
    - `9` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
- Interpretation:
  flatearth has crossed from control-only alignment into treatment alignment.
  MACS2 peak-calling should begin after paired treatment/control BAMs exist for
  individual runs, but no category is score-ready yet.

## 2026-05-26 Remaining MACS2 Control Runtime Checkpoint 15

- Rechecked the four active matched MACS2 control sessions.
- Live Snakemake output:
  - flatearth peak narrow reached about `317 / 642` jobs complete (`49%`)
  - flatearth plateau broad reached about `318 / 642` jobs complete (`50%`)
  - hilly peak narrow and hilly plateau broad session handles remain active,
    but did not emit additional output in this polling window
- Current artifact counts:
  - flatearth peak narrow:
    - `114` treatment read-file sets
    - `80` control read-file sets
    - `14` treatment BAMs
    - `48` control BAMs
    - `0` MACS2 peak BEDs
  - flatearth plateau broad:
    - `114` treatment read-file sets
    - `80` control read-file sets
    - `14` treatment BAMs
    - `48` control BAMs
    - `0` MACS2 peak BEDs
  - hilly peak narrow:
    - `10` treatment read-file sets
    - `10` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - hilly plateau broad:
    - `11` treatment read-file sets
    - `9` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
- Interpretation:
  flatearth treatment alignment is progressing, but MACS2 peak-calling has not
  started yet. No remaining category is score-ready.

## 2026-05-27 Flatearth MACS2 Controls Completed and Scored

- Both flatearth matched MACS2 control batches now completed locally and now
  contain all `128` treatment BAMs, `128` control BAMs, and `128` normalized
  MACS2 peak BEDs:
  - `results_macs2_control_tfclean_flatearth_peak_narrow_128`
  - `results_macs2_control_tfclean_flatearth_plateau_broad_128`
- Scored both flatearth MACS2 controls with `scripts/peak_pr_stats.py`:
  - `analysis_outputs/macs2_control_tfclean_flatearth_peak_narrow_128_stats`
  - `analysis_outputs/macs2_control_tfclean_flatearth_plateau_broad_128_stats`
- Matched flatearth comparison against the already-scored HOMER categories:
  - MACS2 flatearth narrow: best F1 `0.9450` at control `24`
    - precision `0.9420`, recall `0.9479`
  - HOMER flatearth narrow: best F1 `0.5371` at control `16`
    - precision `0.5829`, recall `0.4979`
  - MACS2 flatearth broad: best F1 `0.9233` at control `4`
    - precision `0.8926`, recall `0.9563`
  - HOMER flatearth broad: best F1 `0.6556` at control `1`
    - precision `0.9672`, recall `0.4958`
- Interpretation:
  the flatearth matched-control result reinforces the earlier wavy result:
  HOMER's lower ceiling is not a property of the benchmark alone. On matched
  staged parameter sets, MACS2 substantially outperforms HOMER for both narrow
  and broad categories, driven especially by far stronger recall.

## 2026-05-27 Remaining Hilly MACS2 Control Checkpoint

- Rechecked the two remaining matched MACS2 hilly control categories:
  - `results_macs2_control_tfclean_realistic_peaks_hilly_narrow_128`
  - `results_macs2_control_tfclean_realistic_plateaus_hilly_broad_128`
- Current artifact counts:
  - hilly narrow:
    - `28` treatment read-file sets
    - `28` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
  - hilly broad:
    - `28` treatment read-file sets
    - `28` control read-file sets
    - `0` treatment BAMs
    - `0` control BAMs
    - `0` MACS2 peak BEDs
- Interpretation:
  hilly MACS2 controls are still in read generation and are not score-ready.
  This is a good place to continue using heartbeat-style deferred checks rather
  than frequent polling.

## 2026-05-27 Hilly MACS2 Idle-State Diagnosis and Resume

- Rechecked the same two hilly MACS2 control categories after resuming `/goals`
  work.
- Artifact counts were unchanged from the earlier checkpoint:
  - hilly narrow: `28` treatment read sets, `28` control read sets, `0`
    treatment BAMs, `0` control BAMs, `0` MACS2 peak BEDs
  - hilly broad: `28` treatment read sets, `28` control read sets, `0`
    treatment BAMs, `0` control BAMs, `0` MACS2 peak BEDs
- Latest file modification times in the hilly result trees were from
  `2026-05-26 19:18 PDT`, which is much older than the current morning check on
  `2026-05-27`.
- The matching Snakemake logs:
  - `chipseq_pipeline_v2/.snakemake/log/2026-05-26T183656.154568.snakemake.log`
  - `chipseq_pipeline_v2/.snakemake/log/2026-05-26T183656.154845.snakemake.log`
  stop around `52/642` steps (`8%`) and do not contain a final completion
  footer.
- Interpretation:
  the remaining hilly control runs appear to have gone idle or lost their
  driver process during read generation rather than continuing slowly in the
  background.
- Next action chosen:
  relaunch both hilly MACS2 controls with Snakemake `--rerun-incomplete` so
  the existing partial outputs can be resumed instead of discarded.
- Resume proof:
  both relaunches rebuilt the DAG cleanly with reduced remaining work
  (`589` steps each instead of the full `642`) and immediately entered
  `align_bowtie2`, confirming that partial read-generation outputs were reused
  as intended.
  - hilly narrow resumed into control alignment and began at
    `1/589` completed steps
  - hilly broad resumed into control alignment and began at
    `3/589` completed steps

## 2026-05-27 Hilly MACS2 Resume Milestone

- Follow-up checkpoint on the resumed hilly MACS2 sessions shows both are still
  actively progressing rather than idling again.
- Runtime evidence from the live Snakemake sessions:
  - hilly narrow reached `27/589` completed steps and was writing treatment
    BAMs
  - hilly broad reached `29/589` completed steps and was writing treatment BAMs
- Artifact-state evidence now matches the live session output:
  - remaining treatment `reads_R1.fasta` intermediates dropped from `28` to
    `13` per category as alignments consumed them
  - latest modified files in both trees are fresh treatment BAM/BAM index pairs
- Interpretation:
  the `--rerun-incomplete` resume solved the idle-state problem and moved both
  hilly controls into the alignment phase. They are still far from score-ready,
  so the lowest-token next step remains deferred monitoring until BAM and peak
  production are substantially further along.

## 2026-05-27 Hilly MACS2 Mid-Resume Checkpoint

- Later runtime check confirms both resumed hilly controls are continuing to
  make forward progress without another stall.
- Current artifact counts:
  - hilly narrow:
    - `53` Bowtie2 BAMs
    - `0` MACS2 peak BEDs
    - `0` remaining treatment `reads_R1.fasta` intermediates
  - hilly broad:
    - `53` Bowtie2 BAMs
    - `0` MACS2 peak BEDs
    - `1` remaining treatment `reads_R1.fasta` intermediate
- Live Snakemake session evidence:
  - hilly narrow reached `53/589` completed steps (`9%`)
  - hilly broad reached `53/589` completed steps (`9%`)
  - both sessions are still interleaving `align_bowtie2` with some
    `simulate_reads_{treat,con}` reruns, which is expected under
    `--rerun-incomplete` because missing partial outputs are being rebuilt
    instead of starting from a clean slate
- Interpretation:
  the hilly controls are healthy but still materially distant from the first
  peak-calling milestone. The next useful check should wait for much higher BAM
  completion or the first MACS2 peak BEDs, rather than another near-term poll.

## 2026-05-27 Hilly MACS2 First Peak-Calling Milestone

- Follow-up check reached the next meaningful transition point:
  hilly broad has now begun emitting MACS2 peak calls.
- Current artifact counts:
  - hilly narrow:
    - `53` Bowtie2 BAMs
    - `0` MACS2 peak BEDs
    - `2` remaining treatment `reads_R1.fasta` intermediates
  - hilly broad:
    - `53` Bowtie2 BAMs
    - `1` MACS2 peak BED
    - `3` remaining treatment `reads_R1.fasta` intermediates
- First confirmed broad peak output:
  - `results_macs2_control_tfclean_realistic_plateaus_hilly_broad_128/0014/peaks/macs2/0014_peaks.bed`
- Live-session evidence showed broad entering `call_peaks_macs2` while narrow
  was still interleaving alignment with recovery of a few incomplete simulated
  read outputs.
- Interpretation:
  the resumed hilly controls are approaching the scoreable phase, but they are
  still far from the `128` peak-BED threshold needed for valid category-level
  scoring. The next low-token checkpoint should wait for substantially more
  peak output rather than another immediate poll.

## 2026-05-27 Heartbeat-First Monitoring Rule Recorded

- Added an explicit long-run monitoring rule to the repo root `AGENTS.md`:
  for long local or cluster jobs, do one early checkpoint to estimate the next
  meaningful gate, then defer with heartbeat-style resumes instead of frequent
  polling.
- Synced that rule back into the active `/goals` strategy so this wave uses the
  same default behavior.
- Interpretation:
  future progress on the remaining hilly MACS2 controls should come from later
  milestone checks, not manual near-term status passes unless a real decision
  point is likely to have changed.

## 2026-05-27 Staged HOMER-vs-MACS2 Comparison Helper Added

- Added a reusable comparison script:
  `chipseq_pipeline_v2/scripts/compare_peakcaller_staged_categories.py`
- Purpose:
  read paired HOMER and MACS2 `*_128_stats` roots, extract each category's best
  F1/precision/recall/control setting, compute paired F1 gains, aggregate by
  broad/narrow mode, and emit a staged decision summary that can be refreshed
  as remaining categories finish.
- Also updated:
  `chipseq_pipeline_v2/scripts/README.md`
- Validation:
  - script compiled with `python -m py_compile`
  - script ran successfully to
    `chipseq_pipeline_v2/analysis_outputs/peakcaller_staged_comparison_current`
- Current paired comparison output covers the four already-complete paired
  categories and leaves the two hilly MACS2 categories marked missing.
- Current partial result:
  - completed paired categories: `4`
  - missing paired categories: `2`
  - mean MACS2-minus-HOMER best F1 gain across completed categories: `0.3303`
  - broad gain: `0.2576`
  - narrow gain: `0.4030`
  - MACS2 wins: `4`
  - HOMER wins: `0`
- Interpretation:
  the full six-category caller comparison no longer depends on ad hoc manual
  synthesis. Once the remaining hilly MACS2 stats finish, we can refresh a
  consistent six-category HOMER-vs-MACS2 decision artifact directly.

## 2026-05-27 One-Command Staged Refresh Expanded

- Expanded:
  `chipseq_pipeline_v2/refresh_homer_staged_artifacts.sh`
- New behavior:
  the one-command staged refresh now rebuilds all three current artifact
  families in sequence:
  - `analysis_outputs/homer_tfclean_128_staged_report_current`
  - `analysis_outputs/homer_staged_summary_current`
  - `analysis_outputs/peakcaller_staged_comparison_current`
- Docs updated:
  - `chipseq_pipeline_v2/scripts/README.md`
  - `chipseq_pipeline_v2/README.md`
- Validation:
  - ran `bash chipseq_pipeline_v2/refresh_homer_staged_artifacts.sh`
  - confirmed fresh stage-status timestamps for all three output roots:
    - staged report `2026-05-27 12:12:32`
    - staged HOMER summary `2026-05-27 12:12:33`
    - staged HOMER-vs-MACS2 comparison `2026-05-27 12:12:34`
- Interpretation:
  once the remaining hilly MACS2 category stats are score-ready, the full
  staged reporting and caller-comparison refresh path is now a single command
  rather than a manual multi-step sequence.

## 2026-05-27 Full Six-Category MACS2 Control Comparison Completed

- The remaining hilly MACS2 matched-control categories finished and were scored:
  - `analysis_outputs/macs2_control_tfclean_realistic_peaks_hilly_narrow_128_stats`
  - `analysis_outputs/macs2_control_tfclean_realistic_plateaus_hilly_broad_128_stats`
- Best hilly MACS2 results from the counts-based group summaries:
  - hilly narrow: best F1 `0.8870` at control `24`
  - hilly broad: best F1 `0.8565` at control `24`
- Refreshed staged artifacts with the one-command path, then refreshed again
  after fixing a wording bug in
  `scripts/compare_peakcaller_staged_categories.py` so the completed comparison
  artifact no longer claims more categories are still pending.
- Final six-category paired HOMER-vs-MACS2 comparison from
  `analysis_outputs/peakcaller_staged_comparison_current`:
  - completed paired categories: `6`
  - missing paired categories: `0`
  - mean MACS2-minus-HOMER best F1 gain across all six categories: `0.3192`
  - broad mean best F1: HOMER `0.6469` vs MACS2 `0.8952` (gain `0.2483`)
  - narrow mean best F1: HOMER `0.5266` vs MACS2 `0.9167` (gain `0.3900`)
  - MACS2 category wins: `6`
  - HOMER category wins: `0`
- Category-level final paired outcomes:
  - flatearth peak narrow: HOMER `0.5371` vs MACS2 `0.9450`
  - flatearth plateau broad: HOMER `0.6556` vs MACS2 `0.9233`
  - hilly peak narrow: HOMER `0.5228` vs MACS2 `0.8870`
  - hilly plateau broad: HOMER `0.6267` vs MACS2 `0.8565`
  - wavy peak narrow: HOMER `0.5200` vs MACS2 `0.9181`
  - wavy plateau broad: HOMER `0.6584` vs MACS2 `0.9059`
- Interpretation:
  the staged matched-control evidence is now complete and strongly favors MACS2
  over HOMER for both broad and narrow categories. The local benchmark strategy
  should now treat MACS2 as the validated larger-sweep path for the next
  rerun stage unless a new scientific requirement changes the comparison axis.

## 2026-05-27 Local 288-Run MACS2 Rerun Launcher Added

- Added local sequential launcher:
  `chipseq_pipeline_v2/run_balanced_tfclean_288_local.sh`
- Purpose:
  provide a local-Mac-safe entrypoint for the validated six
  `balanced_tfclean_*_288.yaml` MACS2 configs, defaulting to `--dry-run` and
  allowing explicit full execution with `--run`.
- Supporting docs updated:
  - `chipseq_pipeline_v2/scripts/README.md`
  - `chipseq_pipeline_v2/README.md`
- First validation exposed a real local-only issue:
  Snakemake tried to create source-cache state under the user cache directory
  instead of the writable temp home, causing a permissions failure.
- Fix applied:
  force `HOME=/private/tmp/chipseq_snakemake_home` inside the launcher instead
  of only setting it when `HOME` was unset.
- Second validation exposed an execution-safety issue:
  the `balanced_tfclean_*_288.yaml` configs still default to shared
  `results/` output paths, so a real sequential local run would risk output
  collisions across configs.
- Fix applied:
  the launcher now overrides both `result_root` and `params_table` per config
  at runtime, using config-specific roots like
  `results_balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288`.
- Validation after both fixes:
  ran a representative dry-run for
  `configs/balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288.yaml`
  and Snakemake built the full DAG successfully.
  - resolved jobs: `1442`
  - breakdown:
    - `simulate_reads_treat`: `288`
    - `simulate_reads_con`: `288`
    - `align_bowtie2`: `576`
    - `call_peaks_macs2`: `288`
    - `write_params_table`: `1`
    - `all`: `1`
  - representative isolated output root:
    `results_balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288`
- Interpretation:
  the next larger local MACS2 rerun stage now has a concrete launcher that is
  validated on this Mac and avoids shared-result collisions. The next execution decision is whether to start with
  one or more selected `288`-run categories locally or move directly into the
  full six-config local sequence.

## 2026-05-27 Representative Local 288-Run MACS2 Rerun Launched

- Started the next larger local MACS2 rerun stage with one representative
  category:
  `configs/balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288.yaml`
- Launch path:
  `chipseq_pipeline_v2/run_balanced_tfclean_288_local.sh --run ...`
- Early execution proof:
  - resolved isolated output root:
    `results_balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288`
  - resolved DAG size: `1442` jobs
  - running with `4` cores locally
  - first startup checkpoint already advanced to `8/1442` completed steps
  - active phase at checkpoint: simulated read generation for treatment and
    control
- Interpretation:
  the validated post-staged MACS2 rerun path is now in real local execution,
  not just dry-run readiness. Further monitoring should return to the
  heartbeat/milestone pattern rather than repeated near-term checks.

## 2026-05-27 Local 288-Run Scoring/Report Refresh Helper Added

- Added guarded local report refresher:
  `chipseq_pipeline_v2/refresh_balanced_tfclean_288_local_report.sh`
- Purpose:
  score completed isolated local `results_balanced_tfclean_*_288` MACS2 runs
  with `scripts/peak_pr_stats.py`, then rebuild the combined local 288 report
  with `scripts/build_balanced_288_config_report.py`.
- Supporting docs updated:
  - `chipseq_pipeline_v2/scripts/README.md`
  - `chipseq_pipeline_v2/README.md`
- Guard validation:
  ran the helper against the newly launched representative category
  `balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288` and confirmed
  that it refuses early scoring while output is incomplete:
  - `Not score-ready: results_balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288 has 0/288 MACS2 peak BEDs`
- Interpretation:
  the local 288-run MACS2 stage now has both halves of the execution loop in
  place:
  1. isolated local runner
  2. guarded local scorer/report refresher once a config finishes

## 2026-05-27 Local 288 Helper Ergonomics Tightened

- Added canonical-selector support to both local 288 helpers so they accept
  cleaned category names like `wavy_peak_narrow` and `hilly_plateau_broad`,
  not just raw config/config-root names.
- Validation showed two realistic operator cases:
  - a selector-based dry-run against the already-running representative root
    maps correctly, but Snakemake reports the active directory lock because the
    real run is already in progress there
  - a selector-based report refresh against an unstarted category should fail
    cleanly, not with a raw shell error
- Fix applied:
  `refresh_balanced_tfclean_288_local_report.sh` now checks whether the
  selected local result root exists before calling `find`.
- Validation after the fix:
  - `bash chipseq_pipeline_v2/refresh_balanced_tfclean_288_local_report.sh hilly_plateau_broad`
  - clean result:
    `Not score-ready: results_balanced_tfclean_realistic_plateaus_hilly_broad_integrated_288 does not exist yet`
- Docs updated:
  - `chipseq_pipeline_v2/scripts/README.md`
  - `chipseq_pipeline_v2/README.md`
- Interpretation:
  the local 288-run operator path now behaves more predictably in the two
  common partial-state cases:
  1. selected category already running
  2. selected category not started yet

## 2026-05-27 Pair-Level 288 Batch Selectors Added

- Extended both local 288 helpers with higher-level batch selectors:
  - `flatearth_pair`
  - `wavy_pair`
  - `hilly_pair`
  - `all_six`
- Purpose:
  make the next larger MACS2 rerun stage easier to drive in scientifically
  sensible subsets without manually listing raw config filenames.
- Validation:
  - selector-based dry-run through the runner:
    `bash chipseq_pipeline_v2/run_balanced_tfclean_288_local.sh hilly_pair`
    mapped correctly onto the two hilly configs and built the expected dry-run
    DAGs
  - selector-based guarded refresh:
    `bash chipseq_pipeline_v2/refresh_balanced_tfclean_288_local_report.sh wavy_pair`
    failed cleanly on the still-incomplete local wavy narrow root with:
    `Not score-ready: results_balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288 has 0/288 MACS2 peak BEDs`
- Docs updated:
  - `chipseq_pipeline_v2/scripts/README.md`
  - `chipseq_pipeline_v2/README.md`
- Interpretation:
  the local 288-run control loop is now easier to scale from one category to
  a matched pair, then to all six, without losing the cleaned canonical naming
  discipline introduced earlier in the wave.

## 2026-05-27 Heartbeat-First Monitoring Tightened

- Confirmed the repo-level operating rule is present in `AGENTS.md`:
  long-running local or cluster jobs should use one early checkpoint to
  estimate the next meaningful gate, then defer with heartbeat-style resumes
  instead of frequent polling.
- Updated the active thread heartbeat automation
  (`resume-macs2-control-check`) to make that behavior the default for all
  long-running `/goal` work in this thread, not just the earlier MACS2 control
  wave.
- The heartbeat remains active at a 45-minute interval and now explicitly says
  to avoid extra narration when no new evidence or decision point exists.
- Deliberate follow-through for this turn:
  no new manual poll was performed on the active local 288 runs because there
  was no reason to expect a decision-relevant milestone yet.

## 2026-05-27 Local 288 Progress Summarizer Added

- Added `chipseq_pipeline_v2/scripts/summarize_balanced_tfclean_288_progress.py`
  as a low-token checkpoint helper for active local
  `results_balanced_tfclean_*_288` MACS2 roots.
- Purpose:
  replace repeated manual directory spelunking with one compact milestone
  artifact that reports, per selected config root:
  - result-root existence
  - params-table existence
  - run directory count
  - treatment/control read counts
  - treatment/control BAM counts
  - MACS2 peak BED count
  - coarse phase (`not_started`, `read_generation`, `alignment`,
    `peak_calling`, `score_ready`)
- Default artifact root:
  `analysis_outputs/tfclean_balanced_288_local_progress_current/`
- Validation:
  ran
  `python scripts/summarize_balanced_tfclean_288_progress.py wavy_pair`
  inside the local `background_project` environment and generated the current
  progress artifact successfully.
- Current evidence from that validation:
  - `balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288`:
    `285` run dirs, `259` treatment read sets, `250` control read sets,
    `0` BAMs, `0` peak BEDs, phase `read_generation`
  - `balanced_tfclean_realistic_plateaus_wavy_broad_integrated_288`:
    `205` run dirs, `133` treatment read sets, `136` control read sets,
    `0` BAMs, `0` peak BEDs, phase `read_generation`
- Interpretation:
  the new helper gives heartbeat resumes a single durable place to inspect
  whether the local rerun wave has actually crossed into alignment, peak
  calling, or score-ready output, without spending tokens on repeated
  ad hoc status narration.

## 2026-05-27 Local 288 Progress Refresh Wrapper Added

- Added `chipseq_pipeline_v2/refresh_balanced_tfclean_288_local_progress.sh`
  as the shell-level operator entrypoint for the local `288` progress
  summarizer.
- Purpose:
  keep the local rerun loop ergonomically consistent:
  1. `run_balanced_tfclean_288_local.sh` launches
  2. `refresh_balanced_tfclean_288_local_progress.sh` checkpoints progress
  3. `refresh_balanced_tfclean_288_local_report.sh` scores and rebuilds the
     combined report once a selected root is complete
- Validation:
  ran
  `bash chipseq_pipeline_v2/refresh_balanced_tfclean_288_local_progress.sh wavy_pair`
  and refreshed the progress artifact successfully.
- Updated evidence from the refreshed wavy checkpoint:
  - `balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288`:
    `287` run dirs, `276` treatment read sets, `277` control read sets,
    `0` BAMs, `0` peak BEDs, phase `read_generation`
  - `balanced_tfclean_realistic_plateaus_wavy_broad_integrated_288`:
    `232` run dirs, `173` treatment read sets, `171` control read sets,
    `0` BAMs, `0` peak BEDs, phase `read_generation`
- Interpretation:
  both active local wavy roots are still in the same coarse phase despite
  healthy forward motion, so the scientifically correct next move is still to
  defer further checks until a later heartbeat when alignment or peak-calling
  is more plausible.

## 2026-05-27 Full Local 1728-Wave Status Artifact Tightened

- Extended `scripts/summarize_balanced_tfclean_288_progress.py` so the local
  progress artifact now includes:
  - canonical selector
  - `launch_state` (`not_started`, `in_progress`, `score_ready`, `scored`)
  - `recommended_action`
- Purpose:
  make the full local six-config/`1728` MACS2 wave manageable from one compact
  status artifact rather than only a pair-level partial view.
- Validation:
  ran
  `bash chipseq_pipeline_v2/refresh_balanced_tfclean_288_local_progress.sh all_six`
  and refreshed the six-config artifact successfully.
- Current six-config operator state from that artifact:
  - `wavy_peak_narrow`: `in_progress`, phase `alignment`,
    `288` run dirs, `288` treatment read sets, `284` control read sets,
    `4` control BAMs, `0` peak BEDs
  - `wavy_plateau_broad`: `in_progress`, phase `read_generation`,
    `256` run dirs, `207` treatment read sets, `203` control read sets,
    `0` BAMs, `0` peak BEDs
  - `flatearth_peak_narrow`: `not_started`
  - `flatearth_plateau_broad`: `not_started`
  - `hilly_peak_narrow`: `not_started`
  - `hilly_plateau_broad`: `not_started`
- Interpretation:
  the local `1728` wave is no longer implicit. We now have a durable operator
  view that cleanly distinguishes:
  1. roots already consuming local compute
  2. roots still waiting to be launched
  3. roots that would later become score-ready or fully scored
  This keeps the next launch/scoring decisions grounded in one low-token
  artifact rather than repeated ad hoc inspection.

## 2026-05-27 Remaining-Config Local Launcher Added

- Added `chipseq_pipeline_v2/launch_balanced_tfclean_288_remaining_local.sh`.
- Purpose:
  use the refreshed `all_six` local progress artifact to select only configs
  whose `launch_state` is `not_started`, then forward those canonical
  selectors into `run_balanced_tfclean_288_local.sh`.
- Important local-Mac compatibility fix:
  the first version used `mapfile`, which is not available in the older Bash
  shipped on macOS. Replaced it with a portable `while read` loop.
- Validation:
  ran
  `bash chipseq_pipeline_v2/launch_balanced_tfclean_288_remaining_local.sh --dry-run`
  after the portability fix.
- Verified selector set:
  - `flatearth_peak_narrow`
  - `flatearth_plateau_broad`
  - `hilly_peak_narrow`
  - `hilly_plateau_broad`
  This correctly excluded the already-active wavy pair.
- Verified downstream behavior:
  the first remaining config
  `balanced_tfclean_flatearth_peaks_broad_integrated_288`
  resolved the expected full `1442`-job dry-run DAG under its isolated local
  result root, confirming that the helper hands valid selectors to the
  existing local runner.
- Interpretation:
  the local six-config/`1728` wave now has an operator path not just for
  viewing status, but for safely continuing the remaining untouched configs
  without manual selector bookkeeping or accidental overlap with active roots.

## 2026-05-27 Remaining-Launcher Subset Support Added

- Extended `launch_balanced_tfclean_288_remaining_local.sh` so it no longer
  assumes the full `all_six` wave only. It now accepts subset selectors such
  as `flatearth_pair` or `hilly_pair`, refreshes the matching progress
  artifact, and then intersects that subset with the currently `not_started`
  roots.
- Purpose:
  let the local six-config/`1728` wave continue in scientifically coherent
  chunks without having to choose between:
  1. launching everything untouched at once
  2. manually spelling out per-category selectors
- Validation:
  ran
  `bash chipseq_pipeline_v2/launch_balanced_tfclean_288_remaining_local.sh --dry-run flatearth_pair`
  and confirmed the helper selected exactly:
  - `flatearth_peak_narrow`
  - `flatearth_plateau_broad`
  while still excluding the already-active wavy pair and the untouched hilly
  pair.
- Additional checkpoint evidence from the refreshed progress artifact used by
  that validation:
  - `wavy_peak_narrow` remains `in_progress`, now in phase `alignment`, with
    `18` control BAMs and `288` treatment read sets complete
  - `wavy_plateau_broad` remains `in_progress`, still in
    `read_generation`, with `280` run dirs, `236` treatment read sets, and
    `240` control read sets
- Verified downstream behavior:
  the first selected flatearth config again resolved the expected full
  `1442`-job dry-run DAG under its isolated local result root.
- Interpretation:
  the local rerun wave now has a safer staged-launch control surface:
  not just “launch remaining,” but “launch the remaining configs for this
  scientific subset only.”

## 2026-05-27 Remaining-Launcher Limit Control Added

- Extended `launch_balanced_tfclean_288_remaining_local.sh` with `--limit N`.
- Purpose:
  allow bounded continuation of the local six-config/`1728` wave when we want
  to move only the next one or two untouched roots from a selected ordered
  subset, rather than every currently untouched config in that subset.
- Validation:
  ran
  `bash chipseq_pipeline_v2/launch_balanced_tfclean_288_remaining_local.sh --dry-run --limit 1 flatearth_pair`
  and confirmed the helper selected exactly:
  - `flatearth_peak_narrow`
  not the full flatearth pair.
- Verified downstream behavior:
  the bounded selection still resolves through the existing local runner into
  the expected isolated `1442`-job dry-run DAG for the chosen config.
- Interpretation:
  the operator path can now scale the remaining local wave along three axes:
  1. full remaining set
  2. scientific subset like `flatearth_pair`
  3. bounded prefix of that subset via `--limit N`
  That makes later continuation less all-or-nothing while keeping the
  no-overlap safeguards from the progress artifact.

## 2026-05-27 Bounded Flatearth Continuation Launch Started

- Used the bounded remaining-launch helper to advance the local six-config
  MACS2 wave with one additional real config launch rather than another
  dry-run-only refinement:
  `bash chipseq_pipeline_v2/launch_balanced_tfclean_288_remaining_local.sh --run --limit 1 flatearth_pair`
- Launch selection:
  `Remaining selectors: flatearth_peak_narrow`
- Startup checkpoint:
  the helper forwarded into
  `configs/balanced_tfclean_flatearth_peaks_broad_integrated_288.yaml`
  under the isolated local root
  `results_balanced_tfclean_flatearth_peaks_broad_integrated_288`
- Startup validation proof:
  - Snakemake built the expected `1442`-job DAG
  - execution entered real work immediately rather than failing at parse/setup
  - first finished steps appeared during the startup window, reaching
    `4/1442` completed steps while running read-generation jobs
- Context from the pre-launch low-token checkpoint:
  - `wavy_peak_narrow` was already in `alignment` with `13` treatment BAMs and
    `36` control BAMs
  - `wavy_plateau_broad` was also in `alignment` with `1` control BAM after
    completing all treatment reads and nearly all control reads
  - four roots were still cleanly `not_started`, so launching exactly one
    bounded flatearth root was the conservative next extension of the local
    wave
- Interpretation:
  this turn moved the local six-config/`1728` MACS2 path from “operator
  tooling ready” into the next actual compute step, while still respecting the
  heartbeat-first rule and avoiding a larger uncontrolled batch jump.

## 2026-05-27 Progress Artifact Scoping Fixed

- Found a real operator bug in the local progress artifact flow:
  subset refreshes like `flatearth_pair` were writing to the same
  `analysis_outputs/tfclean_balanced_288_local_progress_current/` root as the
  authoritative `all_six` view.
- Why this mattered:
  the full-wave operator summary could be silently replaced by a subset-only
  table, which is dangerous once the local wave contains a mix of active,
  untouched, and later score-ready roots.
- Fix applied:
  `scripts/summarize_balanced_tfclean_288_progress.py` now treats
  `all_six` as the only selector set that writes to the canonical current root.
  Subset calls automatically write selector-scoped roots such as:
  `analysis_outputs/tfclean_balanced_288_local_progress_flatearth_pair/`
- Validation:
  ran both:
  - `python scripts/summarize_balanced_tfclean_288_progress.py all_six`
  - `python scripts/summarize_balanced_tfclean_288_progress.py flatearth_pair`
- Verified behavior:
  - the `all_six` artifact now correctly preserves the full-wave view:
    `3` active roots (`wavy_peak_narrow`, `wavy_plateau_broad`,
    `flatearth_peak_narrow`) and `3` not-started roots
  - the `flatearth_pair` artifact now lives separately and shows only:
    `flatearth_peak_narrow` in progress and `flatearth_plateau_broad`
    not started
- Current full-wave checkpoint from the refreshed authoritative artifact:
  - `wavy_peak_narrow`: phase `alignment`, `23` treatment BAMs,
    `36` control BAMs
  - `wavy_plateau_broad`: phase `alignment`, `10` control BAMs
  - `flatearth_peak_narrow`: phase `read_generation`, `77` run dirs,
    `53` treatment read sets, `49` control read sets
  - `flatearth_plateau_broad`: `not_started`
  - `hilly_peak_narrow`: `not_started`
  - `hilly_plateau_broad`: `not_started`
- Interpretation:
  the low-token checkpoint system is now safer: subset-focused operator work
  no longer corrupts the authoritative full-wave decision surface used by
  heartbeat resumes.

## 2026-05-27 Batch-Level Decision Artifact Added

- Added `chipseq_pipeline_v2/scripts/build_balanced_tfclean_288_decision_summary.py`.
- Purpose:
  build a heartbeat-friendly operator recommendation from the authoritative
  `all_six` progress CSV rather than forcing future resumes to manually infer
  the batch-level next action from per-root counts.
- Default output root:
  `analysis_outputs/tfclean_balanced_288_local_decision_current/`
- Validation:
  ran
  `python scripts/build_balanced_tfclean_288_decision_summary.py`
  in the local `background_project` environment.
- Current recommendation from the generated decision artifact:
  - next action:
    `defer_until_later_heartbeat_before_more_launches`
  - launch-state counts:
    `in_progress=3`, `not_started=3`
  - phase counts:
    `alignment=2`, `read_generation=1`, `not_started=3`
  - active roots:
    `flatearth_peak_narrow`, `wavy_peak_narrow`,
    `wavy_plateau_broad`
  - not-started roots:
    `flatearth_plateau_broad`, `hilly_peak_narrow`,
    `hilly_plateau_broad`
- Interpretation:
  the new artifact confirms that the correct low-token behavior right now is
  not to extend the batch further yet. We have already moved beyond the
  original wavy pair, but the local wave is still busy enough that the next
  resume should defer until a later milestone rather than launch another root
  immediately.

## 2026-05-27 One-Command Decision Refresh Path Added

- Added `chipseq_pipeline_v2/refresh_balanced_tfclean_288_local_decision.sh`.
- Purpose:
  make the heartbeat checkpoint path explicit and low-friction by refreshing:
  1. the authoritative `all_six` progress artifact
  2. the current batch-level decision summary
  in one command.
- Validation:
  ran
  `bash chipseq_pipeline_v2/refresh_balanced_tfclean_288_local_decision.sh`
  successfully.
- Refreshed combined checkpoint from that run:
  - `flatearth_peak_narrow`: still `in_progress`, phase `read_generation`,
    now `204` run dirs, `137` treatment read sets, `141` control read sets
  - `wavy_peak_narrow`: still `in_progress`, phase `alignment`,
    now `48` treatment BAMs and `36` control BAMs
  - `wavy_plateau_broad`: still `in_progress`, phase `alignment`,
    now `31` control BAMs
  - `flatearth_plateau_broad`, `hilly_peak_narrow`, `hilly_plateau_broad`:
    still `not_started`
- Refreshed decision output:
  - next action:
    `defer_until_later_heartbeat_before_more_launches`
  - launch-state counts:
    `in_progress=3`, `not_started=3`
- Interpretation:
  the standard operator checkpoint is now a single command that both updates
  the evidence and recomputes the recommendation, which is exactly the shape
  we want for low-token heartbeat resumes.

## 2026-05-27 Wavy Local 288 Pair Became Score-Ready And Was Refreshed

- Heartbeat checkpoint through
  `refresh_balanced_tfclean_288_local_decision.sh` showed a real decision
  change:
  - `wavy_peak_narrow` became `score_ready`
  - `wavy_plateau_broad` became `score_ready`
  - batch-level next action changed to `run_score_refresh`
- Follow-through in the same turn:
  ran
  `bash chipseq_pipeline_v2/refresh_balanced_tfclean_288_local_report.sh wavy_pair`
  rather than waiting for another heartbeat.
- Result:
  the guarded scorer/report refresh completed successfully and rebuilt
  `analysis_outputs/tfclean_balanced_288_local_current/` for the two wavy
  categories.
- Included categories in the refreshed larger-scale local report:
  - `wavy_peak_narrow`
  - `wavy_plateau_broad`
- Best aggregated F1 values from the refreshed wavy `288` local report:
  - `balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288`:
    best F1 `0.9278` at treatment coverage `10` and control coverage `8`
  - `balanced_tfclean_realistic_plateaus_wavy_broad_integrated_288`:
    best F1 `0.9393` at treatment coverage `20` and control coverage `4`
- Concurrent wave state after the refresh checkpoint:
  - `flatearth_peak_narrow` remained active in `read_generation`
  - `flatearth_plateau_broad`, `hilly_peak_narrow`,
    `hilly_plateau_broad` remained `not_started`
- Operational note:
  the refresh emitted local Matplotlib/font-cache warnings about unwritable
  user cache paths, but these were non-fatal and did not prevent report
  generation.

## 2026-05-28 Flatearth Peak Local 288 Root Became Score-Ready And Was Refreshed

- Heartbeat checkpoint through
  `refresh_balanced_tfclean_288_local_decision.sh` showed the next real gate:
  `flatearth_peak_narrow` became `score_ready`, while the wavy pair were
  already `scored`.
- Follow-through in the same turn:
  ran
  `bash chipseq_pipeline_v2/refresh_balanced_tfclean_288_local_report.sh flatearth_peak_narrow`
  and completed the guarded report refresh successfully.
- Refreshed larger-scale local report state:
  `analysis_outputs/tfclean_balanced_288_local_current/README.md` now includes:
  - `flatearth_peak_narrow`
- Best aggregated F1 from the refreshed flatearth peak `288` local report:
  - `balanced_tfclean_flatearth_peaks_broad_integrated_288`:
    best F1 `0.9373` at treatment coverage `10` and control coverage `24`
- Decision layer follow-up:
  reran
  `bash chipseq_pipeline_v2/refresh_balanced_tfclean_288_local_decision.sh`
  so the batch-level recommendation reflects the new scored state.
- Updated batch decision after that refresh:
  - next action:
    `launch_when_capacity_allows`
  - launch-state counts:
    `scored=3`, `not_started=3`
  - remaining not-started roots:
    `flatearth_plateau_broad`, `hilly_peak_narrow`,
    `hilly_plateau_broad`
- Interpretation:
  the local six-config wave has now converted the original wavy pair plus one
  bounded flatearth continuation into scored larger-scale outputs. The next
  meaningful move is no longer to score, but to launch one of the remaining
  untouched roots when we want to extend the batch again.

## 2026-05-28 Bounded Flatearth Plateau Continuation Launch Started

- Heartbeat followed the updated decision-layer recommendation
  `launch_when_capacity_allows` by extending the local wave with one more
  bounded untouched root:
  `bash chipseq_pipeline_v2/launch_balanced_tfclean_288_remaining_local.sh --run --limit 1 flatearth_pair`
- Launch selection:
  `Remaining selectors: flatearth_plateau_broad`
- Startup checkpoint:
  the helper forwarded into
  `configs/balanced_tfclean_flatearth_plateaus_broad_integrated_288.yaml`
  under the isolated local root
  `results_balanced_tfclean_flatearth_plateaus_broad_integrated_288`
- Startup validation proof:
  - Snakemake built the expected `1442`-job DAG
  - execution entered real work immediately in read generation
  - first finished steps appeared during the startup window, reaching
    `4/1442` completed steps
- Interpretation:
  the local six-config/`1728` MACS2 wave now has four roots that are no longer
  untouched:
  - `wavy_peak_narrow` scored
  - `wavy_plateau_broad` scored
  - `flatearth_peak_narrow` scored
  - `flatearth_plateau_broad` newly active
  while the remaining untouched roots are now the two hilly categories.

## 2026-05-28 Flatearth Plateau Local 288 Root Became Score-Ready And Was Refreshed

- Heartbeat checkpoint through
  `refresh_balanced_tfclean_288_local_decision.sh` showed
  `flatearth_plateau_broad` had become `score_ready`.
- Follow-through in the same turn:
  ran
  `bash chipseq_pipeline_v2/refresh_balanced_tfclean_288_local_report.sh flatearth_plateau_broad`
  and completed the guarded report refresh successfully.
- Refreshed larger-scale local report state:
  `analysis_outputs/tfclean_balanced_288_local_current/README.md` now includes:
  - `flatearth_plateau_broad`
- Best aggregated F1 from the refreshed flatearth plateau `288` local report:
  - `balanced_tfclean_flatearth_plateaus_broad_integrated_288`:
    best F1 `0.9302` at treatment coverage `20` and control coverage `4`
- Decision layer follow-up:
  reran
  `bash chipseq_pipeline_v2/refresh_balanced_tfclean_288_local_decision.sh`
  so the batch recommendation reflects the new scored state.
- Updated batch decision after that refresh:
  - next action:
    `launch_when_capacity_allows`
  - launch-state counts:
    `scored=4`, `not_started=2`
  - remaining not-started roots:
    `hilly_peak_narrow`, `hilly_plateau_broad`
- Interpretation:
  the local six-config wave now has all wavy and flatearth categories scored
  at `288`, leaving the two hilly categories as the only untouched larger-scale
  local roots.

## 2026-05-28 Bounded Hilly Peak Continuation Launch Started

- Heartbeat followed the updated decision-layer recommendation
  `launch_when_capacity_allows` by extending the local wave into the hilly
  pair with one bounded untouched root:
  `bash chipseq_pipeline_v2/launch_balanced_tfclean_288_remaining_local.sh --run --limit 1 hilly_pair`
- Launch selection:
  `Remaining selectors: hilly_peak_narrow`
- Startup checkpoint:
  the helper forwarded into
  `configs/balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288.yaml`
  under the isolated local root
  `results_balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288`
- Startup validation proof:
  - Snakemake built the expected `1442`-job DAG
  - execution entered real work immediately in read generation
  - early jobs were selected and running during the startup window, including
    treatment and control read-generation tasks
- Interpretation:
  the local six-config/`1728` MACS2 wave now has only one untouched root left:
  `hilly_plateau_broad`.

## 2026-05-28 Hilly Peak Local 288 Root Became Score-Ready And Was Refreshed

- Heartbeat checkpoint through
  `refresh_balanced_tfclean_288_local_decision.sh` showed
  `hilly_peak_narrow` had become `score_ready`.
- Follow-through in the same turn:
  ran
  `bash chipseq_pipeline_v2/refresh_balanced_tfclean_288_local_report.sh hilly_peak_narrow`
  and completed the guarded report refresh successfully.
- Refreshed larger-scale local report state:
  `analysis_outputs/tfclean_balanced_288_local_current/README.md` now includes:
  - `hilly_peak_narrow`
- Best aggregated F1 from the refreshed hilly peak `288` local report:
  - `balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288`:
    best F1 `0.9212` at treatment coverage `10` and control coverage `16`
- Decision layer follow-up:
  reran
  `bash chipseq_pipeline_v2/refresh_balanced_tfclean_288_local_decision.sh`
  so the batch recommendation reflects the new scored state.
- Updated batch decision after that refresh:
  - next action:
    `launch_when_capacity_allows`
  - launch-state counts:
    `scored=5`, `not_started=1`
  - remaining not-started root:
    `hilly_plateau_broad`
- Interpretation:
  the local six-config wave now has every larger-scale category scored except
  the final hilly plateau root.

## 2026-05-28 Final Local 288 Completion Audit And Realstudy Follow-Up

- Heartbeat-driven final local MACS2 `288` scoring completed for the last root:
  `hilly_plateau_broad`.
- Final hilly plateau score refresh proof:
  `bash chipseq_pipeline_v2/refresh_balanced_tfclean_288_local_report.sh hilly_plateau_broad`
  completed successfully.
- Final hilly plateau best aggregated F1 from
  `analysis_outputs/balanced_tfclean_realistic_plateaus_hilly_broad_integrated_288/group_summary.csv`:
  `0.8590` at control coverage `4`.
- Refreshed authoritative full-wave local decision artifact:
  `analysis_outputs/tfclean_balanced_288_local_decision_current/README.md`
  now reports:
  - `next action: complete_or_compare`
  - `launch-state counts: scored=6`
  - `phase counts: score_ready=6`
- Refreshed authoritative full-wave local report artifact:
  `analysis_outputs/tfclean_balanced_288_local_current/README.md`
  now includes all six categories:
  - `flatearth_peak_narrow`
  - `flatearth_plateau_broad`
  - `wavy_peak_narrow`
  - `hilly_peak_narrow`
  - `wavy_plateau_broad`
  - `hilly_plateau_broad`
- Completion-audit conclusion for the local controlled wave:
  the replacement-caller retry objective is now fully evidenced through
  matched staged `128` HOMER-vs-MACS2 comparison plus the all-six local MACS2
  `288` rerun wave.
- Follow-up on the still-open HPC realstudy front:
  the old cluster failure root cause remains relevant in the current files.
  `chipseq_pipeline_v2_realstudy/config.yaml` still used directory-style
  Bowtie2 index values (`indexes/ce11/bowtie2_index`) while the local full-run
  override already used basename-prefix values (`.../ce11`).
- Repair applied:
  updated the default realstudy production config to use explicit Bowtie2
  basename prefixes:
  - `indexes/ce11/bowtie2_index/ce11`
  - `indexes/mm10/bowtie2_index/mm10`
- Documentation repair applied in
  `chipseq_pipeline_v2_realstudy/docs/chips_realsim_workflow.md` so the
  required `bowtie2 -x` basename-prefix form is explicit and the old
  directory-only interpretation is less likely to recur.
- Validation after the realstudy config/doc repair:
  `cd chipseq_pipeline_v2 && python -m py_compile scripts/*.py`
  passed in `background_project`.
- Realstudy cluster-shape dry-run on the patched config still does not finish
  locally, but the remaining failure moved to missing production reference
  assets (`references/ce11/genome.fa`) rather than the previously documented
  Bowtie2 prefix mismatch. That means the local evidence gap is now staged
  assets, not another discovered config-shape bug.
