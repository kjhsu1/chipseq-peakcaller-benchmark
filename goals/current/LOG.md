# Log

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
