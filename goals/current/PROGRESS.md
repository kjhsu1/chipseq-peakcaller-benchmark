# Progress

## ChIPs Ontology Analysis and Cluster Sweep Wave

- [x] New isolated branch/worktree created:
  `codex/chips-ontology-analysis`.
- [x] Existing `/goals` history preserved; new work is appended, not replacing
  previous waves.
- [x] Current ChIPs realstudy workflow inspected.
- [x] Existing ontology helper scripts located.
- [x] Latest successful 128-core Slurm sweep script patterns inspected.
- [x] Scientific validity report written before implementation.
- [x] Experimental plan adjusted from the report.
- [x] Ontology downstream Snakemake step implemented.
- [x] Ontology dry-run target passes.
- [x] 128-core / 2TB-RAM Slurm script added.
- [x] Matching shell submission script added.
- [x] Compile/tests/dry-runs complete.

## ChIPs Realstudy Bug-Fix Wave

- [x] Bug hunt completed on the ChIPs-based realstudy path.
- [x] Local realstudy tests rerun and confirmed green: `17 passed`.
- [x] Local ChIPs dry-run rechecked and confirmed broken at the ingest/alignment
  handoff.
- [x] Ingest download outputs connected to alignment inputs in the DAG.
- [x] Processed-BAM and FASTQ alignment paths separated cleanly.
- [x] Parse-time stale metadata loading removed or neutralized.
- [x] Local reference/index override path unified across the workflow.
- [x] Assembly-aware MACS2 genome-size handling added.
- [x] ChIPs learn-model replicate/peak-input logic made internally consistent.
- [x] Manifest bookkeeping corrected for actual local-file presence in runtime logic.
- [x] Run ID scheme hardened against future parameter growth.
- [x] Local dry-run upgraded from failing to passing for the ChIPs path.
- [x] End-to-end local-runnability status re-evaluated and documented clearly.

- [x] Goal initialized.
- [x] Current repo state inspected.
- [x] Simulator literature/web pass completed.
- [x] Low-utility PMF/visual helper scripts archived from active path.
- [x] Controlled-sweep repair context documented.
- [x] Simulator comparison and recommendation documented.
- [x] Revised 11-step plan documented.
- [x] ChIPs workflow assets added.
- [x] EPIC2 sweep assets added.
- [x] AGENTS.md updated for `/goals` operation and repo hygiene.
- [x] Active scripts README/index added.
- [x] Local unit and compile verification completed.
- [x] Snakemake dry-run verification completed.
- [x] EPIC2 expected run tables generated and checked.

## Current Center of Gravity

The current center of gravity is the new ChIPs realstudy bug-fix wave. The
Python tests are green, and the realstudy ChIPs Snakemake path is now local
dry-run clean. ChIPs `2.4` is now installed in `background_project`, and a tiny
WCE/control ChIPs smoke target runs locally. Full production execution is still
cluster-oriented because it requires full real-study downloads and production
reference/index assets.

## ChIPs Environment and Tiny Smoke Wave

- [x] Confirmed `background_project` did not expose `chips` at wave start.
- [x] Installed or exposed ChIPs inside `background_project`.
- [x] Added or selected a tiny local smoke config.
- [x] Re-ran realstudy compile checks.
- [x] Re-ran realstudy tests.
- [x] Ran a tiny ChIPs Snakemake dry-run.
- [x] Ran one tiny ChIPs output target, if safe.
- [x] Updated the local-runnability assessment.

## Local and Slurm ChIPs Execution Wave

- [x] Inspected current ChIPs launch/config files.
- [x] Added local full-run config for staged local assets.
- [x] Added local Snakemake launcher.
- [x] Updated local-vs-Slurm documentation.
- [x] Re-ran compile and tests.
- [x] Validated local launcher dry-run behavior.
- [x] Rechecked tiny smoke execution or output.
- [x] Logged final local and cluster run commands.
