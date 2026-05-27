# Chipseq Pipeline Realstudy

This sibling pipeline is the prototype real-study-conditioned benchmark track.

Execution notes:

```bash
source ../snakemake_stuff/setup.sh
pytest tests
```

Current scope:

- explicit study selection manifests with accepted, rejected, fallback, and pending rows
- resumable real-data ingestion scaffolding
- reference-intensity model builders that use estimated or truth-proxy language only
- coverage-sweep sampling metadata and ontology-based evaluation scaffolding
- ChIPs-based realistic simulation workflow assets that learn from selected
  real-study BAM/peak inputs, simulate treatment and WCE control FASTQs, and
  hand outputs back to the alignment/peak-calling workflow shape
- production validation and reproducibility hooks for the ChIPs ontology run,
  including a post-run validation report, runtime/resource capture, and
  reproducibility package under the final output root
- production real-study input handling now defaults to pre-staged local inputs:
  populate the manifest-listed BAM/FASTQ files under `data/raw/` before launch
- workflow entrypoints that mirror the controlled pipeline structure without
  forcing a shared-framework refactor in this first pass

Key docs:

- `docs/chips_realsim_workflow.md`

Production validation command:

```bash
source ../snakemake_stuff/setup.sh
python scripts/validate_chips_ontology_production_run.py \
  --run-table metadata/prototype_run_table.csv \
  --output-root /path/to/chips_run_output_root \
  --write-report /path/to/chips_run_output_root/reproducibility/validation_report.md
```

Local and cluster execution:

- local Mac validation now supports tests, full Snakemake dry-run shape for the
  ChIPs path, and a tiny one-rule ChIPs WCE/control smoke run
- the local `background_project` environment has been validated with ChIPs
  `2.4`
- full local execution is available through `run_chips_realsim_local.sh` once
  the real-study FASTQs/BAMs plus full reference/index assets are staged locally
- full ChIPs execution should run on the cluster through
  `slurm/submit_chips_realsim.sh`

Local dry-run:

```bash
./run_chips_realsim_local.sh
```

Local full run after staging assets:

```bash
MODE=run CORES=4 ./run_chips_realsim_local.sh
```
