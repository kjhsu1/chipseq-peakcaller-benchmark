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
- workflow entrypoints that mirror the controlled pipeline structure without
  forcing a shared-framework refactor in this first pass

Key docs:

- `docs/chips_realsim_workflow.md`

Blocker boundary:

- local Mac validation stops before full ChIPs/Snakemake execution if those
  tools are not installed locally
- full ChIPs execution should run on the cluster through
  `slurm/submit_chips_realsim.sh`
