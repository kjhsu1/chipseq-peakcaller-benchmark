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

Blocker boundary:

- implementation stops before live ENCODE or NCBI fetches are executed
- those steps require plain network access at minimum
