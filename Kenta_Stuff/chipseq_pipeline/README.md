# Chipseq Pipeline

This directory contains a Snakemake workflow for simulating CHIP-seq reads,
aligning them and calling peaks.

## Layout
- `Snakefile.py` – entry point for the workflow.
- `config.yaml` – configuration options.
- `envs/` – Conda environment definitions.
- `rules/` – Snakemake rule files.
- `scripts/` – helper scripts used by the rules.
- `tests/` – unit tests for pipeline utilities.

## Running Tests
```bash
source ../snakemake_stuff/setup.sh
pytest tests
```

