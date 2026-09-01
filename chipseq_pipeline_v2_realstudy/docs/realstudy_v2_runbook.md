# Realstudy v2 Future-Operator Runbook

These commands were documented but not executed on the implementation branch.

```bash
eval "$(/opt/anaconda3/bin/conda shell.zsh hook)"
conda activate background_project
cd chipseq_pipeline_v2_realstudy
python -m py_compile scripts/*.py
pytest tests/test_realstudy_v2_*.py
HOME=/private/tmp/chipseq_snakemake_home ./run_realstudy_v2_local.sh
```

Next, stage all seven FASTQs at the exact paths in
`manifests/realstudy_v2_files.csv`, plus the configured ce11 FASTA, Bowtie2 index
basename, and GFF. Local-only staging is the default. An operator can explicitly
enable registered downloads with
`--config realstudy_require_local_inputs=false` when policy permits.

Execute locally with `MODE=run CORES=8 ./run_realstudy_v2_local.sh`, or set
`PIPELINE_DIR` and a durable `OUTPUT_ROOT` then run
`bash slurm/submit_realstudy_v2.sh` on the cluster.

The final target is `validation/final_integrity_report.json`. It passes only
when the 44-run database, constraints, CSV checksums, two decisions, and all six
figure packages are complete. A zero-line MACS2 peak file is valid; a missing
file or nonzero command exit is a failure.
