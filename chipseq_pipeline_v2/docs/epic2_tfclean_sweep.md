# EPIC2 TF-Clean Comparison Sweep

## Purpose

This sweep asks whether the repaired control-depth conclusions are specific to
MACS2 or also appear under EPIC2, a broad-domain-oriented caller. The EPIC2
configs mirror the six repaired TF-clean configs and only switch
`peakcaller_list` from `macs2` to `epic2`.

## Scope

- Six categories.
- 288 runs per category.
- Total resolved design: `6 x 288 = 1728` runs.
- Local Mac validation is limited to config parsing, run-table generation, path
  expansion, and Snakemake dry-runs.
- Actual EPIC2 execution is cluster-only.

## Configs

- `configs/epic2_tfclean_flatearth_peaks_broad_integrated_288.yaml`
- `configs/epic2_tfclean_flatearth_plateaus_broad_integrated_288.yaml`
- `configs/epic2_tfclean_realistic_peaks_wavy_narrow_integrated_288.yaml`
- `configs/epic2_tfclean_realistic_peaks_hilly_narrow_integrated_288.yaml`
- `configs/epic2_tfclean_realistic_plateaus_wavy_broad_integrated_288.yaml`
- `configs/epic2_tfclean_realistic_plateaus_hilly_broad_integrated_288.yaml`

## Cluster Launch

From `chipseq_pipeline_v2` on the cluster:

```bash
source ../snakemake_stuff/setup.sh
bash slurm/submit_epic2_tfclean_1728_series.sh
```

The submitter chains the six category jobs. Each job runs one 288-run config
through `slurm/epic2_tfclean_288_singlejob.sbatch`.

## After Cluster Execution

Regenerate per-category performance tables and plots with:

```bash
python scripts/peak_pr_stats.py --input-dir results --params-csv results/params/epic2_tfclean_<category>_run_params.csv --out-dir analysis_outputs/epic2_tfclean_<category>
python scripts/build_balanced_288_config_report.py --input-root analysis_outputs --out-root analysis_outputs/epic2_tfclean_1728_report
```

Final EPIC2 deliverables are the computed performance table, six PR/recall/F1
vs control-coverage plots, and a short MACS2-vs-EPIC2 comparison note.

## Local Validation

The intended local checks are:

```bash
eval "$(/opt/anaconda3/bin/conda shell.zsh hook)"
conda activate background_project
python -m py_compile scripts/*.py
HOME=/private/tmp/chipseq_snakemake_home snakemake -s Snakefile.py \
  --configfile configs/epic2_tfclean_realistic_peaks_wavy_narrow_integrated_288.yaml \
  --dry-run
```

The dry-run validates config parsing, run-table expansion, EPIC2 output paths,
and non-EPIC2 workflow plumbing. It does not run EPIC2 locally.

Generate expected EPIC2 run tables locally with:

```bash
HOME=/private/tmp/chipseq_snakemake_home snakemake \
  results/params/epic2_tfclean_realistic_peaks_wavy_narrow_integrated_288_run_params.csv \
  -s Snakefile.py \
  --configfile configs/epic2_tfclean_realistic_peaks_wavy_narrow_integrated_288.yaml \
  --cores 1
```
