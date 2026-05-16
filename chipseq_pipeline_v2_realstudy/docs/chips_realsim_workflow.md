# ChIPs Realistic Simulation Workflow

## Purpose

This workflow replaces the naive real-pileup-to-PMF idea with a published
ChIP-seq simulator path. Real studies are used for parameter learning and
truth-proxy/reference regions, while ChIPs generates new treatment and WCE
control reads.

## Default Study

The default ChIPs target is the locally downloaded clean worm broad-mark study:

- `geo_gse67028_celegans_h3k9me2_adult`

The ENCODE SKN-1 study remains a narrow reference asset, but it is not included
in the default ChIPs target list until local BAM inputs and peak calls are ready.

## Workflow Shape

1. Align selected downloaded FASTQs into `data/aligned/{study_id}/{role}/`.
2. Call MACS2 peaks on the selected real-study treatment/control BAMs.
3. Learn a ChIPs model from the selected treatment BAM and peak BED.
4. Generate treatment FASTQs from ChIPs peak mode.
5. Generate control FASTQs from ChIPs WCE mode.
6. Align simulated FASTQs with Bowtie2.
7. Call MACS2 peaks on simulated treatment/control BAMs.

## Inputs

- `metadata/data_manifest.csv`
- `metadata/prototype_run_table.csv`
- `references/{assembly}/genome.fa`
- `indexes/{assembly}/bowtie2_index`
- downloaded FASTQs or BAMs listed in the manifest

Reference and index paths are configured under `chips.reference_fasta_by_assembly`
and `chips.bowtie2_index_by_assembly`. The default config uses production-style
paths such as `references/ce11/genome.fa`. The local dry-run override
`configs/chips_local_dryrun.yaml` points to the existing toy `ce11_1pct`
reference only to validate workflow shape on a Mac; it is not the production
reference for scientific execution.

## Outputs

- `analysis_outputs/chips_models/{study_id}/{study_id}.json`
- `results_chips/{run_id}/treat/reads_R1.fastq`
- `results_chips/{run_id}/treat/reads_R2.fastq`
- `results_chips/{run_id}/con/reads_R1.fastq`
- `results_chips/{run_id}/con/reads_R2.fastq`
- `results_chips/{run_id}/bowtie2/{cond}/aligned.sorted.bam`
- `results_chips/{run_id}/peaks/macs2/{run_id}_peaks.bed`

## Cluster Launch

From `chipseq_pipeline_v2_realstudy` on the cluster:

```bash
source ../snakemake_stuff/setup.sh
bash slurm/submit_chips_realsim.sh
```

The Slurm wrapper enables `enable_chips_targets=true` at runtime. Local Mac
validation should use Snakemake dry-runs and tests; full ChIPs execution is
cluster-oriented.

Local dry-run:

```bash
eval "$(/opt/anaconda3/bin/conda shell.zsh hook)"
conda activate background_project
HOME=/private/tmp/chipseq_snakemake_home snakemake -s Snakefile.py \
  --configfile config.yaml configs/chips_local_dryrun.yaml \
  --config enable_chips_targets=true \
  --dry-run
```
