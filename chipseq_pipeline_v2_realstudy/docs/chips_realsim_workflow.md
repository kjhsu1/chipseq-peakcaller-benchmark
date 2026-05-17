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

For full local execution, stage assets under ignored local paths and use
`configs/chips_local_full.yaml`:

- `data/raw/geo_gse67028_celegans_h3k9me2_adult/SRR1917669.fastq.gz`
- `data/raw/geo_gse67028_celegans_h3k9me2_adult/SRR1917670.fastq.gz`
- `data/raw/geo_gse67028_celegans_h3k9me2_adult/SRR1917671.fastq.gz`
- `data/local/references/ce11/genome.fa`
- `data/local/references/ce11/genome.fa.fai`
- `data/local/indexes/ce11/bowtie2/ce11.*.bt2`

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
validation should use Snakemake dry-runs, tests, and the tiny ChIPs smoke config
before attempting larger runs. The workflow shape now dry-runs locally,
including ingest downloads, real-study BAM merging, ChIPs learning, simulated
read generation, alignment, and MACS2 peak calling.

Validated local pieces:

- ChIPs `2.4` is available from `background_project` after source-build
  installation.
- `configs/chips_tiny_smoke.yaml` resolves one tiny run and can generate a WCE
  control FASTQ pair without downloading real-study files.

Full local execution still requires real Bowtie2 index files at the configured
local paths and network access for any missing real-study downloads. Until those
are present, full ChIPs execution remains cluster-oriented.

## Local Launch

Local execution uses the same Snakemake workflow, without Slurm:

```bash
./run_chips_realsim_local.sh
```

The local launcher defaults to `MODE=dry-run`, `CORES=4`, and:

```text
CONFIG_FILES="config.yaml configs/chips_local_full.yaml"
```

After staging the required local assets, run the same workflow locally with:

```bash
MODE=run CORES=4 ./run_chips_realsim_local.sh
```

To target a specific output instead of the full ChIPs target set:

```bash
TARGETS="results_chips/<run_id>/con/reads_R1.fastq results_chips/<run_id>/con/reads_R2.fastq" \
  ./run_chips_realsim_local.sh
```

Local dry-run:

```bash
eval "$(/opt/anaconda3/bin/conda shell.zsh hook)"
conda activate background_project
HOME=/private/tmp/chipseq_snakemake_home snakemake -s Snakefile.py \
  --configfile config.yaml configs/chips_local_dryrun.yaml \
  --config enable_chips_targets=true \
  --dry-run
```

Tiny local ChIPs smoke run:

```bash
eval "$(/opt/anaconda3/bin/conda shell.zsh hook)"
conda activate background_project
HOME=/private/tmp/chipseq_snakemake_home snakemake -s Snakefile.py \
  --configfile config.yaml configs/chips_tiny_smoke.yaml \
  --cores 1 \
  results_chips/geo_gse67028_celegans_h3k9me2_adult_t0p001_c0p001_s11_fl150_rl38_albowtie2_pcmacs2_mbroad/con/reads_R1.fastq \
  results_chips/geo_gse67028_celegans_h3k9me2_adult_t0p001_c0p001_s11_fl150_rl38_albowtie2_pcmacs2_mbroad/con/reads_R2.fastq
```
## Ontology Downstream Analysis

Set `enable_chips_ontology_targets=true` to extend the ChIPs real-study DAG through downstream ontology analysis. The ontology branch scores each simulated ChIPs run against the real-study template peak BED that was used for `chips simreads -t bed`, classifies each template region with `scripts/classify_regions.py`, combines all classified tables, and writes aggregate summaries under `analysis_outputs/chips_ontology/summary/`.

This is an empirical truth-template recovery analysis. ChIPs provides learned simulation from real ChIP-seq signal and can generate both treatment-like and whole-cell-extract control reads, but the workflow does not assume ChIPs exposes a closed-form intrinsic PMF for each parameter set.

For cluster validation at the 128-core/2TB scale, use:

```bash
slurm/submit_chips_realsim_ontology_128cpu_2tb.sh
```

The submission script follows the latest successful cluster pattern: copy the pipeline to a per-job work directory, remove Python bytecode/Snakemake state, symlink high-volume outputs into the archived result directory, then run Snakemake in the `background_project` Conda environment.
