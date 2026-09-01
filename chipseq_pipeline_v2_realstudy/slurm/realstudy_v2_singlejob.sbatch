#!/usr/bin/env bash
#SBATCH --job-name=realstudy-v2
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=7-00:00:00
#SBATCH --output=slurm_logs/realstudy-v2-%j.out
#SBATCH --error=slurm_logs/realstudy-v2-%j.err

set -euo pipefail

PIPELINE_DIR="${PIPELINE_DIR:?Set PIPELINE_DIR to the staged chipseq_pipeline_v2_realstudy directory}"
OUTPUT_ROOT="${OUTPUT_ROOT:?Set OUTPUT_ROOT to a durable cluster output directory}"
CORES="${SLURM_CPUS_PER_TASK:-32}"

source "${PIPELINE_DIR}/../snakemake_stuff/setup.sh"
cd "${PIPELINE_DIR}"
mkdir -p slurm_logs "${OUTPUT_ROOT}"

snakemake -s Snakefile.py \
  --configfile config.yaml configs/realstudy_v2.yaml \
  --config realstudy_v2_output_root="${OUTPUT_ROOT}" \
  --cores "${CORES}" --rerun-incomplete --printshellcmds
