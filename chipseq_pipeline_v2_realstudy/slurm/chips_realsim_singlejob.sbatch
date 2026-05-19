#!/bin/bash
#SBATCH --job-name=chips_realsim
#SBATCH --account=publicgrp
#SBATCH --partition=high
#SBATCH --output=/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2_realstudy/slurm_jobs/%j/%x.out
#SBATCH --error=/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2_realstudy/slurm_jobs/%j/%x.err
#SBATCH --chdir=/home/kjhsu/Code/Korflab/chipseq-peakcaller-benchmark/chipseq_pipeline_v2_realstudy
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4-00:00:00

set -euo pipefail

CONDA_BIN="${CONDA_BIN:-/cvmfs/hpc.ucdavis.edu/sw/conda/root/bin/conda}"
ENV_PREFIX="${ENV_PREFIX:-/home/kjhsu/.conda/envs/background_project}"
PIPE_SRC="${PIPE_SRC:-/home/kjhsu/Code/Korflab/chipseq-peakcaller-benchmark/chipseq_pipeline_v2_realstudy}"
DATA_ROOT="${DATA_ROOT:-/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2_realstudy}"
CONFIG_REL="${CONFIG_REL:-config.yaml}"

RUN_TS="$(date +%Y%m%d_%H%M%S)"
JOB_ROOT="${DATA_ROOT}/work/${SLURM_JOB_ID}/pipeline"
DEST_ROOT="${DATA_ROOT}/chips_runs/chips_realsim_${RUN_TS}"
LOG_ROOT="${DATA_ROOT}/logs"
LOG_FILE="${LOG_ROOT}/chips_realsim_${RUN_TS}.log"

mkdir -p "${DATA_ROOT}/slurm_jobs/${SLURM_JOB_ID}" "${JOB_ROOT}" "${DEST_ROOT}/results_chips" "${LOG_ROOT}"
mkdir -p "/tmp/kjhsu_snakemake_cache/${SLURM_JOB_ID}"

echo "[$(date)] START ChIPs realstudy workflow" | tee -a "${LOG_FILE}"
echo "results=${DEST_ROOT}" | tee -a "${LOG_FILE}"

export XDG_CACHE_HOME="/tmp/kjhsu_snakemake_cache/${SLURM_JOB_ID}"
export PYTHONDONTWRITEBYTECODE=1

cp -a "${PIPE_SRC}/." "${JOB_ROOT}/"
find "${JOB_ROOT}" \( -type d -name '__pycache__' -o -type f -name '*.pyc' \) -exec rm -rf {} +
rm -rf "${JOB_ROOT}/results_chips" "${JOB_ROOT}/.snakemake"
ln -sfn "${DEST_ROOT}/results_chips" "${JOB_ROOT}/results_chips"

cd "${JOB_ROOT}"

set +e
"${CONDA_BIN}" run -p "${ENV_PREFIX}" snakemake \
  -j "${SLURM_CPUS_PER_TASK}" \
  -s Snakefile.py \
  --configfile "${CONFIG_REL}" \
  --config enable_chips_targets=true \
  --latency-wait 60 \
  --rerun-incomplete 2>&1 | tee -a "${LOG_FILE}"
exit_code=${PIPESTATUS[0]}
set -e

if [[ ${exit_code} -eq 0 ]]; then
  touch "${DEST_ROOT}/RUN_COMPLETE"
fi

echo "[$(date)] END ChIPs realstudy workflow exit=${exit_code}" | tee -a "${LOG_FILE}"
exit "${exit_code}"
