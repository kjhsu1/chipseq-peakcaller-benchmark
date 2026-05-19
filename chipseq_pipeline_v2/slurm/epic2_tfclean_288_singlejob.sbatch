#!/bin/bash
#SBATCH --job-name=epic2_tfclean_288
#SBATCH --account=publicgrp
#SBATCH --partition=low
#SBATCH --output=/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/slurm_jobs/%j/%x.out
#SBATCH --error=/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/slurm_jobs/%j/%x.err
#SBATCH --chdir=/home/kjhsu/Code/Korflab/chipseq-peakcaller-benchmark/chipseq_pipeline_v2
#SBATCH --cpus-per-task=64
#SBATCH --mem=50G
#SBATCH --time=2-00:00:00

set -euo pipefail

CONDA_BIN="${CONDA_BIN:-/cvmfs/hpc.ucdavis.edu/sw/conda/root/bin/conda}"
ENV_PREFIX="${ENV_PREFIX:-/home/kjhsu/.conda/envs/background_project}"
PIPE_SRC="${PIPE_SRC:-/home/kjhsu/Code/Korflab/chipseq-peakcaller-benchmark/chipseq_pipeline_v2}"
DATA_ROOT="${DATA_ROOT:-/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2}"
RUN_NAME="${RUN_NAME:?Set RUN_NAME to one epic2_tfclean config basename without .yaml}"

RUN_TS="$(date +%Y%m%d_%H%M%S)"
JOB_ROOT="${DATA_ROOT}/work/${SLURM_JOB_ID}/${RUN_NAME}/pipeline"
DEST_ROOT="${DATA_ROOT}/archived_results/${RUN_NAME}_${RUN_TS}"
LOG_ROOT="${DATA_ROOT}/archived_results/sweep_logs"
LOG_FILE="${LOG_ROOT}/${RUN_NAME}_epic2_${RUN_TS}.log"
CONFIG_REL="configs/${RUN_NAME}.yaml"
CONFIG_RUN="${JOB_ROOT}/${CONFIG_REL}"

mkdir -p "${DATA_ROOT}/slurm_jobs/${SLURM_JOB_ID}" "${JOB_ROOT}" "${DEST_ROOT}" "${LOG_ROOT}"
mkdir -p "/tmp/kjhsu_snakemake_cache/${SLURM_JOB_ID}"

echo "[$(date)] START ${RUN_NAME}" | tee -a "${LOG_FILE}"
echo "results=${DEST_ROOT}" | tee -a "${LOG_FILE}"

export XDG_CACHE_HOME="/tmp/kjhsu_snakemake_cache/${SLURM_JOB_ID}"
export PYTHONDONTWRITEBYTECODE=1

cp -a "${PIPE_SRC}/." "${JOB_ROOT}/"
find "${JOB_ROOT}" \( -type d -name '__pycache__' -o -type f -name '*.pyc' \) -exec rm -rf {} +
rm -rf "${JOB_ROOT}/results" "${JOB_ROOT}/archived_results" "${JOB_ROOT}/.snakemake"
ln -sfn "${DEST_ROOT}" "${JOB_ROOT}/results"
ln -sfn "${DATA_ROOT}/archived_results" "${JOB_ROOT}/archived_results"

cd "${JOB_ROOT}"

set +e
"${CONDA_BIN}" run -p "${ENV_PREFIX}" snakemake \
  -j "${SLURM_CPUS_PER_TASK}" \
  -s Snakefile.py \
  --configfile "${CONFIG_RUN}" \
  --latency-wait 60 \
  --rerun-incomplete 2>&1 | tee -a "${LOG_FILE}"
exit_code=${PIPESTATUS[0]}
set -e

if [[ ${exit_code} -eq 0 ]]; then
  touch "${DEST_ROOT}/RUN_COMPLETE"
fi

echo "[$(date)] END ${RUN_NAME} exit=${exit_code}" | tee -a "${LOG_FILE}"
exit "${exit_code}"
