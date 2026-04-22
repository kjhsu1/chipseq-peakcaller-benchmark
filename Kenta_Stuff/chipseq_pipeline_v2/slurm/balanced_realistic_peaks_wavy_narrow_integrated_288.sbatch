#!/bin/bash
#SBATCH --job-name=real_wavy_peak288
#SBATCH --account=ikorfgrp
#SBATCH --partition=high
#SBATCH --output=/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/slurm_jobs/%j/%x.out
#SBATCH --error=/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/slurm_jobs/%j/%x.err
#SBATCH --chdir=/home/kjhsu/Code/Korflab/chipseq-peakcaller-benchmark/Kenta_Stuff/chipseq_pipeline_v2
#SBATCH --cpus-per-task=4
#SBATCH --mem=15G
#SBATCH --time=12:00:00

set -euo pipefail

CONDA_BIN="/cvmfs/hpc.ucdavis.edu/sw/conda/root/bin/conda"
ENV_PREFIX="/home/kjhsu/.conda/envs/background_project"
PIPE_SRC="/home/kjhsu/Code/Korflab/chipseq-peakcaller-benchmark/Kenta_Stuff/chipseq_pipeline_v2"
DATA_ROOT="/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2"
RUN_NAME="balanced_realistic_peaks_wavy_narrow_integrated_288"
CONFIG_REL="configs/${RUN_NAME}.yaml"

JOB_ROOT="${DATA_ROOT}/work/${SLURM_JOB_ID}"
PIPE_RUN="${JOB_ROOT}/pipeline"
RUN_RESULTS="${JOB_ROOT}/results"
LOG_ROOT="${DATA_ROOT}/archived_results/sweep_logs"
TS="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${LOG_ROOT}/${RUN_NAME}_${TS}.log"
CONFIG_RUN="${PIPE_RUN}/${CONFIG_REL}"

mkdir -p "${DATA_ROOT}/slurm_jobs/${SLURM_JOB_ID}" "${JOB_ROOT}" "${LOG_ROOT}"
mkdir -p "/tmp/kjhsu_snakemake_cache/${SLURM_JOB_ID}"

echo "[$(date)] START ${RUN_NAME}"
echo "job_root=${JOB_ROOT}"
echo "pipe_src=${PIPE_SRC}"
echo "env_prefix=${ENV_PREFIX}"

if [[ ! -x "${CONDA_BIN}" ]]; then
  echo "Missing conda binary: ${CONDA_BIN}" >&2
  exit 1
fi

if [[ ! -d "${ENV_PREFIX}" ]]; then
  echo "Missing conda env: ${ENV_PREFIX}" >&2
  exit 1
fi

if ! compgen -G "${PIPE_SRC}/data/indexes/ce11_1pct/bowtie2_index/ce11_1pct*.bt2" > /dev/null; then
  echo "Missing Bowtie2 index files under ${PIPE_SRC}/data/indexes/ce11_1pct/bowtie2_index" >&2
  exit 1
fi

mkdir -p "${PIPE_RUN}"
cp -a "${PIPE_SRC}/." "${PIPE_RUN}/"
find "${PIPE_RUN}" \( -type d -name '__pycache__' -o -type f -name '*.pyc' \) -exec rm -rf {} +
rm -rf "${PIPE_RUN}/results" "${PIPE_RUN}/archived_results" "${PIPE_RUN}/.snakemake"
mkdir -p "${RUN_RESULTS}"
ln -sfn "${RUN_RESULTS}" "${PIPE_RUN}/results"
ln -sfn "${DATA_ROOT}/archived_results" "${PIPE_RUN}/archived_results"

export XDG_CACHE_HOME="/tmp/kjhsu_snakemake_cache/${SLURM_JOB_ID}"
export PYTHONDONTWRITEBYTECODE=1

cd "${PIPE_RUN}"

set +e
"${CONDA_BIN}" run -p "${ENV_PREFIX}" snakemake \
  -j "${SLURM_CPUS_PER_TASK}" \
  -s Snakefile.py \
  --configfile "${CONFIG_RUN}" \
  --latency-wait 60 \
  --rerun-incomplete 2>&1 | tee "${LOG_FILE}"
exit_code=${PIPESTATUS[0]}
set -e

if [[ ${exit_code} -eq 0 ]]; then
  dest="${DATA_ROOT}/archived_results/${RUN_NAME}_${TS}"
  cp -a "${RUN_RESULTS}" "${dest}"
  echo "archived -> ${dest}" | tee -a "${LOG_FILE}"
else
  echo "partial results remain at ${RUN_RESULTS}" | tee -a "${LOG_FILE}"
fi

echo "[$(date)] END ${RUN_NAME} exit=${exit_code}" | tee -a "${LOG_FILE}"

exit "${exit_code}"
