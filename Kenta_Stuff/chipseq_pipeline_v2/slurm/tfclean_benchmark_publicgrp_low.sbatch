#!/bin/bash
#SBATCH --job-name=tfclean_bench
#SBATCH --account=publicgrp
#SBATCH --partition=low
#SBATCH --output=/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/slurm_jobs/%j/%x.out
#SBATCH --error=/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/slurm_jobs/%j/%x.err
#SBATCH --chdir=/home/kjhsu/Code/Korflab/chipseq-peakcaller-benchmark/Kenta_Stuff/chipseq_pipeline_v2
#SBATCH --mem=15G
#SBATCH --time=1-00:00:00

set -euo pipefail

CONDA_BIN="/cvmfs/hpc.ucdavis.edu/sw/conda/root/bin/conda"
ENV_PREFIX="/home/kjhsu/.conda/envs/background_project"
PIPE_SRC="/home/kjhsu/Code/Korflab/chipseq-peakcaller-benchmark/Kenta_Stuff/chipseq_pipeline_v2"
DATA_ROOT="/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2"
LOG_ROOT="${DATA_ROOT}/archived_results/sweep_logs"
RUN_TS="$(date +%Y%m%d_%H%M%S)"
JOB_ROOT="${DATA_ROOT}/work/${SLURM_JOB_ID}"
RUN_NAME="${RUN_NAME:?RUN_NAME is required}"
BENCH_TAG="${BENCH_TAG:?BENCH_TAG is required}"

PIPE_RUN="${JOB_ROOT}/${RUN_NAME}/pipeline"
RUN_RESULTS="${JOB_ROOT}/${RUN_NAME}/results"
DEST_ROOT="${DATA_ROOT}/archived_results/${RUN_NAME}_${BENCH_TAG}_${RUN_TS}"
CONFIG_REL="configs/${RUN_NAME}.yaml"
CONFIG_RUN="${PIPE_RUN}/${CONFIG_REL}"
LOG_FILE="${LOG_ROOT}/${RUN_NAME}_${BENCH_TAG}_${RUN_TS}.log"

mkdir -p "${DATA_ROOT}/slurm_jobs/${SLURM_JOB_ID}" "${JOB_ROOT}" "${LOG_ROOT}" "${RUN_RESULTS}"
mkdir -p "/tmp/kjhsu_snakemake_cache/${SLURM_JOB_ID}"

echo "[$(date)] START ${RUN_NAME} ${BENCH_TAG} cpus=${SLURM_CPUS_PER_TASK}" | tee -a "${LOG_FILE}"
echo "run_results=${RUN_RESULTS}" | tee -a "${LOG_FILE}"
echo "dest_root=${DEST_ROOT}" | tee -a "${LOG_FILE}"

mkdir -p "${PIPE_RUN}"
cp -a "${PIPE_SRC}/." "${PIPE_RUN}/"
find "${PIPE_RUN}" \( -type d -name '__pycache__' -o -type f -name '*.pyc' \) -exec rm -rf {} +
rm -rf "${PIPE_RUN}/results" "${PIPE_RUN}/archived_results" "${PIPE_RUN}/.snakemake"
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
  --rerun-incomplete 2>&1 | tee -a "${LOG_FILE}"
exit_code=${PIPESTATUS[0]}
set -e

if [[ ${exit_code} -eq 0 ]]; then
  mv "${RUN_RESULTS}" "${DEST_ROOT}"
  touch "${DEST_ROOT}/RUN_COMPLETE"
  echo "archived -> ${DEST_ROOT}" | tee -a "${LOG_FILE}"
else
  echo "partial results remain at ${RUN_RESULTS}" | tee -a "${LOG_FILE}"
fi

echo "[$(date)] END ${RUN_NAME} ${BENCH_TAG} exit=${exit_code}" | tee -a "${LOG_FILE}"

exit "${exit_code}"
