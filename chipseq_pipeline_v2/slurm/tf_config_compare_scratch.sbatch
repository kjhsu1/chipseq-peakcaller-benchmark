#!/bin/bash
#SBATCH --job-name=tf_cfg_cmp
#SBATCH --account=ikorfgrp
#SBATCH --partition=high
#SBATCH --output=/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/slurm_jobs/%j/%x.out
#SBATCH --error=/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/slurm_jobs/%j/%x.err
#SBATCH --chdir=/home/kjhsu/Code/Korflab/chipseq-peakcaller-benchmark/chipseq_pipeline_v2
#SBATCH --cpus-per-task=4
#SBATCH --mem=15G
#SBATCH --time=12:00:00

set -euo pipefail

CONDA_BIN="/cvmfs/hpc.ucdavis.edu/sw/conda/root/bin/conda"
ENV_PREFIX="/home/kjhsu/.conda/envs/background_project"
PIPE_SRC="/home/kjhsu/Code/Korflab/chipseq-peakcaller-benchmark/chipseq_pipeline_v2"
DATA_ROOT="/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2"
LOG_ROOT="${DATA_ROOT}/archived_results/sweep_logs"
RUN_NAME="${RUN_NAME:?RUN_NAME is required}"
BENCH_TAG="${BENCH_TAG:?BENCH_TAG is required}"
LOCAL_ROOT="/scratch/kjhsu_chipseq_pipeline_v2_compare/${SLURM_JOB_ID}"
PIPE_RUN="${LOCAL_ROOT}/pipeline"
RUN_RESULTS="${LOCAL_ROOT}/results"
CONFIG_REL="configs/${RUN_NAME}.yaml"
CONFIG_RUN="${PIPE_RUN}/${CONFIG_REL}"
TS="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${LOG_ROOT}/${RUN_NAME}_${BENCH_TAG}_${TS}.log"
SUMMARY_FILE="${LOG_ROOT}/${RUN_NAME}_${BENCH_TAG}_${TS}.summary.txt"

mkdir -p "${DATA_ROOT}/slurm_jobs/${SLURM_JOB_ID}" "${LOG_ROOT}" "${LOCAL_ROOT}" "${RUN_RESULTS}"
mkdir -p "/tmp/kjhsu_snakemake_cache/${SLURM_JOB_ID}"

echo "[$(date)] START ${RUN_NAME} ${BENCH_TAG}" | tee -a "${LOG_FILE}"
echo "local_root=${LOCAL_ROOT}" | tee -a "${LOG_FILE}"
echo "cpus=${SLURM_CPUS_PER_TASK}" | tee -a "${LOG_FILE}"

if [[ ! -x "${CONDA_BIN}" ]]; then
  echo "Missing conda binary: ${CONDA_BIN}" | tee -a "${LOG_FILE}" >&2
  exit 1
fi

if [[ ! -d "${ENV_PREFIX}" ]]; then
  echo "Missing conda env: ${ENV_PREFIX}" | tee -a "${LOG_FILE}" >&2
  exit 1
fi

if ! compgen -G "${PIPE_SRC}/data/indexes/ce11_1pct/bowtie2_index/ce11_1pct*.bt2" > /dev/null; then
  echo "Missing Bowtie2 index files under ${PIPE_SRC}/data/indexes/ce11_1pct/bowtie2_index" | tee -a "${LOG_FILE}" >&2
  exit 1
fi

mkdir -p "${PIPE_RUN}"
cp -a "${PIPE_SRC}/." "${PIPE_RUN}/"
find "${PIPE_RUN}" \( -type d -name '__pycache__' -o -type f -name '*.pyc' \) -exec rm -rf {} +
rm -rf "${PIPE_RUN}/results" "${PIPE_RUN}/archived_results" "${PIPE_RUN}/.snakemake"
ln -sfn "${RUN_RESULTS}" "${PIPE_RUN}/results"
ln -sfn "${DATA_ROOT}/archived_results" "${PIPE_RUN}/archived_results"

export XDG_CACHE_HOME="/tmp/kjhsu_snakemake_cache/${SLURM_JOB_ID}"
export PYTHONDONTWRITEBYTECODE=1

cd "${PIPE_RUN}"
start_epoch="$(date +%s)"

set +e
"${CONDA_BIN}" run -p "${ENV_PREFIX}" snakemake \
  -j "${SLURM_CPUS_PER_TASK}" \
  -s Snakefile.py \
  --configfile "${CONFIG_RUN}" \
  --latency-wait 60 \
  --rerun-incomplete 2>&1 | tee -a "${LOG_FILE}"
exit_code=${PIPESTATUS[0]}
set -e

end_epoch="$(date +%s)"
elapsed_seconds="$((end_epoch - start_epoch))"
bam_count="$(find "${RUN_RESULTS}" -type f -name '*.bam' | wc -l)"
peak_count="$(find "${RUN_RESULTS}" -type f -name '*_peaks.narrowPeak' | wc -l)"

{
  echo "run_name=${RUN_NAME}"
  echo "bench_tag=${BENCH_TAG}"
  echo "cpus=${SLURM_CPUS_PER_TASK}"
  echo "exit_code=${exit_code}"
  echo "elapsed_seconds=${elapsed_seconds}"
  echo "bam_count=${bam_count}"
  echo "peak_count=${peak_count}"
} | tee "${SUMMARY_FILE}"

rm -rf "${LOCAL_ROOT}"

echo "[$(date)] END ${RUN_NAME} ${BENCH_TAG} exit=${exit_code} elapsed_seconds=${elapsed_seconds}" | tee -a "${LOG_FILE}"

exit "${exit_code}"
