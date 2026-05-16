#!/bin/bash
#SBATCH --job-name=tfclean_hilly_rem
#SBATCH --account=publicgrp
#SBATCH --partition=low
#SBATCH --output=/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/slurm_jobs/%j/%x.out
#SBATCH --error=/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/slurm_jobs/%j/%x.err
#SBATCH --chdir=/home/kjhsu/Code/Korflab/chipseq-peakcaller-benchmark/chipseq_pipeline_v2
#SBATCH --cpus-per-task=128
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00

set -euo pipefail

CONDA_BIN="/cvmfs/hpc.ucdavis.edu/sw/conda/root/bin/conda"
ENV_PREFIX="/home/kjhsu/.conda/envs/background_project"
PIPE_SRC="/home/kjhsu/Code/Korflab/chipseq-peakcaller-benchmark/chipseq_pipeline_v2"
DATA_ROOT="/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2"
JOB_ROOT="${DATA_ROOT}/work/${SLURM_JOB_ID}"
PIPE_RUN="${JOB_ROOT}/pipeline"
LOG_ROOT="${DATA_ROOT}/archived_results/sweep_logs"
RUN_TS="$(date +%Y%m%d_%H%M%S)"

runs=(
  "balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288"
  "balanced_tfclean_realistic_plateaus_hilly_broad_integrated_288"
)

mkdir -p "${DATA_ROOT}/slurm_jobs/${SLURM_JOB_ID}" "${JOB_ROOT}" "${LOG_ROOT}"
mkdir -p "${PIPE_RUN}"
cp -a "${PIPE_SRC}/." "${PIPE_RUN}/"
find "${PIPE_RUN}" \( -type d -name '__pycache__' -o -type f -name '*.pyc' \) -exec rm -rf {} +
rm -rf "${PIPE_RUN}/results" "${PIPE_RUN}/archived_results" "${PIPE_RUN}/.snakemake"
ln -sfn "${DATA_ROOT}/archived_results" "${PIPE_RUN}/archived_results"

export PYTHONDONTWRITEBYTECODE=1
export XDG_CACHE_HOME="/tmp/kjhsu_snakemake_cache/${SLURM_JOB_ID}"
mkdir -p "${XDG_CACHE_HOME}"

cd "${PIPE_RUN}"

echo "[$(date)] START hilly tfclean remainder cpus=${SLURM_CPUS_PER_TASK} mem=100G"

for RUN_NAME in "${runs[@]}"; do
  DEST="${DATA_ROOT}/archived_results/${RUN_NAME}_${RUN_TS}"
  LOG_FILE="${LOG_ROOT}/${RUN_NAME}_hilly_remainder_${RUN_TS}.log"

  echo "[$(date)] START ${RUN_NAME}"
  rm -rf results
  mkdir -p "${DEST}"
  ln -sfn "${DEST}" results

  "${CONDA_BIN}" run -p "${ENV_PREFIX}" snakemake \
    -j "${SLURM_CPUS_PER_TASK}" \
    -s Snakefile.py \
    --configfile "configs/${RUN_NAME}.yaml" \
    --config retain_bams=true emit_pmf_disabled_marker=false \
    --latency-wait 60 \
    --rerun-incomplete 2>&1 | tee "${LOG_FILE}"

  touch "${DEST}/RUN_COMPLETE"
  echo "[$(date)] DONE ${RUN_NAME} -> ${DEST}"
done

echo "[$(date)] END hilly tfclean remainder"
