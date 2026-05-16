#!/bin/bash
#SBATCH --job-name=flat_broad288
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
RUN_NAME="balanced_flatearth_peaks_broad_integrated_288"
CONFIG_PATH="configs/${RUN_NAME}.yaml"

JOB_ROOT="${DATA_ROOT}/work/${SLURM_JOB_ID}"
PIPE_RUN="${JOB_ROOT}/pipeline"
RUN_RESULTS="${JOB_ROOT}/results"
TS="$(date +%Y%m%d_%H%M%S)"

mkdir -p "${DATA_ROOT}/slurm_jobs/${SLURM_JOB_ID}" "${JOB_ROOT}" "${DATA_ROOT}/archived_results/sweep_logs"
cp -a "${PIPE_SRC}/." "${PIPE_RUN}"
find "${PIPE_RUN}" \( -type d -name '__pycache__' -o -type f -name '*.pyc' \) -exec rm -rf {} +
rm -rf "${PIPE_RUN}/results" "${PIPE_RUN}/archived_results" "${PIPE_RUN}/.snakemake"
mkdir -p "${RUN_RESULTS}"
ln -sfn "${RUN_RESULTS}" "${PIPE_RUN}/results"
ln -sfn "${DATA_ROOT}/archived_results" "${PIPE_RUN}/archived_results"

export XDG_CACHE_HOME="/tmp/kjhsu_snakemake_cache/${SLURM_JOB_ID}"
export PYTHONDONTWRITEBYTECODE=1

cd "${PIPE_RUN}"
"${CONDA_BIN}" run -p "${ENV_PREFIX}" snakemake -j "${SLURM_CPUS_PER_TASK}" -s Snakefile.py --configfile "${CONFIG_PATH}" --latency-wait 60 --rerun-incomplete 2>&1 | tee "${DATA_ROOT}/archived_results/sweep_logs/${RUN_NAME}_${TS}.log"
cp -a "${RUN_RESULTS}" "${DATA_ROOT}/archived_results/${RUN_NAME}_${TS}"
