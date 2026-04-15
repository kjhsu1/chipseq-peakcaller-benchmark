#!/bin/bash
#SBATCH --job-name=real_plateau128
#SBATCH --account=ikorfgrp
#SBATCH --partition=high
#SBATCH --output=/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/slurm_jobs/%j/%x.out
#SBATCH --error=/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/slurm_jobs/%j/%x.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=3:00:00

set -euo pipefail

: """Paths"""

CONDA_BIN="/cvmfs/hpc.ucdavis.edu/sw/conda/root/bin/conda"
ENV_PREFIX="/home/kjhsu/.conda/envs/background_project"
PIPE="/home/kjhsu/Code/Korflab/chipseq-peakcaller-benchmark/Kenta_Stuff/chipseq_pipeline_v2"
DATA="/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2"
CONFIG_PATH="${PIPE}/configs/realistic_plateau_integrated_enrich128.yaml"


: """Setup"""

mkdir -p "${DATA}/slurm_jobs/${SLURM_JOB_ID}"
mkdir -p "${DATA}/results" "${DATA}/archived_results/sweep_logs"
mkdir -p /tmp/kjhsu_snakemake_cache

cd "${PIPE}"

rm -rf results archived_results
ln -sfn "${DATA}/results" results
ln -sfn "${DATA}/archived_results" archived_results

rm -rf "${DATA}/results"/*
ts=$(date +%Y%m%d_%H%M%S)
log="${DATA}/archived_results/sweep_logs/realistic_plateau_integrated_enrich128_${ts}.log"
export XDG_CACHE_HOME=/tmp/kjhsu_snakemake_cache

echo "[$(date)] START ${CONFIG_PATH}" | tee "${log}"


: """Run"""

set +e
"${CONDA_BIN}" run -p "${ENV_PREFIX}" snakemake \
  -j "${SLURM_CPUS_PER_TASK}" \
  -s Snakefile.py \
  --configfile "${CONFIG_PATH}" \
  --latency-wait 60 2>&1 | tee -a "${log}"
exit_code=${PIPESTATUS[0]}
set -e

if [[ ${exit_code} -eq 0 ]]; then
  dest="${DATA}/archived_results/realistic_plateau_integrated_enrich128_${ts}"
  cp -a "${DATA}/results" "${dest}"
  echo "archived -> ${dest}" | tee -a "${log}"
  rm -rf "${DATA}/results"/*
fi

echo "[$(date)] END ${CONFIG_PATH} exit=${exit_code}" | tee -a "${log}"

exit "${exit_code}"
