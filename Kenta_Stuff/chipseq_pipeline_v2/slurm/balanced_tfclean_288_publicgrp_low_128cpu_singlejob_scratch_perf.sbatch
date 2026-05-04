#!/bin/bash
#SBATCH --job-name=tfclean288_publo_perf
#SBATCH --account=publicgrp
#SBATCH --partition=low
#SBATCH --output=/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/slurm_jobs/%j/%x.out
#SBATCH --error=/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2/slurm_jobs/%j/%x.err
#SBATCH --chdir=/home/kjhsu/Code/Korflab/chipseq-peakcaller-benchmark/Kenta_Stuff/chipseq_pipeline_v2
#SBATCH --cpus-per-task=128
#SBATCH --mem=15G
#SBATCH --time=4-00:00:00

set -euo pipefail

CONDA_BIN="/cvmfs/hpc.ucdavis.edu/sw/conda/root/bin/conda"
ENV_PREFIX="/home/kjhsu/.conda/envs/background_project"
PIPE_SRC="/home/kjhsu/Code/Korflab/chipseq-peakcaller-benchmark/Kenta_Stuff/chipseq_pipeline_v2"
DATA_ROOT="/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2"
LOCAL_ROOT="/scratch/kjhsu_chipseq_pipeline_v2/${SLURM_JOB_ID}"
LOG_ROOT="${DATA_ROOT}/archived_results/sweep_logs"
JOB_LOG="${LOG_ROOT}/balanced_tfclean_288_publicgrp_low_perf_${SLURM_JOB_ID}.log"

runs=(
  "balanced_tfclean_flatearth_peaks_broad_integrated_288"
  "balanced_tfclean_flatearth_plateaus_broad_integrated_288"
  "balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288"
  "balanced_tfclean_realistic_plateaus_wavy_broad_integrated_288"
  "balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288"
  "balanced_tfclean_realistic_plateaus_hilly_broad_integrated_288"
)

copy_pids=()

mkdir -p "${DATA_ROOT}/slurm_jobs/${SLURM_JOB_ID}" "${LOG_ROOT}" "${LOCAL_ROOT}"
mkdir -p "/tmp/kjhsu_snakemake_cache/${SLURM_JOB_ID}"

echo "[$(date)] START publicgrp low tfclean perf run" | tee -a "${JOB_LOG}"
echo "local_root=${LOCAL_ROOT}" | tee -a "${JOB_LOG}"
echo "cpus=${SLURM_CPUS_PER_TASK}" | tee -a "${JOB_LOG}"

export XDG_CACHE_HOME="/tmp/kjhsu_snakemake_cache/${SLURM_JOB_ID}"
export PYTHONDONTWRITEBYTECODE=1

for RUN_NAME in "${runs[@]}"; do
  TS="$(date +%Y%m%d_%H%M%S)"
  CONFIG_REL="configs/${RUN_NAME}.yaml"
  PIPE_RUN="${LOCAL_ROOT}/${RUN_NAME}/pipeline"
  RUN_RESULTS="${LOCAL_ROOT}/${RUN_NAME}/results"
  CONFIG_RUN="${PIPE_RUN}/${CONFIG_REL}"
  LOG_FILE="${LOG_ROOT}/${RUN_NAME}_publicgrp_low_perf_${TS}.log"
  DEST="${DATA_ROOT}/archived_results/${RUN_NAME}_${TS}"

  echo "[$(date)] START ${RUN_NAME}" | tee -a "${JOB_LOG}"
  mkdir -p "${PIPE_RUN}" "${RUN_RESULTS}"
  cp -a "${PIPE_SRC}/." "${PIPE_RUN}/"
  find "${PIPE_RUN}" \( -type d -name '__pycache__' -o -type f -name '*.pyc' \) -exec rm -rf {} +
  rm -rf "${PIPE_RUN}/results" "${PIPE_RUN}/archived_results" "${PIPE_RUN}/.snakemake"
  ln -sfn "${RUN_RESULTS}" "${PIPE_RUN}/results"
  ln -sfn "${DATA_ROOT}/archived_results" "${PIPE_RUN}/archived_results"

  cd "${PIPE_RUN}"

  set +e
  "${CONDA_BIN}" run -p "${ENV_PREFIX}" snakemake \
    -j "${SLURM_CPUS_PER_TASK}" \
    -s Snakefile.py \
    --configfile "${CONFIG_RUN}" \
    --config retain_bams=true emit_pmf_disabled_marker=false \
    --latency-wait 60 \
    --rerun-incomplete 2>&1 | tee "${LOG_FILE}"
  exit_code=${PIPESTATUS[0]}
  set -e

  if [[ ${exit_code} -ne 0 ]]; then
    echo "[$(date)] FAIL ${RUN_NAME} exit=${exit_code}" | tee -a "${JOB_LOG}"
    exit "${exit_code}"
  fi

  (
    cp -a "${RUN_RESULTS}" "${DEST}"
    touch "${DEST}/RUN_COMPLETE"
    echo "[$(date)] archived ${RUN_NAME} -> ${DEST}" | tee -a "${JOB_LOG}"
  ) &
  copy_pids+=("$!")

  echo "[$(date)] DONE ${RUN_NAME}; archive pid=${copy_pids[-1]}" | tee -a "${JOB_LOG}"
done

for pid in "${copy_pids[@]}"; do
  wait "${pid}"
done

echo "[$(date)] END publicgrp low tfclean perf run" | tee -a "${JOB_LOG}"
