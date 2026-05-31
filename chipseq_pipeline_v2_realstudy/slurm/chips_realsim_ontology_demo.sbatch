#!/bin/bash
#SBATCH --job-name=chips_ontology_demo
#SBATCH --account=publicgrp
#SBATCH --partition=low
#SBATCH --output=/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2_realstudy/slurm_jobs/%j/%x.out
#SBATCH --error=/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2_realstudy/slurm_jobs/%j/%x.err
#SBATCH --chdir=/home/kjhsu/Code/Korflab/chipseq-peakcaller-benchmark/chipseq_pipeline_v2_realstudy
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=0-08:00:00

set -euo pipefail

CONDA_BIN="${CONDA_BIN:-/cvmfs/hpc.ucdavis.edu/sw/conda/root/bin/conda}"
ENV_PREFIX="${ENV_PREFIX:-/home/kjhsu/.conda/envs/background_project}"
PIPE_SRC="${PIPE_SRC:-${SLURM_SUBMIT_DIR}}"
DATA_ROOT="${DATA_ROOT:-/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2_realstudy}"
CONFIG_FILES="${CONFIG_FILES:-config.yaml configs/chips_demo_pairing.yaml}"

RUN_TS="$(date +%Y%m%d_%H%M%S)"
JOB_ROOT="${DATA_ROOT}/work/${SLURM_JOB_ID}/pipeline"
DEST_ROOT="${DATA_ROOT}/chips_runs/chips_realsim_ontology_demo_${RUN_TS}"
LOG_ROOT="${DATA_ROOT}/logs"
LOG_FILE="${LOG_ROOT}/chips_realsim_ontology_demo_${RUN_TS}.log"
SUMMARY_ROOT="${DEST_ROOT}/analysis_outputs/chips_ontology/summary"
REPRO_ROOT="${DEST_ROOT}/reproducibility"
VALIDATION_REPORT="${REPRO_ROOT}/validation_report.md"
VALIDATION_JSON="${REPRO_ROOT}/validation_summary.json"
RESOURCE_TXT="${SUMMARY_ROOT}/runtime_resource_report.txt"
RESOURCE_CSV="${SUMMARY_ROOT}/runtime_resource_report.csv"

mkdir -p "${DATA_ROOT}/slurm_jobs/${SLURM_JOB_ID}" "${JOB_ROOT}" "${DEST_ROOT}" "${LOG_ROOT}"
mkdir -p "${DEST_ROOT}/results_chips" "${DEST_ROOT}/analysis_outputs/chips_ontology"
mkdir -p "/tmp/kjhsu_snakemake_cache/${SLURM_JOB_ID}"

echo "[$(date)] START ChIPs realstudy ontology demo workflow" | tee -a "${LOG_FILE}"
echo "results=${DEST_ROOT}" | tee -a "${LOG_FILE}"
echo "cpus=${SLURM_CPUS_PER_TASK}" | tee -a "${LOG_FILE}"
echo "mem=${SLURM_MEM_PER_NODE:-64G}" | tee -a "${LOG_FILE}"

export XDG_CACHE_HOME="/tmp/kjhsu_snakemake_cache/${SLURM_JOB_ID}"
export PYTHONDONTWRITEBYTECODE=1

cp -a "${PIPE_SRC}/." "${JOB_ROOT}/"
find "${JOB_ROOT}" \( -type d -name '__pycache__' -o -type f -name '*.pyc' \) -exec rm -rf {} +
rm -rf "${JOB_ROOT}/results_chips" "${JOB_ROOT}/analysis_outputs/chips_ontology" "${JOB_ROOT}/.snakemake"
ln -sfn "${DEST_ROOT}/results_chips" "${JOB_ROOT}/results_chips"
mkdir -p "${JOB_ROOT}/analysis_outputs"
ln -sfn "${DEST_ROOT}/analysis_outputs/chips_ontology" "${JOB_ROOT}/analysis_outputs/chips_ontology"

cd "${JOB_ROOT}"

read -r -a CONFIG_FILE_ARRAY <<< "${CONFIG_FILES}"
SNAKEMAKE_CONFIG_ARGS=()
for config_file in "${CONFIG_FILE_ARRAY[@]}"; do
  SNAKEMAKE_CONFIG_ARGS+=(--configfile "${config_file}")
done

write_resource_reports() {
  mkdir -p "${SUMMARY_ROOT}" "${REPRO_ROOT}"
  {
    echo "slurm_job_id=${SLURM_JOB_ID}"
    echo "job_name=${SLURM_JOB_NAME:-chips_ontology_demo}"
    echo "partition=${SLURM_JOB_PARTITION:-${SLURM_PARTITION:-UNAVAILABLE}}"
    echo "account=${SLURM_JOB_ACCOUNT:-${SLURM_ACCOUNT:-UNAVAILABLE}}"
    echo "cpu_allocation=${SLURM_CPUS_PER_TASK:-UNAVAILABLE}"
    echo "requested_mem=${SLURM_MEM_PER_NODE:-UNAVAILABLE}"
    echo "output_root=${DEST_ROOT}"
    echo "log_file=${LOG_FILE}"
    echo "snakemake_exit_code=${exit_code}"
    if [[ ${exit_code} -eq 0 ]]; then
      echo "completed_job_count=1"
      echo "failed_job_count=0"
    else
      echo "completed_job_count=0"
      echo "failed_job_count=1"
    fi
  } > "${RESOURCE_TXT}"

  if sacct -j "${SLURM_JOB_ID}" --parsable2 \
    --format=JobIDRaw,JobName,Partition,Account,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,MaxDiskRead,MaxDiskWrite,Submit,Start,End \
    > "${RESOURCE_CSV}.tmp" 2>>"${LOG_FILE}"; then
    mv "${RESOURCE_CSV}.tmp" "${RESOURCE_CSV}"
    {
      echo ""
      echo "sacct_csv=${RESOURCE_CSV}"
      echo "sacct_available=true"
      tail -n +1 "${RESOURCE_CSV}"
    } >> "${RESOURCE_TXT}"
  else
    rm -f "${RESOURCE_CSV}.tmp"
    printf '%s\n' "metric,value" "sacct,UNAVAILABLE" > "${RESOURCE_CSV}"
    {
      echo ""
      echo "sacct_available=false"
      echo "sacct_reason=UNAVAILABLE"
    } >> "${RESOURCE_TXT}"
  fi
  cp "${RESOURCE_TXT}" "${REPRO_ROOT}/runtime_resource_report.txt"
  cp "${RESOURCE_CSV}" "${REPRO_ROOT}/runtime_resource_report.csv"
}

package_reproducibility() {
  mkdir -p "${REPRO_ROOT}/configs" "${REPRO_ROOT}/metadata" "${REPRO_ROOT}/slurm"
  cp "config.yaml" "${REPRO_ROOT}/config.yaml"
  cp "configs/chips_demo_pairing.yaml" "${REPRO_ROOT}/configs/chips_demo_pairing.yaml"
  cp "metadata/chips_demo_pairing_run_table.csv" "${REPRO_ROOT}/metadata/chips_demo_pairing_run_table.csv"
  cp "metadata/prototype_run_table.summary.json" "${REPRO_ROOT}/metadata/prototype_run_table.summary.json"
  cp "slurm/submit_chips_realsim_ontology_demo.sh" "${REPRO_ROOT}/slurm/submit_chips_realsim_ontology_demo.sh"
  cp "slurm/chips_realsim_ontology_demo.sbatch" "${REPRO_ROOT}/slurm/chips_realsim_ontology_demo.sbatch"
  if git rev-parse HEAD > "${REPRO_ROOT}/commit_hash.txt" 2>>"${LOG_FILE}"; then
    :
  else
    printf '%s\n' "UNAVAILABLE" > "${REPRO_ROOT}/commit_hash.txt"
  fi
  printf '%s\n' "${DEST_ROOT}" > "${REPRO_ROOT}/output_root.txt"
  printf '%s\n' "${LOG_FILE}" > "${REPRO_ROOT}/log_file.txt"
  if [[ -f "${SUMMARY_ROOT}/failure_mode_summary.md" ]]; then
    cp "${SUMMARY_ROOT}/failure_mode_summary.md" "${REPRO_ROOT}/top_conclusions.md"
  else
    printf '%s\n' "# Top Conclusions" "" "UNAVAILABLE: summary file not generated." > "${REPRO_ROOT}/top_conclusions.md"
  fi
}

run_validation() {
  mkdir -p "${REPRO_ROOT}"
  if "${CONDA_BIN}" run -p "${ENV_PREFIX}" python scripts/validate_chips_ontology_production_run.py \
    --run-table metadata/chips_demo_pairing_run_table.csv \
    --output-root "${DEST_ROOT}" \
    --summary-dir "${SUMMARY_ROOT}" \
    --write-report "${VALIDATION_REPORT}" >> "${LOG_FILE}" 2>&1; then
    validation_exit=0
  else
    validation_exit=$?
  fi
  if [[ -f "${VALIDATION_REPORT}" && -f "${VALIDATION_JSON}" ]]; then
    :
  elif [[ -f "${VALIDATION_REPORT}" && -f "${REPRO_ROOT}/validation_summary.json" ]]; then
    mv "${REPRO_ROOT}/validation_summary.json" "${VALIDATION_JSON}"
  fi
}

set +e
"${CONDA_BIN}" run -p "${ENV_PREFIX}" snakemake \
  -j "${SLURM_CPUS_PER_TASK}" \
  -s Snakefile.py \
  "${SNAKEMAKE_CONFIG_ARGS[@]}" \
  --config enable_chips_targets=true enable_chips_ontology_targets=true \
  --latency-wait 60 \
  --rerun-incomplete 2>&1 | tee -a "${LOG_FILE}"
exit_code=${PIPESTATUS[0]}
set -e
resource_exit=0
validation_exit=0
package_exit=0

set +e
write_resource_reports || resource_exit=$?
package_reproducibility || package_exit=$?
run_validation
set -e

final_exit=${exit_code}
if [[ ${final_exit} -eq 0 && ${resource_exit} -ne 0 ]]; then
  final_exit=${resource_exit}
fi
if [[ ${final_exit} -eq 0 && ${package_exit} -ne 0 ]]; then
  final_exit=${package_exit}
fi
if [[ ${final_exit} -eq 0 && ${validation_exit} -ne 0 ]]; then
  final_exit=${validation_exit}
fi

if [[ ${final_exit} -eq 0 ]]; then
  touch "${DEST_ROOT}/RUN_COMPLETE"
  echo "archived -> ${DEST_ROOT}" | tee -a "${LOG_FILE}"
else
  touch "${DEST_ROOT}/RUN_FAILED"
fi

echo "validation_exit=${validation_exit}" | tee -a "${LOG_FILE}"
echo "resource_exit=${resource_exit}" | tee -a "${LOG_FILE}"
echo "package_exit=${package_exit}" | tee -a "${LOG_FILE}"
echo "[$(date)] END ChIPs realstudy ontology demo workflow exit=${final_exit}" | tee -a "${LOG_FILE}"
exit "${final_exit}"
