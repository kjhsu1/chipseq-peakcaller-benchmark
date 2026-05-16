#!/bin/bash

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/home/kjhsu/Code/Korflab/chipseq-peakcaller-benchmark/chipseq_pipeline_v2}"
SBATCH_SCRIPT="${SBATCH_SCRIPT:-slurm/epic2_tfclean_288_singlejob.sbatch}"

runs=(
  "epic2_tfclean_flatearth_peaks_broad_integrated_288"
  "epic2_tfclean_flatearth_plateaus_broad_integrated_288"
  "epic2_tfclean_realistic_peaks_wavy_narrow_integrated_288"
  "epic2_tfclean_realistic_plateaus_wavy_broad_integrated_288"
  "epic2_tfclean_realistic_peaks_hilly_narrow_integrated_288"
  "epic2_tfclean_realistic_plateaus_hilly_broad_integrated_288"
)

cd "${REPO_ROOT}"

prev_jid=""
for run_name in "${runs[@]}"; do
  if [[ -z "${prev_jid}" ]]; then
    jid=$(sbatch --parsable --job-name="${run_name}" --export=ALL,RUN_NAME="${run_name}" "${SBATCH_SCRIPT}")
  else
    jid=$(sbatch --dependency=afterok:${prev_jid} --parsable --job-name="${run_name}" --export=ALL,RUN_NAME="${run_name}" "${SBATCH_SCRIPT}")
  fi
  echo "Submitted ${run_name} as ${jid}"
  prev_jid="${jid}"
done

echo "Submitted EPIC2 TF-clean series: 6 configs x 288 runs = 1728 runs"
