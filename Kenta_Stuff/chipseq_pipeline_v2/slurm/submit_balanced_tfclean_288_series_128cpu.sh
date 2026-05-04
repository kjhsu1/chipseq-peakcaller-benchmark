#!/bin/bash

set -euo pipefail

REPO_ROOT="/home/kjhsu/Code/Korflab/chipseq-peakcaller-benchmark/Kenta_Stuff/chipseq_pipeline_v2"
SBATCH_SCRIPT="slurm/balanced_tfclean_288_128cpu.sbatch"

runs=(
  "balanced_tfclean_flatearth_peaks_broad_integrated_288"
  "balanced_tfclean_flatearth_plateaus_broad_integrated_288"
  "balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288"
  "balanced_tfclean_realistic_plateaus_wavy_broad_integrated_288"
  "balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288"
  "balanced_tfclean_realistic_plateaus_hilly_broad_integrated_288"
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
