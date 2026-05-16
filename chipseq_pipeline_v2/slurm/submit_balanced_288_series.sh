#!/bin/bash

set -euo pipefail

REPO_ROOT="/home/kjhsu/Code/Korflab/chipseq-peakcaller-benchmark/chipseq_pipeline_v2"

scripts=(
  "slurm/balanced_flatearth_peaks_broad_integrated_288.sbatch"
  "slurm/balanced_flatearth_plateaus_broad_integrated_288.sbatch"
  "slurm/balanced_realistic_plateaus_wavy_broad_integrated_288.sbatch"
  "slurm/balanced_realistic_peaks_wavy_narrow_integrated_288.sbatch"
  "slurm/balanced_realistic_plateaus_hilly_broad_integrated_288.sbatch"
  "slurm/balanced_realistic_peaks_hilly_narrow_integrated_288.sbatch"
)

cd "${REPO_ROOT}"

prev_jid=""
for script in "${scripts[@]}"; do
  if [[ -z "${prev_jid}" ]]; then
    jid=$(sbatch --parsable "${script}")
  else
    jid=$(sbatch --dependency=afterok:${prev_jid} --parsable "${script}")
  fi
  echo "Submitted ${script} as ${jid}"
  prev_jid="${jid}"
done
