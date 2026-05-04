#!/bin/bash

set -euo pipefail

REPO_ROOT="/home/kjhsu/Code/Korflab/chipseq-peakcaller-benchmark/Kenta_Stuff/chipseq_pipeline_v2"
SBATCH_SCRIPT="slurm/tfclean_benchmark_publicgrp_low.sbatch"
RUN_NAME="balanced_tfclean_flatearth_peaks_broad_integrated_288"

cd "${REPO_ROOT}"

for cpus in 8 16 32; do
  tag="bench_${cpus}cpu"
  jid=$(sbatch --parsable --cpus-per-task="${cpus}" --job-name="${RUN_NAME}_${tag}" --export=ALL,RUN_NAME="${RUN_NAME}",BENCH_TAG="${tag}" "${SBATCH_SCRIPT}")
  echo "Submitted ${RUN_NAME} ${tag} as ${jid}"
done
