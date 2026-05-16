#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/.."

old_job=$(
  sbatch --parsable \
    --export=ALL,RUN_NAME=balanced_flatearth_peaks_broad_integrated_288,BENCH_TAG=oldcfg_scratch4c \
    slurm/tf_config_compare_scratch.sbatch
)
echo "old=${old_job}"

new_job=$(
  sbatch --parsable \
    --export=ALL,RUN_NAME=balanced_tfclean_flatearth_peaks_broad_integrated_288,BENCH_TAG=tfclean_scratch4c \
    slurm/tf_config_compare_scratch.sbatch
)
echo "new=${new_job}"
