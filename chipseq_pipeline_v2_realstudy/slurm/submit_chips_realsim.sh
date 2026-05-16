#!/bin/bash

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/home/kjhsu/Code/Korflab/chipseq-peakcaller-benchmark/chipseq_pipeline_v2_realstudy}"
SBATCH_SCRIPT="${SBATCH_SCRIPT:-slurm/chips_realsim_singlejob.sbatch}"

cd "${REPO_ROOT}"

job_id=$(sbatch --parsable "${SBATCH_SCRIPT}")
echo "Submitted ChIPs realstudy workflow as ${job_id}"
