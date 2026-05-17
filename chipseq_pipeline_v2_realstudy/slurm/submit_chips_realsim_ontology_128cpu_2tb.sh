#!/bin/bash

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/home/kjhsu/Code/Korflab/chipseq-peakcaller-benchmark/chipseq_pipeline_v2_realstudy}"
SBATCH_SCRIPT="${SBATCH_SCRIPT:-slurm/chips_realsim_ontology_128cpu_2tb.sbatch}"

cd "${REPO_ROOT}"

job_id=$(sbatch --parsable "${SBATCH_SCRIPT}")
echo "Submitted ChIPs realstudy ontology 128-core/2TB workflow as ${job_id}"
