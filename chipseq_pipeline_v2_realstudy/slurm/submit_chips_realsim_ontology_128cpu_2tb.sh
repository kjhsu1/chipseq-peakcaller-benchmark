#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${REPO_ROOT:-$(cd "${SCRIPT_DIR}/.." && pwd)}"
DATA_ROOT="${DATA_ROOT:-/quobyte/ikorfgrp/home/kjhsu/results/chipseq_pipeline_v2_realstudy}"
STAGE_ROOT="${STAGE_ROOT:-${DATA_ROOT}/submission_sources}"
SBATCH_SCRIPT="${SBATCH_SCRIPT:-slurm/chips_realsim_ontology_128cpu_2tb.sbatch}"
RUN_TS="$(date +%Y%m%d_%H%M%S)"
STAGED_REPO="${STAGE_ROOT}/chips_realsim_ontology_submit_${RUN_TS}"
STAGED_SBATCH_SCRIPT="${STAGED_REPO}/${SBATCH_SCRIPT}"
EXPORT_VARS="ALL,PIPE_SRC=${STAGED_REPO}"

cd "${REPO_ROOT}"

mkdir -p "${STAGED_REPO}"
cp -a "${REPO_ROOT}/." "${STAGED_REPO}/"

job_id=$(sbatch --parsable --chdir "${STAGED_REPO}" --export "${EXPORT_VARS}" "${STAGED_SBATCH_SCRIPT}")
echo "Submitted ChIPs realstudy ontology workflow as ${job_id}"
