#!/usr/bin/env bash
set -euo pipefail

PIPELINE_DIR="${PIPELINE_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
OUTPUT_ROOT="${OUTPUT_ROOT:?Set OUTPUT_ROOT before submission}"
export PIPELINE_DIR OUTPUT_ROOT

mkdir -p "${PIPELINE_DIR}/slurm_logs"
sbatch --export=ALL "${PIPELINE_DIR}/slurm/realstudy_v2_singlejob.sbatch"
