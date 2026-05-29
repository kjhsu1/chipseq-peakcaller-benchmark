#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"

eval "$(/opt/anaconda3/bin/conda shell.zsh hook)"
conda activate background_project

export MPLCONFIGDIR="${MPLCONFIGDIR:-/private/tmp/mplconfig_homer_stage}"
mkdir -p "$MPLCONFIGDIR"

cd "$ROOT_DIR"

python scripts/build_homer_staged_report.py \
  --output-dir analysis_outputs/homer_tfclean_128_staged_report_current

python scripts/summarize_homer_staged_categories.py \
  --output-dir analysis_outputs/homer_staged_summary_current

python scripts/compare_peakcaller_staged_categories.py \
  --output-dir analysis_outputs/peakcaller_staged_comparison_current
