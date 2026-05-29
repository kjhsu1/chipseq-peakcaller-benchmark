#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"

eval "$(/opt/anaconda3/bin/conda shell.zsh hook)"
conda activate background_project

cd "$ROOT_DIR"

python scripts/summarize_balanced_tfclean_288_progress.py all_six
python scripts/build_balanced_tfclean_288_decision_summary.py
