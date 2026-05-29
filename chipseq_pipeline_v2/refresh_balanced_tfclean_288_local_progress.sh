#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"

SELECTORS=("$@")
if [[ "${#SELECTORS[@]}" -eq 0 ]]; then
  SELECTORS=("all_six")
fi

eval "$(/opt/anaconda3/bin/conda shell.zsh hook)"
conda activate background_project

cd "$ROOT_DIR"

python scripts/summarize_balanced_tfclean_288_progress.py "${SELECTORS[@]}"
