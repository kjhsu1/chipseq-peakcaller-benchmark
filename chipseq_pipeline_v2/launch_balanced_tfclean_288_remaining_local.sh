#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
MODE="dry-run"
SELECTORS=()
LIMIT=""

while [[ "$#" -gt 0 ]]; do
  case "${1:-}" in
    --run)
      MODE="run"
      shift
      ;;
    --dry-run)
      shift
      ;;
    --limit)
      LIMIT="${2:-}"
      shift 2
      ;;
    *)
      break
      ;;
  esac
done

if [[ "$#" -gt 0 ]]; then
  SELECTORS=("$@")
else
  SELECTORS=("all_six")
fi

eval "$(/opt/anaconda3/bin/conda shell.zsh hook)"
conda activate background_project

cd "$ROOT_DIR"

python scripts/summarize_balanced_tfclean_288_progress.py "${SELECTORS[@]}" >/dev/null

REMAINING_SELECTORS=()
while IFS= read -r selector; do
  [[ -n "$selector" ]] && REMAINING_SELECTORS+=("$selector")
done < <(
  python - "$LIMIT" <<'PY'
import csv
import sys
from pathlib import Path

limit_arg = sys.argv[1]
limit = int(limit_arg) if limit_arg else None
summary_path = Path("analysis_outputs/tfclean_balanced_288_local_progress_current/progress_summary.csv")
rows = list(csv.DictReader(summary_path.open(encoding="utf-8")))
selected = 0
for row in rows:
    if row["launch_state"] == "not_started":
        print(row["canonical_selector"])
        selected += 1
        if limit is not None and selected >= limit:
            break
PY
)

if [[ "${#REMAINING_SELECTORS[@]}" -eq 0 ]]; then
  echo "No not-started local balanced_tfclean_288 configs remain."
  exit 0
fi

echo "Remaining selectors: ${REMAINING_SELECTORS[*]}"

if [[ "$MODE" == "dry-run" ]]; then
  bash run_balanced_tfclean_288_local.sh --dry-run "${REMAINING_SELECTORS[@]}"
else
  bash run_balanced_tfclean_288_local.sh --run "${REMAINING_SELECTORS[@]}"
fi
