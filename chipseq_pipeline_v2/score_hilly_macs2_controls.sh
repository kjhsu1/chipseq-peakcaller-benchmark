#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"

eval "$(/opt/anaconda3/bin/conda shell.zsh hook)"
conda activate background_project

cd "$ROOT_DIR"

score_category() {
  local results_dir="$1"
  local params_csv="$2"
  local output_dir="$3"

  local peak_count
  peak_count="$(find "$results_dir" -path '*/peaks/macs2/*_peaks.bed' -type f | wc -l | tr -d ' ')"

  if [[ "$peak_count" != "128" ]]; then
    echo "Not score-ready: $results_dir has $peak_count/128 MACS2 peak BEDs" >&2
    return 1
  fi

  python scripts/peak_pr_stats.py \
    --results-dir "$results_dir" \
    --params-csv "$params_csv" \
    --output-dir "$output_dir"
}

score_category \
  "results_macs2_control_tfclean_realistic_peaks_hilly_narrow_128" \
  "results_macs2_control_tfclean_realistic_peaks_hilly_narrow_128/params/macs2_control_tfclean_realistic_peaks_hilly_narrow_integrated_128_run_params.csv" \
  "analysis_outputs/macs2_control_tfclean_realistic_peaks_hilly_narrow_128_stats"

score_category \
  "results_macs2_control_tfclean_realistic_plateaus_hilly_broad_128" \
  "results_macs2_control_tfclean_realistic_plateaus_hilly_broad_128/params/macs2_control_tfclean_realistic_plateaus_hilly_broad_integrated_128_run_params.csv" \
  "analysis_outputs/macs2_control_tfclean_realistic_plateaus_hilly_broad_128_stats"
