#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"

resolve_config_name() {
  case "$1" in
    all_six)
      printf '%s\n' \
        "balanced_tfclean_flatearth_peaks_broad_integrated_288" \
        "balanced_tfclean_flatearth_plateaus_broad_integrated_288" \
        "balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288" \
        "balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288" \
        "balanced_tfclean_realistic_plateaus_wavy_broad_integrated_288" \
        "balanced_tfclean_realistic_plateaus_hilly_broad_integrated_288"
      ;;
    flatearth_pair)
      printf '%s\n' \
        "balanced_tfclean_flatearth_peaks_broad_integrated_288" \
        "balanced_tfclean_flatearth_plateaus_broad_integrated_288"
      ;;
    wavy_pair)
      printf '%s\n' \
        "balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288" \
        "balanced_tfclean_realistic_plateaus_wavy_broad_integrated_288"
      ;;
    hilly_pair)
      printf '%s\n' \
        "balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288" \
        "balanced_tfclean_realistic_plateaus_hilly_broad_integrated_288"
      ;;
    flatearth_peak_narrow)
      echo "balanced_tfclean_flatearth_peaks_broad_integrated_288"
      ;;
    flatearth_plateau_broad)
      echo "balanced_tfclean_flatearth_plateaus_broad_integrated_288"
      ;;
    wavy_peak_narrow)
      echo "balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288"
      ;;
    hilly_peak_narrow)
      echo "balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288"
      ;;
    wavy_plateau_broad)
      echo "balanced_tfclean_realistic_plateaus_wavy_broad_integrated_288"
      ;;
    hilly_plateau_broad)
      echo "balanced_tfclean_realistic_plateaus_hilly_broad_integrated_288"
      ;;
    *)
      echo "$1"
      ;;
  esac
}

CONFIGS=(
  "balanced_tfclean_flatearth_peaks_broad_integrated_288"
  "balanced_tfclean_flatearth_plateaus_broad_integrated_288"
  "balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288"
  "balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288"
  "balanced_tfclean_realistic_plateaus_wavy_broad_integrated_288"
  "balanced_tfclean_realistic_plateaus_hilly_broad_integrated_288"
)

if [[ "$#" -gt 0 ]]; then
  CONFIGS=()
  for selector in "$@"; do
    while IFS= read -r resolved; do
      [[ -n "$resolved" ]] && CONFIGS+=("$resolved")
    done < <(resolve_config_name "$selector")
  done
fi

eval "$(/opt/anaconda3/bin/conda shell.zsh hook)"
conda activate background_project

cd "$ROOT_DIR"

completed_stats_dirs=()

score_config() {
  local config_name="$1"
  local results_dir="results_${config_name}"
  local params_csv="${results_dir}/params/${config_name}_run_params.csv"
  local output_dir="analysis_outputs/${config_name}"
  local peak_count

  if [[ ! -d "$results_dir" ]]; then
    echo "Not score-ready: ${results_dir} does not exist yet" >&2
    return 1
  fi

  peak_count="$(find "$results_dir" -path '*/peaks/macs2/*_peaks.bed' -type f | wc -l | tr -d ' ')"
  if [[ "$peak_count" != "288" ]]; then
    echo "Not score-ready: ${results_dir} has ${peak_count}/288 MACS2 peak BEDs" >&2
    return 1
  fi

  python scripts/peak_pr_stats.py \
    --results-dir "$results_dir" \
    --params-csv "$params_csv" \
    --output-dir "$output_dir"

  completed_stats_dirs+=("$output_dir")
}

for config_name in "${CONFIGS[@]}"; do
  score_config "$config_name"
done

python scripts/build_balanced_288_config_report.py \
  --output-dir analysis_outputs/tfclean_balanced_288_local_current \
  --category-map configs/category_name_map.yaml \
  --input-dirs "${completed_stats_dirs[@]}"
