#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"

MODE="dry-run"
CORES="${CORES:-4}"

resolve_config() {
  case "$1" in
    all_six)
      printf '%s\n' \
        "configs/balanced_tfclean_flatearth_peaks_broad_integrated_288.yaml" \
        "configs/balanced_tfclean_flatearth_plateaus_broad_integrated_288.yaml" \
        "configs/balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288.yaml" \
        "configs/balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288.yaml" \
        "configs/balanced_tfclean_realistic_plateaus_wavy_broad_integrated_288.yaml" \
        "configs/balanced_tfclean_realistic_plateaus_hilly_broad_integrated_288.yaml"
      ;;
    flatearth_pair)
      printf '%s\n' \
        "configs/balanced_tfclean_flatearth_peaks_broad_integrated_288.yaml" \
        "configs/balanced_tfclean_flatearth_plateaus_broad_integrated_288.yaml"
      ;;
    wavy_pair)
      printf '%s\n' \
        "configs/balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288.yaml" \
        "configs/balanced_tfclean_realistic_plateaus_wavy_broad_integrated_288.yaml"
      ;;
    hilly_pair)
      printf '%s\n' \
        "configs/balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288.yaml" \
        "configs/balanced_tfclean_realistic_plateaus_hilly_broad_integrated_288.yaml"
      ;;
    flatearth_peak_narrow)
      echo "configs/balanced_tfclean_flatearth_peaks_broad_integrated_288.yaml"
      ;;
    flatearth_plateau_broad)
      echo "configs/balanced_tfclean_flatearth_plateaus_broad_integrated_288.yaml"
      ;;
    wavy_peak_narrow)
      echo "configs/balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288.yaml"
      ;;
    hilly_peak_narrow)
      echo "configs/balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288.yaml"
      ;;
    wavy_plateau_broad)
      echo "configs/balanced_tfclean_realistic_plateaus_wavy_broad_integrated_288.yaml"
      ;;
    hilly_plateau_broad)
      echo "configs/balanced_tfclean_realistic_plateaus_hilly_broad_integrated_288.yaml"
      ;;
    *)
      echo "$1"
      ;;
  esac
}

if [[ "${1:-}" == "--run" ]]; then
  MODE="run"
  shift
elif [[ "${1:-}" == "--dry-run" ]]; then
  shift
fi

CONFIGS=(
  "configs/balanced_tfclean_flatearth_peaks_broad_integrated_288.yaml"
  "configs/balanced_tfclean_flatearth_plateaus_broad_integrated_288.yaml"
  "configs/balanced_tfclean_realistic_peaks_wavy_narrow_integrated_288.yaml"
  "configs/balanced_tfclean_realistic_peaks_hilly_narrow_integrated_288.yaml"
  "configs/balanced_tfclean_realistic_plateaus_wavy_broad_integrated_288.yaml"
  "configs/balanced_tfclean_realistic_plateaus_hilly_broad_integrated_288.yaml"
)

if [[ "$#" -gt 0 ]]; then
  CONFIGS=()
  for selector in "$@"; do
    while IFS= read -r resolved; do
      [[ -n "$resolved" ]] && CONFIGS+=("$resolved")
    done < <(resolve_config "$selector")
  done
fi

eval "$(/opt/anaconda3/bin/conda shell.zsh hook)"
conda activate background_project

export HOME="/private/tmp/chipseq_snakemake_home"
mkdir -p "$HOME"

cd "$ROOT_DIR"

for config_path in "${CONFIGS[@]}"; do
  config_name="$(basename "$config_path" .yaml)"
  result_root="results_${config_name}"
  params_table="${result_root}/params/${config_name}_run_params.csv"

  echo "=== ${MODE}: ${config_path} -> ${result_root} ==="
  if [[ "$MODE" == "dry-run" ]]; then
    snakemake \
      -s Snakefile.py \
      --cores "$CORES" \
      --configfile "$config_path" \
      --config result_root="$result_root" params_table="$params_table" \
      --dry-run
  else
    snakemake \
      -s Snakefile.py \
      --cores "$CORES" \
      --configfile "$config_path" \
      --config result_root="$result_root" params_table="$params_table"
  fi
done
