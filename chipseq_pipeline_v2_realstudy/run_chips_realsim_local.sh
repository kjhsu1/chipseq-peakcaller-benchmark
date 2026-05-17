#!/bin/bash

set -euo pipefail

PIPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${PIPE_DIR}"

CONDA_ENV="${CONDA_ENV:-background_project}"
CONDA_BIN="${CONDA_BIN:-/opt/anaconda3/bin/conda}"
MODE="${MODE:-dry-run}"
CORES="${CORES:-4}"
LATENCY_WAIT="${LATENCY_WAIT:-30}"
SNAKEMAKE_HOME="${SNAKEMAKE_HOME:-/private/tmp/chipseq_snakemake_home}"
CONFIG_FILES="${CONFIG_FILES:-config.yaml configs/chips_local_full.yaml}"
TARGETS="${TARGETS:-}"
CHECK_ASSETS="${CHECK_ASSETS:-true}"

if [[ -z "${CONDA_PREFIX:-}" || "$(basename "${CONDA_PREFIX}")" != "${CONDA_ENV}" ]]; then
  eval "$("${CONDA_BIN}" shell.bash hook)"
  conda activate "${CONDA_ENV}"
fi

export HOME="${SNAKEMAKE_HOME}"
export PYTHONDONTWRITEBYTECODE=1

config_args=()
for config_file in ${CONFIG_FILES}; do
  config_args+=(--configfile "${config_file}")
done

if [[ "${CHECK_ASSETS}" == "true" && " ${CONFIG_FILES} " == *" configs/chips_local_full.yaml "* ]]; then
  missing=()
  [[ -f data/local/references/ce11/genome.fa ]] || missing+=("data/local/references/ce11/genome.fa")
  [[ -f data/local/references/ce11/genome.fa.fai ]] || missing+=("data/local/references/ce11/genome.fa.fai")
  if ! compgen -G "data/local/indexes/ce11/bowtie2/ce11*.bt2" > /dev/null; then
    missing+=("data/local/indexes/ce11/bowtie2/ce11*.bt2")
  fi

  if [[ "${#missing[@]}" -gt 0 ]]; then
    printf 'Local ChIPs full-run assets are not staged yet.\n' >&2
    printf 'Missing required local asset(s):\n' >&2
    for path in "${missing[@]}"; do
      printf '  - %s\n' "${path}" >&2
    done
    printf '\nStage those assets, or use configs/chips_tiny_smoke.yaml for the tiny smoke path.\n' >&2
    printf 'To bypass this preflight intentionally, set CHECK_ASSETS=false.\n' >&2
    exit 3
  fi
fi

cmd=(
  snakemake
  -s Snakefile.py
  -j "${CORES}"
  "${config_args[@]}"
  --config enable_chips_targets=true
  --latency-wait "${LATENCY_WAIT}"
  --rerun-incomplete
)

if [[ "${MODE}" == "dry-run" ]]; then
  cmd+=(--dry-run)
elif [[ "${MODE}" != "run" ]]; then
  echo "MODE must be 'dry-run' or 'run', got '${MODE}'" >&2
  exit 2
fi

for target in ${TARGETS}; do
  cmd+=("${target}")
done

printf 'Running local ChIPs realstudy command:\n'
printf '  %q' "${cmd[@]}"
printf '\n'

"${cmd[@]}"
