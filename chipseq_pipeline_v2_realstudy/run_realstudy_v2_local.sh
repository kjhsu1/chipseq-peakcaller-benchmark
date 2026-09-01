#!/usr/bin/env bash
set -euo pipefail

MODE="${MODE:-dry-run}"
CORES="${CORES:-4}"
SNAKEMAKE_HOME="${SNAKEMAKE_HOME:-/private/tmp/chipseq_snakemake_home}"
PIPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${PIPE_DIR}"

if [[ -f ../snakemake_stuff/setup.sh ]]; then
  source ../snakemake_stuff/setup.sh
fi

if [[ "${CONDA_PREFIX:-}" == *"/quobyte/"* || -z "${CONDA_PREFIX:-}" ]]; then
  if [[ -x /opt/anaconda3/bin/conda ]]; then
    eval "$(/opt/anaconda3/bin/conda shell.bash hook)"
    conda activate background_project
  fi
fi

export HOME="${SNAKEMAKE_HOME}"
export PYTHONDONTWRITEBYTECODE=1

command=(snakemake -s Snakefile.py --configfile config.yaml configs/realstudy_v2.yaml --cores "${CORES}" --rerun-incomplete)
if [[ "${MODE}" == "dry-run" ]]; then
  command+=(--dry-run --printshellcmds)
elif [[ "${MODE}" != "run" ]]; then
  echo "MODE must be dry-run or run" >&2
  exit 2
fi

"${command[@]}" "$@"
