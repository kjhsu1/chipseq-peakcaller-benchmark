#!/usr/bin/env bash

set -euo pipefail

source ../snakemake_stuff/setup.sh

snakemake -s Snakefile.py --configfile configs/hilly_density_lite_map10.yaml --cores 8
snakemake -s Snakefile.py --configfile configs/hilly_density_lite_map20.yaml --cores 8
