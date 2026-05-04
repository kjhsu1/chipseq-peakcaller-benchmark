#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/.."
sbatch --parsable slurm/balanced_tfclean_288_ikorfgrp_high_25cpu_singlejob_scratch.sbatch
