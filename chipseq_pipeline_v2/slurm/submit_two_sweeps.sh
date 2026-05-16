#!/bin/bash

set -euo pipefail

: """Paths"""

REPO_ROOT="/home/kjhsu/Code/Korflab/chipseq-peakcaller-benchmark/chipseq_pipeline_v2"
PEAKS_SCRIPT="${REPO_ROOT}/slurm/realistic_peaks_integrated_enrich128.sbatch"
PLATEAU_SCRIPT="${REPO_ROOT}/slurm/realistic_plateau_integrated_enrich128.sbatch"


: """Submit"""

jid1=$(sbatch --parsable "${PEAKS_SCRIPT}")
echo "Submitted realistic_peaks_integrated_enrich128 as job ${jid1}"

jid2=$(sbatch --dependency=afterok:${jid1} --parsable "${PLATEAU_SCRIPT}")
echo "Submitted realistic_plateau_integrated_enrich128 as job ${jid2} after ${jid1}"
