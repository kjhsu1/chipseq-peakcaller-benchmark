#!/usr/bin/env bash
# Script will complement the codex "Setup Script"
# In AGENTS.md, tell it to source this before anything

CONDA="/cvmfs/hpc.ucdavis.edu/sw/conda/root/bin/conda"
ENV="$HOME/.conda/envs/background_project"

eval "$("$CONDA" shell.bash hook)"
conda activate "$ENV"
