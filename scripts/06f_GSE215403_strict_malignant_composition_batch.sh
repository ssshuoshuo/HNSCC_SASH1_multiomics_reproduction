#!/bin/bash
set -euo pipefail

PROJECT_DIR="/exports/eddie/scratch/s2469905/HNSCC_SASH1_reproduction"

source /etc/profile.d/modules.sh
module load R/4.5

export R_LIBS_USER="/exports/eddie/scratch/s2469905/R_libs"
export TMPDIR="/exports/eddie/scratch/s2469905/tmp"

mkdir -p "$R_LIBS_USER"
mkdir -p "$TMPDIR"
mkdir -p "$PROJECT_DIR/logs"

cd "$PROJECT_DIR"

Rscript scripts/06f_malignant_cell_composition_check.R \
  > logs/logs_06f_malignant_cell_composition_check.txt 2>&1
