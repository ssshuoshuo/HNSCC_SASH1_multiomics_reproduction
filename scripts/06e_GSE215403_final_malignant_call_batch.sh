#!/bin/bash
set -euo pipefail

PROJECT_DIR="/exports/eddie/scratch/s2469905/HNSCC_SASH1_reproduction"

source /etc/profile.d/modules.sh
module load R/4.5

export R_LIBS_USER="/exports/eddie/scratch/s2469905/R_libs"
export TMPDIR="/exports/eddie/scratch/s2469905/tmp"
export BIOMART_CACHE="/exports/eddie/scratch/s2469905/biomart_cache"

mkdir -p "$R_LIBS_USER"
mkdir -p "$TMPDIR"
mkdir -p "$BIOMART_CACHE"
mkdir -p "$PROJECT_DIR/logs"

cd "$PROJECT_DIR"

Rscript scripts/06e_finalize_malignant_call.R \
  > logs/logs_06e_finalize_malignant_call.txt 2>&1
