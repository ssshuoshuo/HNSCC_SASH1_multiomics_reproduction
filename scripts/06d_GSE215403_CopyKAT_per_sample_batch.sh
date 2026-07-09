#!/bin/bash
set -euo pipefail

PROJECT_DIR="/exports/eddie/scratch/s2469905/HNSCC_SASH1_reproduction"

# Grid Engine batch shell 需要显式加载 Environment Modules。
# 手动初始化后，module load 才能使用。

if [ -f /etc/profile.d/modules.sh ]; then
  source /etc/profile.d/modules.sh
elif [ -f /usr/share/Modules/init/bash ]; then
  source /usr/share/Modules/init/bash
else
  echo "ERROR: Environment Modules initialization script not found."
  exit 1
fi

module load R/4.5
export R_LIBS_USER="/exports/eddie/scratch/s2469905/R_libs"
export TMPDIR="/exports/eddie/scratch/s2469905/tmp"
export BIOMART_CACHE="/exports/eddie/scratch/s2469905/biomart_cache"

mkdir -p "$R_LIBS_USER"
mkdir -p "$TMPDIR"
mkdir -p "$BIOMART_CACHE"
mkdir -p "$PROJECT_DIR/logs"

cd "$PROJECT_DIR"

echo "============================================================"
echo "06d CopyKAT batch job started"
date
hostname
echo "============================================================"

Rscript scripts/06d_CopyKAT_malignant_call.R \
  > logs/logs_06d_CopyKAT_batch.txt 2>&1

echo "============================================================"
echo "06d CopyKAT batch job finished"
date
echo "============================================================"
