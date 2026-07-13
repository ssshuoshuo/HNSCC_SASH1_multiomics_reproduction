#!/usr/bin/env bash
set -euo pipefail

cd /Users/yaoshuo/Desktop/HNSCC_SASH1_reproduction
mkdir -p results/logs

Rscript scripts/17_scTenifold_virtual_KO_OE.R \
  > results/logs/17_scTenifold_virtual_KO_OE.log \
  2>&1

echo "Script 17 completed."
echo "Log: results/logs/17_scTenifold_virtual_KO_OE.log"
