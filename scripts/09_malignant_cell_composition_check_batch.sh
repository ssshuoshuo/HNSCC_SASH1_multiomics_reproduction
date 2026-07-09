#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/.."

mkdir -p results/tables

Rscript scripts/09_malignant_cell_composition_check.R 2>&1 | tee results/tables/09_malignant_cell_composition_check_run.log
