#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/.."

mkdir -p results/tables

Rscript scripts/07_CopyKAT_malignant_call.R 2>&1 | tee results/tables/07_CopyKAT_malignant_call_run.log
