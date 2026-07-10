#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

reconstruct_one() {
  local target_rds="$1"
  local parts_dir="${target_rds}.parts"

  if [ ! -d "$parts_dir" ]; then
    echo "Missing parts directory: $parts_dir"
    return 1
  fi

  if ls "${parts_dir}"/*.part_* >/dev/null 2>&1; then
    cat "${parts_dir}"/*.part_* > "$target_rds"
    echo "Reconstructed: $target_rds"
  else
    echo "No part files found in: $parts_dir"
    return 1
  fi
}

reconstruct_one "results/objects/15_Figure4_like_scRNA_core_gene_expression_Seurat.rds"
reconstruct_one "results/objects/15b_Figure4_like_scRNA_cluster_core_gene_expression_Seurat.rds"
reconstruct_one "results/objects/16_Figure5_like_malignant_trajectory_core_gene_Seurat.rds"

echo "All split RDS objects reconstructed."
