#!/bin/bash
set -euo pipefail

# Reconstruct split files after cloning this repository and running:
# git lfs install
# git lfs pull

while IFS= read -r original_file; do
  part_dir="${original_file}.parts"
  base_name="$(basename "$original_file")"

  if [ -d "$part_dir" ]; then
    echo "Reconstructing $original_file"
    cat "$part_dir/${base_name}".part_* > "$original_file"
  else
    echo "Missing part directory: $part_dir"
  fi
done < split_file_manifest/oversized_files.txt
