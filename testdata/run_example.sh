#!/usr/bin/env bash
set -euo pipefail
testdata_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(cd -- "$testdata_dir/.." && pwd)"
output="${1:-$project_dir/results_example}"
bash "$project_dir/hgt_pipeline.sh" \
  -i "$testdata_dir/mags" -o "$output" \
  --taxonomy "$testdata_dir/taxonomy.tsv" \
  --blast-results "$testdata_dir/blastn.tsv" \
  --ani-results "$testdata_dir/fastani.tsv" \
  --orf-tsv "$testdata_dir/orfs.tsv" \
  --arg-hits "$testdata_dir/arg_hits.tsv" \
  --arg-map "$testdata_dir/arg_map.tsv"
"${PYTHON:-python3}" "$testdata_dir/check_output.py" "$output"
