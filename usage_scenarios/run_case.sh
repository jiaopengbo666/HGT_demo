#!/usr/bin/env bash
set -euo pipefail
task_project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
config="$task_project_dir/usage_scenarios/server_config.sh"
dry_run=false
scenario=""
usage() {
  cat <<'HELP'
Usage: bash usage_scenarios/run_case.sh CASE [--config FILE] [--dry-run]
  auto          MAGs + taxonomy + ARG database; compute all remaining steps.
  orfs          Reuse matched Prodigal GFF + proteins; annotate ARGs automatically.
  arg_hits      Reuse Prodigal GFF + proteins + raw ARG hits; compute DNA/ANI.
  ani           Reuse FastANI; compute BLASTn, ORFs and ARG annotation.
  blastn        Reuse raw BLASTn; compute ANI, ORFs and ARG annotation.
  reuse_all     Reuse BLASTn + ANI + Prodigal GFF/proteins + ARG hits.
  reuse_tsv     As reuse_all, but with normalized ORF TSV instead of GFF.
  filtered_hits As reuse_all, but trust already-filtered standard ARG hits.
  hgt_only      MAGs + taxonomy; no ARG assessment or ORF prediction.

Edit server_config.sh or provide --config. --dry-run prints without executing.
Output: OUTPUT_BASE/CASE (must be new or empty).
HELP
}
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    --config) [[ $# -ge 2 ]] || { echo 'Missing --config FILE' >&2; exit 2; }; config="$2"; shift 2 ;;
    --dry-run) dry_run=true; shift ;;
    -*) echo "Unknown option: $1" >&2; exit 2 ;;
    *) [[ -z "$scenario" ]] || { echo 'Provide one CASE only' >&2; exit 2; }; scenario="$1"; shift ;;
  esac
done
[[ -n "$scenario" ]] || { usage; exit 2; }
[[ -f "$config" ]] || { echo "Config not found: $config" >&2; exit 2; }
source "$config"
options=(-i "$MAG_DIR" -o "$OUTPUT_BASE/$scenario" --taxonomy "$TAXONOMY_TSV" -t "$THREADS"
         --min-fragment-length "${MIN_FRAGMENT_LENGTH:-500}")
case "$scenario" in
  auto|ani|blastn|orfs)
    options+=(--arg-db "$SARG_FASTA" --annotator "${ANNOTATOR:-blastp}") ;;
  arg_hits|reuse_all|reuse_tsv) options+=(--arg-hits "$ARG_HITS") ;;
  filtered_hits) options+=(--arg-hits "$FILTERED_ARG_HITS" --arg-hits-filtered) ;;
  hgt_only) options+=(--hgt-only) ;;
  *) echo "Unknown CASE: $scenario" >&2; usage; exit 2 ;;
esac
case "$scenario" in
  orfs|arg_hits|reuse_all|filtered_hits) options+=(--orf-gff "$ORF_GFF" --orf-proteins "$ORF_PROTEINS") ;;
  reuse_tsv) options+=(--orf-tsv "$ORF_TSV" --orf-proteins "$ORF_PROTEINS") ;;
esac
case "$scenario" in ani|reuse_all|reuse_tsv|filtered_hits) options+=(--ani-results "$ANI_RESULTS") ;; esac
case "$scenario" in blastn|reuse_all|reuse_tsv|filtered_hits) options+=(--blast-results "$BLASTN_RESULTS") ;; esac
if [[ "$scenario" != hgt_only ]]; then
  options+=(--arg-identity "${ARG_IDENTITY:-80}" --arg-query-cover "${ARG_QUERY_COVER:-75}" --arg-evalue "${ARG_EVALUE:-1e-10}")
  [[ -z "${ARG_MAP:-}" ]] || options+=(--arg-map "$ARG_MAP")
fi
case "${RUNNER:-server}" in
  server) command=(bash "$task_project_dir/run_server_existing_envs.sh") ;;
  python) command=("${PYTHON_BIN:-python3}" "$task_project_dir/hgt_demo.py") ;;
  *) echo 'RUNNER must be server or python' >&2; exit 2 ;;
esac
if "$dry_run"; then
  printf '%q ' "${command[@]}" "${options[@]}"; printf '\n'
else
  exec "${command[@]}" "${options[@]}"
fi
