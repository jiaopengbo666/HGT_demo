#!/usr/bin/env bash
# Portable defaults use detached testdata. Edit paths for your own genomes/database.
MAG_DIR="$task_project_dir/testdata/mags"
TAXONOMY_TSV="$task_project_dir/testdata/taxonomy.tsv"
SARG_FASTA="$task_project_dir/SARG.faa" # Set the path to your real ARG protein database.
PREPARED_INPUTS="$task_project_dir/testdata"
ORF_GFF="$PREPARED_INPUTS/genes.gff"
ORF_TSV="$PREPARED_INPUTS/orfs.tsv"
ORF_PROTEINS="$PREPARED_INPUTS/genes.faa"
ARG_HITS="$PREPARED_INPUTS/arg_hits.tsv"
FILTERED_ARG_HITS="$PREPARED_INPUTS/arg_hits.prefiltered.tsv"
ANI_RESULTS="$PREPARED_INPUTS/fastani.tsv"
BLASTN_RESULTS="$PREPARED_INPUTS/blastn.tsv"
ARG_MAP="$PREPARED_INPUTS/arg_map.tsv"
OUTPUT_BASE="$task_project_dir/scenario_results"
THREADS=4
ANNOTATOR=blastp
ARG_IDENTITY=80
ARG_QUERY_COVER=75
ARG_EVALUE=1e-10
MIN_FRAGMENT_LENGTH=500
RUNNER=python # python: active environment; server: separately configured existing environments.
PYTHON_BIN="${PYTHON:-python3}"
