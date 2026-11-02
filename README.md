# HGT_demo v1.5

Find candidate HGT fragments between MAGs and annotate the ARG ORFs they contain.

## Installation

Linux or WSL:

```bash
conda env create -f environment.yml
conda activate hgt-demo
```

Dependencies: Python >=3.10, BLAST+, FastANI, Prodigal and DIAMOND.

## Run

```bash
bash hgt_pipeline.sh \
  -i MAGs/ -o results/ \
  --taxonomy taxonomy.tsv \
  -t 8
```

Provide one nucleotide FASTA per MAG and a GTDB-Tk taxonomy table covering all input MAGs. MAG and contig IDs must be unique; MAG paths must be free of whitespace. Use a new output directory for each run.

The bundled `database/SARG3.2.1-S.fasta` is the default ARG database (12,746 proteins). BLASTp is the default aligner; select DIAMOND with `--annotator diamond`. ARG labels come from `data/SARG_names.tsv`.

Use `--arg-db /path/to/other.fasta` for another database and `--arg-map` for its labels. Use `--hgt-only` for DNA fragment analysis without ARG annotation.

Existing results can be supplied with `--blast-results`, `--ani-results`, `--orf-gff`/`--orf-tsv`, `--orf-proteins` and `--arg-hits`. See [input formats](docs/INPUT_FORMATS.md). Run `bash hgt_pipeline.sh -h` for all options.

## Database cache

The first ARG annotation run builds an index; later runs reuse it from `~/.cache/hgt_demo/arg_databases/` (or `$XDG_CACHE_HOME/hgt_demo/arg_databases/`). BLASTp and DIAMOND indexes are stored separately. A change to the protein FASTA content creates a new index.

Use `--arg-db-cache /path/to/cache` to choose the cache directory, or `--rebuild-arg-db` to rebuild. Each run records the index path and reuse status in `04_arg_annotation/database.json`.

## Default thresholds

- DNA: alignment length >=500 bp, identity >99%, and ANI <95% or different species.
- ARGs: protein identity >=80%, query coverage >=75%, E-value <=1e-10.

All thresholds are configurable. Reciprocal alignments and alignments contained on both endpoints are deduplicated.

## Output

`05_summary/HGT_summary.tsv` contains one fragment pair per row, with separate A/B ORF and ARG lists. The same folder contains both fragment sequences and overlapping ARG ORF sequences. Intermediate results are grouped by analysis step. See [output fields](docs/OUTPUTS.md).

## Test data

Place the separate `testdata/` folder next to `hgt_demo.py`, then run:

```bash
bash testdata/run_example.sh
```

The synthetic example uses precomputed results and requires only Python >=3.10. Expected output: **2 fragment pairs, 2 summary rows and PASS**.

## Reference

Method adapted from [Mahmud et al., Nature Microbiology (2024)](https://doi.org/10.1038/s41564-024-01639-4). Reference scripts: Bejan Mahmud and Luke Diorio-Toth. [Software references](REFERENCES.md).

The associated manuscript is in preparation and has not yet been published.
