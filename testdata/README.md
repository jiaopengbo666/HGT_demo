# testdata

Two synthetic MAGs with precomputed BLASTn, ANI, ORFs and ARG hits. Requires Python >=3.10.

```bash
bash testdata/run_example.sh
```

Expected: 2 fragment pairs, 2 summary rows, 4 fragment FASTA records and 3 overlapping ARG ORF sequences. The runner compares outputs with expected/ and prints PASS.

| Files | Contents |
|---|---|
| mags/ | Two 2,000-bp nucleotide sequences |
| taxonomy.tsv | Species assignments |
| blastn.tsv / fastani.tsv | Nucleotide hits / ANI |
| orfs.tsv | ORF coordinates |
| genes.gff / genes.faa | Matching GFF and proteins |
| arg_hits.tsv | ARG hits with qlen in column 13 |
| arg_hits.12col.tsv | Hits for use with --orf-proteins |
| arg_hits.prefiltered.tsv | Hits for --arg-hits-filtered |
| arg_map.tsv | Example ARG labels |
| expected/ | Expected table and sequences |

Use a new output name to repeat: `bash testdata/run_example.sh results_example_02`.
