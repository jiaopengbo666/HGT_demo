# HGT summary

2 MAGs; 2 fragment pairs; 2 summary rows; 1 pairs with a complete ARG on at least one endpoint.

HGT_summary.tsv has one fragment pair per row. A/B ORF and annotation lists are semicolon-separated, ordered by ORF start, end and ID within each endpoint.

## Column groups

- DNA: MAGs, contigs, fragment coordinates/strand, BLASTn metrics, ANI and species.
- ORFs: predicted coordinates, strand, overlap length/fraction and containment for A/B.
- ARGs: per-ORF status, reference ID, gene, subtype and type; arg_a/b summarize each endpoint.
- Comparison: shared_args lists common annotations; shared_complete_args requires full containment on both endpoints.
- Flags: contains_complete_arg is Yes if either endpoint fully contains an ARG ORF; has_arg_overlap includes partial overlap.

Coordinates are 1-based inclusive. ORF coverage is the coordinate intersection divided by full ORF length; contained requires the ORF to lie entirely within the fragment. NA marks missing values; not_assessed marks skipped annotation.

hgt_fragments.fna contains both DNA endpoints. hgt_arg_orfs.fna contains full ARG ORF sequences with any fragment overlap. Protein alignment metrics are in ../04_arg_annotation/arg_best_hits.tsv; parameters and counts are in ../run_summary.json.
