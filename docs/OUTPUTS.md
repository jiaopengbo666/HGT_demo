# Output fields

`05_summary/HGT_summary.tsv`: one fragment pair per row.

| Group | Fields |
|---|---|
| Fragment | fragment_id, mag_a/b, contig_id_MAGs_a/b |
| DNA alignment | fragment_start/end_a/b, fragment_strand_a/b, mapping_identity_pct, mapping_alignment_length_bp, mapping_evalue, mapping_bitscore |
| HGT evidence | ani_pct, species_a/b, hgt_evidence |
| ORFs | orf_count_a/b, orf_id/start/end/strand_a/b |
| ORF coverage | orf_overlap_bp_a/b, orf_overlap_fraction_a/b, orf_relationship_a/b |
| ARG annotation | arg_status/subject_id/gene/subtype/type_a/b |
| ARG sets | arg_a, arg_b, shared_args, shared_complete_args, same_arg_in_both_fragments |
| Flags | arg_annotation_status, has_arg_overlap, contains_complete_arg_a/b, contains_complete_arg |

Within each endpoint, semicolon-separated ORF and annotation fields follow the same order (start, end, ID). NA keeps missing values in position. A and B are annotated independently.

Coordinates are 1-based inclusive. A is forward-oriented; B strand is relative to A. ORF coverage is the coordinate intersection divided by full ORF length. `contained` requires the entire ORF to lie inside the fragment; otherwise an intersection is `partial_overlap`.

`arg_a/b` list the endpoint ARGs. `shared_args` compares subtypes, using subject IDs when subtype is unavailable; `undetermined` marks an unresolved comparison. `shared_complete_args` requires full ORF containment on both endpoints. `contains_complete_arg` is Yes when either endpoint fully contains an ARG ORF. `non-ARG` means no passing hit; `not_assessed` means annotation was skipped.

| Path | Contents |
|---|---|
| 00_inputs/ | Input paths, taxonomy and mapping tables |
| 01_blastn/ | Raw nucleotide alignments, selected and redundant fragments |
| 02_fastani/ | ANI inputs and results |
| 03_prodigal/ | ORF coordinates and sequences |
| 04_arg_annotation/ | Raw, filtered and best hits; ARG overlap details |
| 05_summary/hgt_fragments.fna | Both DNA endpoints |
| 05_summary/hgt_arg_orfs.fna | Full ARG ORF sequences with any fragment overlap |
| pipeline.log / run_summary.json | Progress, parameters and counts |

ARG protein identity, E-value and coverage are available in `04_arg_annotation/arg_best_hits.tsv`.
