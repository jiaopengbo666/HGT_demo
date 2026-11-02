# Input formats

Tables are tab-separated. ORF coordinates are 1-based inclusive, with strand `+` or `-`.

| Input | Format |
|---|---|
| MAGs | One .fa/.fna/.fasta per MAG; unique MAG and contig IDs; paths free of whitespace |
| Taxonomy | Header: mag_id, species; or GTDB-Tk user_genome, classification |
| ORF GFF | CDS records with ID attributes |
| ORF TSV | Header: orf_id, contig_id, start, end, strand |
| ORF proteins | Protein FASTA from the same prediction as the coordinates |
| BLASTn / ARG hits | Standard 12-column outfmt6, headerless |
| FastANI | Query, reference, ANI, matched fragments, total fragments; headerless |
| ARG map | Header: subject_id, gene, subtype, type |
| ARG database | Protein FASTA; bundled SARG3.2.1-S by default, --arg-db overrides |

Standard outfmt6 column order:

```text
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
```

ARG hits may include qlen as column 13. For 12-column hits, supply `--orf-proteins` to calculate query coverage. `--arg-hits-filtered` imports user-filtered hits and skips the ARG thresholds.

Supply all ORF coordinates from the same prediction batch as the ARG query proteins. Prodigal GFF/protein IDs are matched through the FASTA header's ID attribute. The normalized ORF TSV uses protein IDs directly.

Taxonomy example:

```text
mag_id	species
MAG_A	s__Species A
MAG_B	s__Species B
```

ORF TSV example:

```text
orf_id	contig_id	start	end	strand
orf_A	contig_A	120	419	+
```
