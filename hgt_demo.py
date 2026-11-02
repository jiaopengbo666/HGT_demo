#!/usr/bin/env python3
"""Find candidate HGT fragments and annotate overlapping ARG ORFs."""
import argparse
import csv
import json
import math
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

VERSION = "1.5.0"
DEFAULT_ARG_MAP = Path(__file__).resolve().parent / "data" / "SARG_names.tsv"


def fasta(path):
    name, description, seq = None, None, []
    with open(path, encoding="utf-8-sig") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, description, "".join(seq).upper()
                description = line[1:]
                if not description.split():
                    raise ValueError(f"Empty FASTA identifier: {path}")
                name, seq = description.split()[0], []
            else:
                if name is None:
                    raise ValueError(f"Sequence before FASTA header: {path}")
                seq.append(line)
    if name is not None:
        yield name, description, "".join(seq).upper()


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTRYMKBDHVN", "TGCAYRKMVHDBN"))[::-1]


def extract(seq, start, end, strand="+"):
    if not 1 <= start <= end <= len(seq):
        raise ValueError(f"Invalid sequence interval {start}..{end} for length {len(seq)}")
    region = seq[start - 1:end]
    return revcomp(region) if strand == "-" else region


def write_fasta(handle, name, seq):
    handle.write(f">{name}\n")
    for pos in range(0, len(seq), 80):
        handle.write(seq[pos:pos + 80] + "\n")


def tsv(path, fields, rows):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def table(path):
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Empty table: {path}")
        return list(reader)


def finite(value, label):
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"Non-finite {label}")
    return number


def load_mags(directory):
    files = sorted(p.resolve() for p in Path(directory).iterdir()
                   if p.is_file() and p.suffix.lower() in (".fa", ".fna", ".fasta"))
    if len(files) < 2:
        raise ValueError("Provide at least two MAG nucleotide FASTA files (.fa/.fna/.fasta).")
    sequences, owners, paths = {}, {}, {}
    for path in files:
        if path.stem in paths:
            raise ValueError(f"MAG stems must be unique: {path.stem}")
        if any(c.isspace() for c in str(path)):
            raise ValueError("MAG paths must not contain whitespace (FastANI list format).")
        paths[path.stem] = path
        count = 0
        for name, _, seq in fasta(path):
            if name in sequences:
                raise ValueError(f"Duplicate contig ID {name}; IDs must be globally unique.")
            if not seq or set(seq) - set("ACGTRYMKSWBDHVNU"):
                raise ValueError(f"Invalid nucleotide sequence: {name}")
            sequences[name], owners[name] = seq, path.stem
            count += 1
        if not count:
            raise ValueError(f"Empty MAG FASTA: {path}")
    return sequences, owners, paths


def resolve_mag(value, paths):
    # Resolve a MAG stem, filename or path.
    name = str(value).replace("\\", "/").rsplit("/", 1)[-1]
    if name in paths:
        return name
    matches = [mag for mag, path in paths.items() if path.name == name]
    if len(matches) == 1:
        return matches[0]
    raise ValueError(f"Unknown MAG identifier: {value}")


def load_taxonomy(path, paths):
    result = {}
    for row in table(path):
        value = row.get("mag_id") or row.get("user_genome")
        if not value or not ("species" in row or "classification" in row):
            raise ValueError("Taxonomy requires mag_id/species or GTDB user_genome/classification columns.")
        mag = resolve_mag(value, paths)
        species = row.get("species", "").strip()
        if "classification" in row:
            species = next((s.strip()[3:] for s in row["classification"].split(";")
                            if s.strip().startswith("s__")), "")
        species = species.removeprefix("s__").strip()
        if species.lower() in ("", "na", "n/a", "unknown", "unclassified", "none", "nan"):
            species = ""
        if mag in result and result[mag] != species:
            raise ValueError(f"Conflicting taxonomy for {mag}")
        result[mag] = species
    missing = sorted(set(paths) - set(result))
    if missing:
        raise ValueError("GTDB-Tk taxonomy is missing input MAGs: " + ", ".join(missing))
    return result


def load_ani(path, paths):
    values = defaultdict(list)
    with open(path, encoding="utf-8-sig") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) < 3:
                raise ValueError("FastANI table requires query, reference and ANI columns.")
            a, b = resolve_mag(cols[0], paths), resolve_mag(cols[1], paths)
            value = finite(cols[2], "ANI")
            if not 0 <= value <= 100:
                raise ValueError("ANI must be in percent, between 0 and 100.")
            if a != b:
                values[tuple(sorted((a, b)))].append(value)
    # Use the higher directional ANI estimate.
    return {key: max(nums) for key, nums in values.items()}


def load_alignments(path, sequences, owners, min_length, identity):
    records = {}
    with open(path, encoding="utf-8-sig") as handle:
        for n, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip().split("\t")
            if len(cols) < 12:
                raise ValueError(f"BLAST row {n}: expected standard outfmt 6 (12 columns).")
            q, s = cols[:2]
            if q not in sequences or s not in sequences:
                raise ValueError(f"BLAST row {n}: unknown contig ID {q} or {s}")
            ident, length = finite(cols[2], "identity"), int(cols[3])
            evalue, score = finite(cols[10], "E-value"), finite(cols[11], "bitscore")
            if not 0 <= ident <= 100 or length < 1 or evalue < 0:
                raise ValueError(f"Invalid BLAST row {n}")
            q1, q2, s1, s2 = map(int, cols[6:10])
            qa, qb, sa, sb = min(q1, q2), max(q1, q2), min(s1, s2), max(s1, s2)
            if not (1 <= qa <= qb <= len(sequences[q]) and 1 <= sa <= sb <= len(sequences[s])):
                raise ValueError(f"BLAST row {n}: coordinates outside contig")
            if owners[q] == owners[s] or length < min_length or ident <= identity:
                continue
            orientation = "+" if (q2 - q1) * (s2 - s1) >= 0 else "-"
            endpoints = sorted(((q, qa, qb), (s, sa, sb)))
            key = (tuple(endpoints), orientation)
            a, b = endpoints
            row = dict(contig_a=a[0], start_a=a[1], end_a=a[2],
                       contig_b=b[0], start_b=b[1], end_b=b[2],
                       mag_a=owners[a[0]], mag_b=owners[b[0]], strand_b=orientation,
                       identity=ident, alignment_length=length, evalue=evalue, bitscore=score)
            # Group reciprocal hits by coordinates.
            if key not in records or (score, ident, -evalue) > (
                    records[key]["bitscore"], records[key]["identity"], -records[key]["evalue"]):
                records[key] = row
    return [records[key] for key in sorted(records)]


def load_orfs(path, sequences):
    result = {}
    if path.suffix.lower() in (".gff", ".gff3"):
        rows = []
        with open(path, encoding="utf-8-sig") as handle:
            for line in handle:
                if line.startswith("##FASTA"):
                    break
                if not line.strip() or line.startswith("#"):
                    continue
                cols = line.rstrip().split("\t")
                if len(cols) != 9:
                    raise ValueError("GFF must have nine tab-separated columns.")
                if cols[2] != "CDS":
                    continue
                attrs = dict(item.split("=", 1) for item in cols[8].split(";") if "=" in item)
                rows.append(dict(orf_id=unquote(attrs.get("ID", "")), contig_id=cols[0],
                                 start=cols[3], end=cols[4], strand=cols[6]))
    else:
        rows = table(path)
    for row in rows:
        if not all(key in row for key in ("orf_id", "contig_id", "start", "end", "strand")):
            raise ValueError("ORF TSV requires orf_id, contig_id, start, end, strand.")
        oid, cid = row["orf_id"], row["contig_id"]
        if not oid or oid in result or cid not in sequences:
            raise ValueError(f"Unknown contig, empty or duplicate ORF ID: {oid}")
        a, b = int(row["start"]), int(row["end"])
        if not 1 <= a <= b <= len(sequences[cid]) or row["strand"] not in ("+", "-"):
            raise ValueError(f"Invalid ORF coordinates/strand: {oid}")
        result[oid] = dict(orf_id=oid, contig_id=cid, start=a, end=b, strand=row["strand"])
    return result


def progress(message, log_path=None, console=True):
    line = f"[{datetime.now().isoformat(timespec='seconds')}] {message}"
    if console:
        print(line, flush=True)
    if log_path is not None:
        with log_path.open("a", encoding="utf-8") as log:
            log.write(line + "\n")


def run_command(command, work):
    if not shutil.which(command[0]):
        progress(f"[TOOL] FAILED: required executable not found: {command[0]}", work.parent / "pipeline.log")
        raise ValueError(f"Required executable not found: {command[0]}")
    logs = work / "logs"
    logs.mkdir(exist_ok=True)
    with open(logs / "commands.jsonl", "a", encoding="utf-8") as log:
        log.write(json.dumps(command) + "\n")
    started = time.monotonic()
    pipeline_log = work.parent / "pipeline.log"
    progress(f"[TOOL] START {command[0]} | detailed log: {logs / 'tools.log'}", pipeline_log, console=False)
    with open(logs / "tools.log", "a", encoding="utf-8") as log:
        log.write(f"\n[{datetime.now().isoformat(timespec='seconds')}] COMMAND: {json.dumps(command)}\n")
        log.flush()
        try:
            subprocess.run(command, stdout=log, stderr=log, check=True)
        except (OSError, subprocess.CalledProcessError) as error:
            progress(f"[TOOL] FAILED {command[0]}: {error}", pipeline_log)
            raise
    progress(f"[TOOL] DONE {command[0]} ({time.monotonic() - started:.1f} s)", pipeline_log, console=False)


def reconcile_orf_ids(orfs, proteins):
    """Match Prodigal GFF and protein IDs through header attributes."""
    if not proteins:
        return orfs
    result, used = dict(orfs), set()
    for oid, description, _ in fasta(proteins):
        if oid in used:
            raise ValueError(f"Duplicate protein ID: {oid}")
        used.add(oid)
        if oid in orfs:
            continue
        match = re.search(r"\bID=([^;\s]+)", description)
        alias = match.group(1) if match else None
        if alias not in orfs or alias not in result:
            raise ValueError(f"Cannot map protein ID {oid} to ORF coordinates.")
        result[oid] = dict(result.pop(alias), orf_id=oid)
    return result


def generate_alignments(paths, args, work):
    databases = work / "db"
    databases.mkdir()
    for index, (mag, path) in enumerate(paths.items()):
        progress(f"[01.1] Building BLASTn database {index + 1}/{len(paths)}: {mag}", work.parent / "pipeline.log", console=False)
        run_command(["makeblastdb", "-in", str(path), "-dbtype", "nucl", "-parse_seqids",
                     "-out", str(databases / str(index))], work)
    result = work / "blastn.tsv"
    tsv(work / "database_manifest.tsv", ["database_index", "mag_id", "fasta_path"],
        [dict(database_index=index, mag_id=mag, fasta_path=str(path))
         for index, (mag, path) in enumerate(paths.items())])
    progress("[01.1] Nucleotide database construction complete", work.parent / "pipeline.log")
    progress("[01.2] Starting pairwise BLASTn alignments", work.parent / "pipeline.log")
    with result.open("w", encoding="utf-8") as merged:
        for qi, (qmag, query) in enumerate(paths.items()):
            for si, (smag, _) in enumerate(paths.items()):
                if qmag == smag:
                    continue
                pair = work / f"pair_{qi}_{si}.tsv"
                progress(f"[01.2] BLASTn {qmag} vs {smag}", work.parent / "pipeline.log", console=False)
                run_command(["blastn", "-query", str(query), "-db", str(databases / str(si)),
                             "-outfmt", "6", "-perc_identity", str(args.fragment_identity),
                             "-num_threads", str(args.threads), "-out", str(pair)], work)
                with pair.open() as source:
                    shutil.copyfileobj(source, merged)
    return result


def generate_ani(paths, args, work):
    manifest = work / "MAG_paths.txt"
    manifest.write_text("".join(str(p) + "\n" for p in paths.values()), encoding="utf-8")
    result = work / "fastani.tsv"
    run_command(["fastANI", "--ql", str(manifest), "--rl", str(manifest),
                 "-o", str(result), "-t", str(args.threads)], work)
    return result


def predict_orfs(sequences, work):
    contigs = work / "contigs.fna"
    with contigs.open("w") as handle:
        for cid, seq in sequences.items():
            write_fasta(handle, cid, seq)
    gff, proteins = work / "genes.gff", work / "genes.faa"
    run_command(["prodigal", "-i", str(contigs), "-p", "meta", "-f", "gff",
                 "-o", str(gff), "-a", str(proteins), "-d", str(work / "genes.fna")], work)
    orfs = load_orfs(gff, sequences)
    # Prodigal GFF IDs (e.g. 1_1) differ from FASTA IDs (e.g. contig_1).
    normalized = {}
    for oid, description, _ in fasta(proteins):
        match = re.search(r"\bID=([^;\s]+)", description)
        if not match or match.group(1) not in orfs:
            raise ValueError(f"Cannot map Prodigal protein to GFF: {oid}")
        normalized[oid] = dict(orfs[match.group(1)], orf_id=oid)
    if len(normalized) != len(orfs):
        raise ValueError("Prodigal GFF/protein counts do not match.")
    return normalized, proteins


def annotate(proteins, args, work):
    database_dir = work / "database"
    database_dir.mkdir(exist_ok=True)
    db = database_dir / "arg_db"
    hits = work / "arg_hits.tsv"
    fields = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"
    if args.annotator == "diamond":
        run_command(["diamond", "makedb", "--in", str(args.arg_db), "--db", str(db)], work)
        run_command(["diamond", "blastp", "--query", str(proteins), "--db", str(db),
                     "--out", str(hits), "--threads", str(args.threads),
                     "--evalue", str(args.arg_evalue), "--max-target-seqs", "5",
                     "--outfmt", "6"] + fields.split(), work)
    else:
        run_command(["makeblastdb", "-in", str(args.arg_db), "-dbtype", "prot", "-out", str(db)], work)
        run_command(["blastp", "-query", str(proteins), "-db", str(db), "-out", str(hits),
                     "-num_threads", str(args.threads), "-evalue", str(args.arg_evalue),
                     "-max_target_seqs", "5", "-outfmt", "6 " + fields], work)
    return hits


def load_arg_hits(path, orfs, proteins, args, audit_dir=None):
    lengths = {}
    if proteins:
        for oid, _, seq in fasta(proteins):
            if oid in lengths:
                raise ValueError(f"Duplicate protein ID: {oid}")
            lengths[oid] = len(seq)
    mapping = {}
    if args.arg_map:
        for row in table(args.arg_map):
            if not row.get("subject_id") or row["subject_id"] in mapping:
                raise ValueError("ARG map requires unique subject_id values.")
            mapping[row["subject_id"]] = row
    best = {}
    passing = []
    with open(path, encoding="utf-8-sig") as handle:
        for n, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip().split("\t")
            if len(cols) not in (12, 13):
                raise ValueError("ARG hits require standard outfmt 6: 12 columns, optionally followed by qlen.")
            oid, subject = cols[:2]
            if oid not in orfs:
                raise ValueError(f"ARG query {oid} absent from ORF coordinates; supply matching coordinates.")
            ident, aln = finite(cols[2], "ARG identity"), int(cols[3])
            ev, score = finite(cols[10], "ARG E-value"), finite(cols[11], "ARG bitscore")
            if not 0 <= ident <= 100 or aln <= 0 or ev < 0:
                raise ValueError(f"Invalid ARG hit row {n}")
            qlen = int(cols[12]) if len(cols) == 13 else lengths.get(oid)
            if qlen is not None and qlen <= 0:
                raise ValueError(f"Invalid qlen for {oid}")
            # Query coverage uses the protein coordinate span.
            span = abs(int(cols[7]) - int(cols[6])) + 1
            if qlen is not None and (min(int(cols[6]), int(cols[7])) < 1 or
                                     max(int(cols[6]), int(cols[7])) > qlen):
                raise ValueError(f"ARG query coordinates exceed qlen: {oid}")
            cov = 100 * span / qlen if qlen else None
            if not args.arg_hits_filtered:
                if cov is None:
                    raise ValueError("ARG coverage needs qlen (column 13) or --orf-proteins; "
                                     "use --arg-hits-filtered only for already filtered hits.")
                if ident < args.arg_identity or ev > args.arg_evalue or cov < args.arg_query_cover:
                    continue
            meta = mapping.get(subject, {})
            hit = dict(subject_id=subject, arg_gene=meta.get("gene", subject),
                       arg_subtype=meta.get("subtype", ""), arg_type=meta.get("type", ""),
                       arg_identity=ident, arg_evalue=ev, arg_query_cover=cov if cov is not None else "NA",
                       bitscore=score)
            if audit_dir is not None:
                passing.append(dict(orf_id=oid, **hit))
            # Rank passing hits by E-value, bitscore, identity and subject ID.
            rank = (ev, -score, -ident, subject)
            if oid not in best or rank < best[oid][0]:
                best[oid] = rank, hit
    if audit_dir is not None:
        tsv(audit_dir / "arg_hits.filtered.tsv",
            ["orf_id", "subject_id", "arg_gene", "arg_subtype", "arg_type", "arg_identity",
             "arg_evalue", "arg_query_cover", "bitscore"], passing)
    return {oid: value[1] for oid, value in best.items()}


FRAGMENT_FIELDS = ["fragment_id", "mag_a", "mag_b", "contig_a", "start_a", "end_a",
                   "contig_b", "start_b", "end_b", "strand_b", "identity", "alignment_length",
                   "evalue", "bitscore", "ani", "species_a", "species_b", "evidence", "status"]
ARG_FIELDS = ["fragment_id", "side", "mag_id", "partner_mag", "contig_id", "fragment_start",
              "fragment_end", "orf_id", "orf_start", "orf_end", "orf_strand", "overlap_bp",
              "orf_overlap_fraction", "relationship", "subject_id", "arg_gene", "arg_subtype",
              "arg_type", "arg_identity", "arg_evalue", "arg_query_cover"]


def orf_overlap(fragment_start, fragment_end, orf_start, orf_end):
    """Return overlap length, ORF coverage and containment."""
    overlap = max(0, min(fragment_end, orf_end) - max(fragment_start, orf_start) + 1)
    contained = fragment_start <= orf_start and orf_end <= fragment_end
    fraction = 1.0 if contained else min(overlap / (orf_end - orf_start + 1), math.nextafter(1.0, 0.0))
    return overlap, fraction, "contained" if contained else ("partial_overlap" if overlap else "outside")


def remove_contained_fragments(rows):
    """Remove alignments contained on both endpoints."""
    retained, redundant = [], []
    groups = defaultdict(list)
    for row in sorted(rows, key=lambda r: (-sum(r[f"end_{s}"] - r[f"start_{s}"] + 1 for s in ("a", "b")),
                                          -r["alignment_length"], -r["bitscore"], r["fragment_id"])):
        key = (row["contig_a"], row["contig_b"], row["strand_b"])
        parent = next((p for p in groups[key] if
                       all(p[f"start_{s}"] <= row[f"start_{s}"] <= row[f"end_{s}"] <= p[f"end_{s}"]
                           for s in ("a", "b"))), None)
        if parent:
            row["status"] = "redundant_contained"
            redundant.append(dict(discarded_fragment_id=row["fragment_id"], retained_fragment_id=parent["fragment_id"]))
        else:
            retained.append(row)
            groups[key].append(row)
    return sorted(retained, key=lambda r: r["fragment_id"]), redundant


MAPPING_COLUMNS = ["fragment_id", "mag_a", "mag_b", "contig_id_MAGs_a", "contig_id_MAGs_b",
                   "fragment_start_a", "fragment_end_a", "fragment_start_b", "fragment_end_b",
                   "fragment_strand_a", "fragment_strand_b", "mapping_identity_pct",
                   "mapping_alignment_length_bp", "mapping_evalue", "mapping_bitscore", "ani_pct",
                   "species_a", "species_b", "hgt_evidence"]
ORF_COLUMNS = ["orf_id", "orf_start", "orf_end", "orf_strand", "orf_overlap_bp",
               "orf_overlap_fraction", "orf_relationship"]
ANNOTATION_COLUMNS = ["arg_status", "arg_subject_id", "arg_gene", "arg_subtype", "arg_type"]
CONCLUSION_COLUMNS = ["arg_annotation_status", "has_arg_overlap", "contains_complete_arg_a",
                      "contains_complete_arg_b", "contains_complete_arg"]
REGION_ARG_COLUMNS = ["arg_a", "arg_b", "shared_args", "shared_complete_args", "same_arg_in_both_fragments"]
MASTER_COLUMNS = (MAPPING_COLUMNS + REGION_ARG_COLUMNS + ["orf_count_a", "orf_count_b"] +
                  [f"{field}_{side}" for side in ("a", "b") for field in ORF_COLUMNS] +
                  [f"{field}_{side}" for side in ("a", "b") for field in ANNOTATION_COLUMNS] +
                  CONCLUSION_COLUMNS)


def endpoint_arg_summary(features, annotation_enabled):
    """Compare endpoint ARG sets by subtype or reference ID."""
    if not annotation_enabled:
        return {key: "not_assessed" for key in REGION_ARG_COLUMNS}
    endpoints = {side: {} for side in ("a", "b")}
    contained = {side: set() for side in ("a", "b")}
    unmapped = {side: False for side in ("a", "b")}
    for feature in features:
        if feature["arg_status"] != "ARG":
            continue
        side = feature["orf_side"]
        subtype, subject = feature["arg_subtype"], feature["arg_subject_id"]
        known_subtype = subtype not in ("", "NA")
        # Different reference sequences may represent the same SARG subtype.
        key = ("subtype", subtype) if known_subtype else ("subject_id", subject)
        unmapped[side] |= not known_subtype
        endpoints[side][key] = feature["arg_gene"]
        if feature["orf_relationship"] == "contained":
            contained[side].add(key)
    shared = endpoints["a"].keys() & endpoints["b"].keys()
    complete_shared = contained["a"] & contained["b"]
    def labels(side, keys=None):
        return "; ".join(sorted({label for key, label in endpoints[side].items()
                                 if keys is None or key in keys})) or "None"
    same = "Yes" if shared else ("undetermined" if endpoints["a"] and endpoints["b"] and any(unmapped.values()) else "No")
    return dict(arg_a=labels("a"), arg_b=labels("b"), shared_args=labels("a", shared),
                shared_complete_args=labels("a", complete_shared), same_arg_in_both_fragments=same)


def build_master_summary(selected, orfs, hits, annotation_enabled):
    """Summarize each fragment pair with separate A/B ORF lists."""
    by_contig = defaultdict(list)
    for orf in orfs.values():
        by_contig[orf["contig_id"]].append(orf)
    summaries = []
    for fragment in selected:
        common = dict(fragment_id=fragment["fragment_id"], mag_a=fragment["mag_a"], mag_b=fragment["mag_b"],
                      contig_id_MAGs_a=fragment["contig_a"], contig_id_MAGs_b=fragment["contig_b"],
                      fragment_strand_a="+", fragment_strand_b=fragment["strand_b"],
                      mapping_identity_pct=fragment["identity"], mapping_alignment_length_bp=fragment["alignment_length"],
                      mapping_evalue=fragment["evalue"], mapping_bitscore=fragment["bitscore"],
                      ani_pct=fragment["ani"], species_a=fragment["species_a"], species_b=fragment["species_b"],
                      hgt_evidence=fragment["evidence"], arg_annotation_status="assessed" if annotation_enabled else "not_assessed")
        any_arg = False
        complete = {}
        features = []
        for side in ("a", "b"):
            start, end = fragment[f"start_{side}"], fragment[f"end_{side}"]
            common[f"fragment_start_{side}"], common[f"fragment_end_{side}"] = start, end
            complete[side] = False
            count = 0
            for orf in sorted(by_contig[fragment[f"contig_{side}"]], key=lambda o: (o["start"], o["end"], o["orf_id"])):
                overlap, fraction, relation = orf_overlap(start, end, orf["start"], orf["end"])
                if not overlap:
                    continue
                count += 1
                hit = hits.get(orf["orf_id"])
                any_arg |= bool(hit)
                complete[side] |= bool(hit and relation == "contained")
                features.append(dict(orf_side=side, orf_mag_id=fragment[f"mag_{side}"],
                                     orf_id=orf["orf_id"], orf_start=orf["start"], orf_end=orf["end"],
                                     orf_strand=orf["strand"], orf_overlap_bp=overlap, orf_overlap_fraction=fraction,
                                     orf_relationship=relation,
                                     arg_status=("ARG" if hit else "non-ARG") if annotation_enabled else "not_assessed",
                                     arg_subject_id=hit["subject_id"] if hit else "NA",
                                     arg_gene=hit["arg_gene"] if hit else "NA",
                                     arg_subtype=(hit["arg_subtype"] or "NA") if hit else "NA",
                                     arg_type=(hit["arg_type"] or "NA") if hit else "NA"))
            common[f"orf_count_{side}"] = count if orfs or annotation_enabled else "NA"
            common[f"contains_complete_arg_{side}"] = ("Yes" if complete[side] else "No") if annotation_enabled else "not_assessed"
        common["has_arg_overlap"] = ("Yes" if any_arg else "No") if annotation_enabled else "not_assessed"
        common["contains_complete_arg"] = ("Yes" if any(complete.values()) else "No") if annotation_enabled else "not_assessed"
        common.update(endpoint_arg_summary(features, annotation_enabled))
        for side in ("a", "b"):
            side_features = [feature for feature in features if feature["orf_side"] == side]
            for field in ORF_COLUMNS + ANNOTATION_COLUMNS:
                common[f"{field}_{side}"] = "; ".join(str(feature[field]) for feature in side_features) or "NA"
        summaries.append(common)
    return summaries


def analyze(args):
    progress("[INPUT] Reading MAGs and taxonomy")
    sequences, owners, paths = load_mags(args.input)
    taxonomy = load_taxonomy(args.taxonomy, paths)
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output directory must be new or empty.")
    output.mkdir(parents=True, exist_ok=True)
    pipeline_log = output / "pipeline.log"
    progress(f"[INPUT] Loaded {len(paths)} MAGs / {len(sequences)} contigs; HGT_demo {VERSION}", pipeline_log)
    stage_dirs = {name: output / name for name in
                  ("00_inputs", "01_blastn", "02_fastani", "03_prodigal", "04_arg_annotation", "05_summary")}
    for folder in stage_dirs.values():
        folder.mkdir()
    inputs, nucleotide, ani_dir, orf_dir, arg_dir, results = stage_dirs.values()
    sources = []
    def retain(source, destination, role):
        shutil.copyfile(source, destination)
        sources.append(dict(role=role, original_path=str(source), retained_path=str(destination.relative_to(output))))
        return destination
    shutil.copyfile(args.taxonomy, inputs / "taxonomy.tsv")
    tsv(inputs / "MAG_files.tsv", ["mag_id", "fasta_path"],
        [dict(mag_id=mag, fasta_path=str(path)) for mag, path in paths.items()])
    if args.arg_map:
        retain(args.arg_map, inputs / "arg_map.tsv", "arg_map")
        provenance = args.arg_map.with_name("SARG_names.source.json")
        if args.arg_map_source == "bundled_SARG" and provenance.is_file():
            retain(provenance, inputs / "arg_map.source.json", "arg_map_provenance")
    progress("[01] Reusing provided BLASTn results" if args.blast_results else "[01] Starting nucleotide database construction and BLASTn", pipeline_log)
    blast = (retain(args.blast_results, nucleotide / "blastn.tsv", "provided_blastn") if args.blast_results
             else generate_alignments(paths, args, nucleotide))
    progress("[01] BLASTn stage complete", pipeline_log)
    progress("[02] Reusing provided FastANI results" if args.ani_results else "[02] Starting FastANI", pipeline_log)
    ani_file = (retain(args.ani_results, ani_dir / "fastani.tsv", "provided_ani") if args.ani_results
                else generate_ani(paths, args, ani_dir))
    ani = load_ani(ani_file, paths)
    progress(f"[02] ANI stage complete: {len(ani)} MAG pairs with estimates", pipeline_log)
    candidates = load_alignments(blast, sequences, owners, args.min_fragment_length, args.fragment_identity)
    selected = []
    for n, row in enumerate(candidates, 1):
        row["fragment_id"] = f"HGT{n:06d}"
        pair = tuple(sorted((row["mag_a"], row["mag_b"])))
        value = ani.get(pair)
        sa, sb = taxonomy.get(row["mag_a"], ""), taxonomy.get(row["mag_b"], "")
        evidence = []
        if value is not None and value < args.max_ani:
            evidence.append("low_ani")
        if sa and sb and sa != sb:
            evidence.append("different_species")
        row.update(ani=value if value is not None else "NA", species_a=sa or "NA", species_b=sb or "NA",
                   evidence=";".join(evidence) or "none", status="selected" if evidence else
                   ("insufficient_evidence" if value is None or not (sa and sb) else "not_selected"))
        if evidence:
            selected.append(row)
    selected, redundant = remove_contained_fragments(selected)
    tsv(nucleotide / "redundant_fragment_pairs.tsv", ["discarded_fragment_id", "retained_fragment_id"], redundant)
    progress(f"[HGT] {len(candidates)} candidate alignments; {len(selected)} selected fragment pairs", pipeline_log)
    orfs, hits, proteins = {}, {}, args.orf_proteins
    coordinates = args.orf_coordinates
    if coordinates:
        coordinates = retain(coordinates, orf_dir / ("provided_coordinates" + coordinates.suffix), "provided_orf_coordinates")
    if proteins:
        proteins = retain(proteins, orf_dir / "genes.faa", "provided_orf_proteins")
    if args.arg_db:
        if coordinates:
            progress("[03] Loading provided ORF coordinates and proteins", pipeline_log)
            orfs = reconcile_orf_ids(load_orfs(coordinates, sequences), proteins)
        else:
            progress("[03] Starting Prodigal ORF prediction", pipeline_log)
            orfs, proteins = predict_orfs(sequences, orf_dir)
        progress(f"[03] ORF stage complete: {len(orfs)} ORFs", pipeline_log)
        progress(f"[04] Starting ARG database construction and {args.annotator} annotation", pipeline_log)
        hits = load_arg_hits(annotate(proteins, args, arg_dir), orfs, proteins, args, arg_dir)
    elif args.arg_hits:
        progress("[03] Loading ORFs paired with provided ARG hits", pipeline_log)
        orfs = reconcile_orf_ids(load_orfs(coordinates, sequences), proteins)
        raw_hits = retain(args.arg_hits, arg_dir / "arg_hits.tsv", "provided_arg_hits")
        progress("[04] Filtering provided ARG hits and selecting best hits", pipeline_log)
        hits = load_arg_hits(raw_hits, orfs, proteins, args, arg_dir)
    if orfs:
        tsv(orf_dir / "orf_coordinates.tsv", ["orf_id", "contig_id", "start", "end", "strand"],
            [orfs[oid] for oid in sorted(orfs)])
    if args.arg_db or args.arg_hits:
        tsv(arg_dir / "arg_best_hits.tsv",
            ["orf_id", "mag_id", "contig_id", "start", "end", "strand", "subject_id", "arg_gene",
             "arg_subtype", "arg_type", "arg_identity", "arg_evalue", "arg_query_cover", "bitscore"],
            [dict(**orfs[oid], mag_id=owners[orfs[oid]["contig_id"]], **hits[oid]) for oid in sorted(hits)])
        progress(f"[04] ARG stage complete: {len(hits)} ARG-positive ORFs across all input MAGs", pipeline_log)
    else:
        progress("[03-04] ORF prediction and ARG annotation skipped (HGT-only mode)", pipeline_log)
    tsv(inputs / "imported_files.tsv", ["role", "original_path", "retained_path"], sources)
    by_contig = defaultdict(list)
    for oid in hits:
        by_contig[orfs[oid]["contig_id"]].append(orfs[oid])
    arg_rows, arg_sequences = [], {}
    progress("[05] Associating ARG coordinates with HGT fragments and exporting results", pipeline_log)
    with (results / "hgt_fragments.fna").open("w", encoding="utf-8") as handle:
        for row in selected:
            for side, partner in (("a", "b"), ("b", "a")):
                cid, start, end = row[f"contig_{side}"], row[f"start_{side}"], row[f"end_{side}"]
                strand = row["strand_b"] if side == "b" else "+"
                fid = row["fragment_id"]
                write_fasta(handle, f"{fid}_{side}|mag={row[f'mag_{side}']}|contig={cid}|"
                            f"start={start}|end={end}|strand={strand}", extract(sequences[cid], start, end, strand))
                for orf in sorted(by_contig[cid], key=lambda x: (x["start"], x["orf_id"])):
                    overlap = max(0, min(end, orf["end"]) - max(start, orf["start"]) + 1)
                    if not overlap:
                        continue
                    overlap, fraction, relation = orf_overlap(start, end, orf["start"], orf["end"])
                    arg_rows.append(dict(fragment_id=fid, side=side, mag_id=row[f"mag_{side}"],
                                         partner_mag=row[f"mag_{partner}"], contig_id=cid,
                                         fragment_start=start, fragment_end=end, orf_id=orf["orf_id"],
                                         orf_start=orf["start"], orf_end=orf["end"], orf_strand=orf["strand"],
                                         overlap_bp=overlap, orf_overlap_fraction=fraction,
                                         relationship=relation,
                                         **hits[orf["orf_id"]]))
                    arg_sequences[orf["orf_id"]] = orf
    with (results / "hgt_arg_orfs.fna").open("w", encoding="utf-8") as handle:
        for oid, orf in sorted(arg_sequences.items()):
            write_fasta(handle, oid, extract(sequences[orf["contig_id"]], orf["start"], orf["end"], orf["strand"]))
    tsv(nucleotide / "candidate_alignments.tsv", FRAGMENT_FIELDS, candidates)
    tsv(nucleotide / "selected_fragment_pairs.tsv", FRAGMENT_FIELDS, selected)
    tsv(arg_dir / "fragment_arg_details.tsv", ARG_FIELDS, arg_rows)
    tsv(inputs / "contig_mag_map.tsv", ["contig_id", "mag_id"],
        [dict(contig_id=cid, mag_id=owners[cid]) for cid in sorted(owners)])
    master = build_master_summary(selected, orfs, hits, bool(args.arg_db or args.arg_hits))
    tsv(results / "HGT_summary.tsv", MASTER_COLUMNS, master)
    arg_fragment_ids = {r["fragment_id"] for r in master if r["has_arg_overlap"] == "Yes"}
    complete_arg_fragment_ids = {r["fragment_id"] for r in master if r["contains_complete_arg"] == "Yes"}
    report = dict(version=VERSION, parameters={k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                  input_mags=len(paths), candidate_fragment_pairs=len(candidates),
                  selected_fragment_pairs=len(selected), redundant_contained_pairs=len(redundant), arg_overlap_rows=len(arg_rows),
                  unique_arg_orfs=len(arg_sequences), annotation_enabled=bool(args.arg_db or args.arg_hits),
                  total_input_orfs=len(orfs), total_arg_positive_orfs=len(hits),
                  arg_map_source=args.arg_map_source,
                  arg_positive_fragment_pairs=len(arg_fragment_ids),
                  complete_arg_positive_fragment_pairs=len(complete_arg_fragment_ids), main_summary_rows=len(master),
                  ani_source="provided" if args.ani_results else "computed_with_fastani",
                  taxonomy_required=True,
                  output_layout={name: name + "/" for name in stage_dirs},
                  ani_policy="maximum available directional ANI; unavailable estimates are NA",
                  coordinate_system="1-based inclusive; both fragment endpoints exported; side b oriented to side a",
                  fragment_definition="individual BLAST HSPs; exact reciprocals and HSPs contained on both endpoints removed; partial overlaps not merged",
                  interpretation="candidate HGT regions; A/B denote alignment endpoints")
    (output / "run_summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    (results / "README.md").write_text(f"""# HGT summary

{len(paths)} MAGs; {len(selected)} fragment pairs; {len(master)} summary rows; {len(complete_arg_fragment_ids)} pairs with a complete ARG on at least one endpoint.

HGT_summary.tsv has one fragment pair per row. A/B ORF and annotation lists are semicolon-separated, ordered by ORF start, end and ID within each endpoint.

## Column groups

- DNA: MAGs, contigs, fragment coordinates/strand, BLASTn metrics, ANI and species.
- ORFs: predicted coordinates, strand, overlap length/fraction and containment for A/B.
- ARGs: per-ORF status, reference ID, gene, subtype and type; arg_a/b summarize each endpoint.
- Comparison: shared_args lists common annotations; shared_complete_args requires full containment on both endpoints.
- Flags: contains_complete_arg is Yes if either endpoint fully contains an ARG ORF; has_arg_overlap includes partial overlap.

Coordinates are 1-based inclusive. ORF coverage is the coordinate intersection divided by full ORF length; contained requires the ORF to lie entirely within the fragment. NA marks missing values; not_assessed marks skipped annotation.

hgt_fragments.fna contains both DNA endpoints. hgt_arg_orfs.fna contains full ARG ORF sequences with any fragment overlap. Protein alignment metrics are in ../04_arg_annotation/arg_best_hits.tsv; parameters and counts are in ../run_summary.json.
""", encoding="utf-8")
    progress(f"[05] DONE: {len(selected)} fragment pairs; {len(complete_arg_fragment_ids)} pairs containing complete ARGs; {len(arg_sequences)} distinct overlapping ARG ORFs", pipeline_log)
    progress(f"[OUTPUT] Final results: {results}", pipeline_log)
    tsv(output / "output_files.tsv", ["path", "size_bytes"],
        [dict(path=str(path.relative_to(output)), size_bytes=path.stat().st_size)
         for path in sorted(output.rglob("*")) if path.is_file()])


def parser():
    class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass

    p = argparse.ArgumentParser(
        description="HGT_demo v1.5: identify candidate shared DNA fragments between MAGs and locate ARG ORFs.",
        formatter_class=HelpFormatter,
        epilog="""Examples:
  python hgt_demo.py -i MAGs/ -o results/ --taxonomy taxonomy.tsv --arg-db SARG.faa
  python hgt_demo.py -i MAGs/ -o results/ --taxonomy taxonomy.tsv \
    --blast-results blastn.tsv --ani-results fastani.tsv \
    --arg-hits arg_hits.tsv --orf-gff genes.gff --orf-proteins genes.faa

Inputs:
  Taxonomy: GTDB-Tk user_genome/classification or mag_id/species TSV.
  BLASTn/ARG hits: standard 12-column outfmt6; ARG hits may add qlen in column 13.
  ORFs: GFF CDS records or TSV (orf_id, contig_id, start, end, strand).
  Imported ARG hits require matching ORF coordinates and query lengths
  from qlen or --orf-proteins. --arg-hits-filtered skips ARG thresholds.
  --arg-map accepts subject_id/gene/subtype/type TSV; default: data/SARG_names.tsv.
  MAG paths must be free of whitespace. Coordinates are 1-based inclusive.

Selection defaults:
  DNA: length >=500 bp, identity >99%, ANI <95% or different known species.
  ARGs: identity >=80%, query coverage >=75%, E-value <=1e-10; best passing hit.

Output:
  05_summary/HGT_summary.tsv: one fragment pair per row, separate A/B ORF lists.
  05_summary/: DNA fragment and ARG ORF FASTA exports.
  00_inputs/ to 04_arg_annotation/: inputs and intermediate results.
  pipeline.log and run_summary.json: progress, parameters and counts.
  See docs/INPUT_FORMATS.md and docs/OUTPUTS.md for field definitions.
""")
    p.add_argument("--version", action="version", version=VERSION)
    p.add_argument("-i", "--input", type=Path, required=True, help="MAG nucleotide FASTA directory (at least two MAGs)")
    p.add_argument("-o", "--output", type=Path, required=True, help="New/empty output directory")
    p.add_argument("-t", "--threads", type=int, default=13, help="Threads passed to external alignment and ANI tools")
    p.add_argument("--taxonomy", type=Path, required=True,
                   help="GTDB-Tk summary or species TSV covering all input MAGs")
    p.add_argument("--blast-results", type=Path, help="Reuse nucleotide BLASTn outfmt 6; otherwise compute MAG-to-MAG alignments")
    p.add_argument("--ani-results", type=Path,
                   help="Reuse FastANI results")
    p.add_argument("--min-fragment-length", type=int, default=500, help="Minimum nucleotide alignment length in bp, inclusive")
    p.add_argument("--fragment-identity", type=float, default=99, help="Strict lower bound: identity > this value")
    p.add_argument("--max-ani", type=float, default=95, help="Strict upper bound: ANI < this value")
    mode = p.add_mutually_exclusive_group()
    mode.add_argument("--arg-hits", type=Path, help="Existing ARG protein hits in outfmt 6; requires matching ORF coordinates")
    mode.add_argument("--arg-db", type=Path, help="ARG protein FASTA for automatic annotation")
    p.add_argument("--orf-gff", "--orf-tsv", dest="orf_coordinates", type=Path,
                   help="ORF GFF CDS coordinates or documented TSV; match ARG query IDs exactly")
    p.add_argument("--orf-proteins", type=Path, help="Matching ORF protein FASTA; required when reusing ORFs for automatic annotation")
    p.add_argument("--arg-map", type=Path,
                   help="Custom subject_id/gene/subtype/type TSV; ARG modes otherwise use bundled data/SARG_names.tsv")
    p.add_argument("--arg-hits-filtered", action="store_true", help="Import prefiltered ARG hits (skip ARG thresholds)")
    p.add_argument("--annotator", choices=("blastp", "diamond"), default="blastp",
                   help="ARG protein aligner")
    p.add_argument("--arg-identity", type=float, default=80,
                   help="Minimum ARG protein identity (%%), inclusive")
    p.add_argument("--arg-query-cover", type=float, default=75, help="Minimum query-protein coordinate coverage (%%)")
    p.add_argument("--arg-evalue", type=float, default=1e-10,
                   help="Maximum ARG E-value, inclusive")
    return p


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    if args.threads < 1 or args.min_fragment_length < 1:
        p.error("Threads and minimum fragment length must be positive.")
    for value in (args.fragment_identity, args.max_ani, args.arg_identity, args.arg_query_cover):
        if not math.isfinite(value) or not 0 <= value <= 100:
            p.error("Percentage thresholds must be finite and between 0 and 100.")
    if not math.isfinite(args.arg_evalue) or args.arg_evalue <= 0:
        p.error("ARG E-value must be finite and positive.")
    if args.arg_hits and not args.orf_coordinates:
        p.error("--arg-hits requires matching --orf-gff or --orf-tsv.")
    if args.arg_db and args.orf_coordinates and not args.orf_proteins:
        p.error("Existing ORF coordinates require --orf-proteins for automatic annotation.")
    if args.arg_db and args.orf_proteins and not args.orf_coordinates:
        p.error("Provide matching ORF coordinates with --orf-proteins, or let Prodigal predict both.")
    if args.arg_hits_filtered and not args.arg_hits:
        p.error("--arg-hits-filtered applies only to imported --arg-hits.")
    if not (args.arg_db or args.arg_hits) and (args.orf_coordinates or args.orf_proteins or args.arg_map):
        p.error("ORF/ARG annotation options require --arg-hits or --arg-db.")
    args.arg_map_source = "user" if args.arg_map else "not_used"
    try:
        if (args.arg_db or args.arg_hits) and args.arg_map is None:
            args.arg_map = DEFAULT_ARG_MAP
            args.arg_map_source = "bundled_SARG"
            if not args.arg_map.is_file():
                raise ValueError("Bundled data/SARG_names.tsv is missing. Upload the complete project including data/, "
                                 "or supply a custom --arg-map TSV.")
        if not args.input.is_dir():
            raise ValueError("Input MAG directory does not exist.")
        for key in ("taxonomy", "blast_results", "ani_results", "arg_hits", "arg_db",
                    "orf_coordinates", "orf_proteins", "arg_map"):
            value = getattr(args, key)
            if value is not None:
                if not value.is_file():
                    raise ValueError(f"Input file does not exist: {value}")
                setattr(args, key, value.resolve())
        analyze(args)
    except (ValueError, OSError, KeyError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
