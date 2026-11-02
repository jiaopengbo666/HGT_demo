"""Synthetic regression tests; external tool calls are mocked where indicated."""
import csv
import importlib.util
import io
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import hgt_demo as hgt


def rows(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def blast(q, s, qstart, qend, sstart, send, identity=99.5, length=600):
    return f"{q}\t{s}\t{identity}\t{length}\t0\t0\t{qstart}\t{qend}\t{sstart}\t{send}\t1e-100\t900\n"


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.mags = self.base / "mags"
        self.mags.mkdir()
        self.seq_a = "AACCGGTT" * 250
        self.seq_b = "ACGTTTGC" * 250
        (self.mags / "A.fa").write_text(">ctgA\n" + self.seq_a + "\n")
        (self.mags / "B.fa").write_text(">ctgB\n" + self.seq_b + "\n")
        self.blast = self.base / "blast.tsv"
        self.blast.write_text(blast("ctgA", "ctgB", 101, 700, 800, 201)
                              + blast("ctgB", "ctgA", 800, 201, 101, 700)
                              + blast("ctgA", "ctgB", 1001, 1600, 1001, 1600)
                              + blast("ctgA", "ctgB", 1, 600, 1, 600, identity=99)
                              + blast("ctgA", "ctgB", 1, 499, 1, 499, length=499))
        self.ani = self.base / "ani.tsv"
        self.ani.write_text("A.fa\tB.fa\t94\t10\t12\nB.fa\tA.fa\t94.2\t10\t12\n")
        self.taxonomy = self.base / "taxonomy.tsv"
        self.taxonomy.write_text("mag_id\tspecies\nA\ts__Same species\nB\ts__Same species\n")
        self.orfs = self.base / "orfs.tsv"
        self.orfs.write_text("orf_id\tcontig_id\tstart\tend\tstrand\n"
                             "argA\tctgA\t120\t419\t+\n"
                             "argPartial\tctgA\t650\t949\t+\n"
                             "argB\tctgB\t220\t519\t-\n"
                             "outside\tctgA\t1701\t2000\t+\n")
        self.hits = self.base / "arg.tsv"
        self.hits.write_text("".join(
            f"{oid}\tref1\t90\t100\t0\t0\t1\t100\t1\t100\t1e-30\t200\t100\n"
            for oid in ("argA", "argPartial", "argB", "outside")))
        self.map = self.base / "map.tsv"
        self.map.write_text("subject_id\tgene\tsubtype\ttype\nref1\tgeneX\tsubtypeX\ttypeX\n")

    def args(self, name="out", extra=()):
        return ["-i", str(self.mags), "-o", str(self.base / name),
                "--taxonomy", str(self.taxonomy), "--blast-results", str(self.blast),
                "--ani-results", str(self.ani)] + list(extra)

    def run_ok(self, name="out", extra=()):
        self.assertEqual(hgt.main(self.args(name, extra)), 0)
        return self.base / name

    def test_basic_boundaries_reciprocals_and_multiple_fragments(self):
        out = self.run_ok()
        selected = rows(out / "01_blastn/selected_fragment_pairs.tsv")
        self.assertEqual(len(selected), 2)
        self.assertEqual({r["start_a"] for r in selected}, {"101", "1001"})
        self.assertEqual(selected[0]["ani"], "94.2")
        self.assertEqual(selected[0]["strand_b"], "-")
        self.assertEqual(rows(out / "04_arg_annotation/fragment_arg_details.tsv"), [])

    def test_fragment_extraction_and_reverse_complement(self):
        out = self.run_ok()
        seqs = {name.split("|")[0]: seq for name, _, seq in hgt.fasta(out / "05_summary/hgt_fragments.fna")}
        self.assertEqual(seqs["HGT000001_a"], self.seq_a[100:700])
        self.assertEqual(seqs["HGT000001_b"], hgt.revcomp(self.seq_b[200:800]))
        self.assertEqual(len(seqs), 4)

    def test_arg_overlap_full_partial_outside_and_orf_sequence(self):
        out = self.run_ok(extra=["--arg-hits", str(self.hits), "--orf-tsv", str(self.orfs),
                                "--arg-map", str(self.map)])
        matches = {r["orf_id"]: r for r in rows(out / "04_arg_annotation/fragment_arg_details.tsv")}
        self.assertEqual(set(matches), {"argA", "argB", "argPartial"})
        self.assertEqual(matches["argPartial"]["relationship"], "partial_overlap")
        self.assertEqual(matches["argPartial"]["overlap_bp"], "51")
        self.assertEqual(matches["argB"]["relationship"], "contained")
        self.assertEqual(matches["argA"]["arg_gene"], "geneX")
        proteins = {name: seq for name, _, seq in hgt.fasta(out / "05_summary/hgt_arg_orfs.fna")}
        self.assertEqual(proteins["argB"], hgt.revcomp(self.seq_b[219:519]))
        summary = rows(out / "05_summary/HGT_summary.tsv")[0]
        self.assertEqual(summary["orf_count_a"], "2")
        self.assertEqual(summary["orf_count_b"], "1")
        self.assertEqual(summary["contains_complete_arg"], "Yes")

    def test_missing_ani_is_unknown_not_low(self):
        self.ani.write_text("")
        out = self.run_ok()
        self.assertEqual(rows(out / "01_blastn/selected_fragment_pairs.tsv"), [])
        self.assertEqual(rows(out / "01_blastn/candidate_alignments.tsv")[0]["status"], "insufficient_evidence")

    def test_master_lists_all_orfs_and_independent_endpoint_annotations(self):
        with self.orfs.open("a") as f:
            f.write("nonArgA\tctgA\t450\t599\t-\n")
        self.hits.write_text("argA\trefA\t90\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\t100\n"
                             "argB\trefB\t85\t80\t0\t0\t1\t80\t1\t80\t1e-12\t180\t100\n")
        out = self.run_ok(extra=["--arg-hits", str(self.hits), "--orf-tsv", str(self.orfs)])
        main = rows(out / "05_summary/HGT_summary.tsv")
        self.assertEqual(len(main), 2)
        summary = main[0]
        self.assertEqual(summary["orf_id_a"], "argA; nonArgA; argPartial")
        self.assertEqual(summary["orf_id_b"], "argB")
        self.assertEqual(summary["arg_status_a"], "ARG; non-ARG; non-ARG")
        self.assertEqual(summary["arg_status_b"], "ARG")
        self.assertEqual(summary["arg_subject_id_a"], "refA; NA; NA")
        self.assertEqual(summary["arg_subject_id_b"], "refB")
        self.assertEqual(main[-1]["orf_id_a"], "NA")
        self.assertEqual(main[-1]["orf_count_a"], "0")
        self.assertEqual((summary["fragment_start_a"], summary["fragment_end_a"]), ("101", "700"))
        self.assertEqual(summary["orf_start_a"], "120; 450; 650")
        self.assertEqual(summary["orf_relationship_a"], "contained; contained; partial_overlap")
        self.assertEqual(summary["orf_overlap_fraction_a"], "1.0; 1.0; 0.17")
        self.assertEqual((summary["arg_a"], summary["arg_b"]), ("refA", "refB"))
        self.assertEqual(summary["same_arg_in_both_fragments"], "undetermined")
        self.assertFalse(any(k.startswith("arg_db_") for k in summary))
        self.assertNotIn("orf_side", summary)
        for side in ("a", "b"):
            for field in hgt.ORF_COLUMNS + hgt.ANNOTATION_COLUMNS:
                self.assertEqual(len(summary[f"{field}_{side}"].split("; ")), int(summary[f"orf_count_{side}"]))
        report = json.loads((out / "run_summary.json").read_text())
        self.assertEqual(report["main_summary_rows"], 2)
        self.assertEqual(report["complete_arg_positive_fragment_pairs"], 1)
        best = {r['orf_id']: r for r in rows(out / "04_arg_annotation/arg_best_hits.tsv")}
        self.assertEqual(best['argB']['arg_identity'], '85.0')
        self.assertEqual(best['argB']['arg_query_cover'], '80.0')
        self.assertEqual(set(p.name for p in (out / "05_summary").glob('*.tsv')), {"HGT_summary.tsv"})

    def test_bundled_sarg_names_and_custom_map_override(self):
        self.hits.write_text("argA\tABA08387\t90\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\t100\n")
        opts = ["--arg-hits", str(self.hits), "--orf-tsv", str(self.orfs)]
        out = self.run_ok("sarg_names", opts)
        row = rows(out / "05_summary/HGT_summary.tsv")[0]
        self.assertEqual((row['arg_gene_a'], row['arg_subtype_a'], row['arg_type_a']),
                         ('dfrA1; NA', 'trimethoprim__dfrA1; NA', 'trimethoprim; NA'))
        self.assertEqual(json.loads((out / 'run_summary.json').read_text())['arg_map_source'], 'bundled_SARG')
        custom = self.base / 'custom_map.tsv'
        custom.write_text('subject_id\tgene\tsubtype\ttype\nABA08387\tcustom_gene\tcustom_subtype\tcustom_type\n')
        out = self.run_ok("custom_names", opts + ["--arg-map", str(custom)])
        self.assertEqual(rows(out / "05_summary/HGT_summary.tsv")[0]['arg_gene_a'], 'custom_gene; NA')
        self.assertEqual(json.loads((out / 'run_summary.json').read_text())['arg_map_source'], 'user')

    def test_missing_bundled_map_reports_error(self):
        with patch.object(hgt, 'DEFAULT_ARG_MAP', self.base / 'missing_map.tsv'):
            self.assertEqual(hgt.main(self.args(extra=['--arg-hits', str(self.hits), '--orf-tsv', str(self.orfs)])), 1)

    def test_endpoint_arg_sets_match_subtypes_across_reference_ids(self):
        features = [dict(orf_side='a', arg_status='ARG', arg_subject_id='alleleA',
                         arg_gene='geneX', arg_subtype='type__geneX', orf_relationship='contained'),
                    dict(orf_side='b', arg_status='ARG', arg_subject_id='alleleB',
                         arg_gene='geneX', arg_subtype='type__geneX', orf_relationship='contained')]
        summary = hgt.endpoint_arg_summary(features, True)
        self.assertEqual((summary['arg_a'], summary['arg_b'], summary['shared_args']), ('geneX', 'geneX', 'geneX'))
        self.assertEqual(summary['shared_complete_args'], 'geneX')
        self.assertEqual(summary['same_arg_in_both_fragments'], 'Yes')
        features[0]['orf_relationship'] = 'partial_overlap'
        summary = hgt.endpoint_arg_summary(features, True)
        self.assertEqual(summary['shared_args'], 'geneX')
        self.assertEqual(summary['shared_complete_args'], 'None')
        features[1]['arg_subtype'] = 'type__geneY'
        features[1]['arg_gene'] = 'geneY'
        summary = hgt.endpoint_arg_summary(features, True)
        self.assertEqual(summary['same_arg_in_both_fragments'], 'No')
        self.assertEqual(summary['shared_args'], 'None')

    def test_endpoint_arg_sets_missing_metadata_and_unassessed(self):
        self.assertEqual(hgt.endpoint_arg_summary([], False)['same_arg_in_both_fragments'], 'not_assessed')
        features = [dict(orf_side=s, arg_status='ARG', arg_subject_id='ref'+s,
                         arg_gene='ref'+s, arg_subtype='NA', orf_relationship='contained') for s in ('a', 'b')]
        self.assertEqual(hgt.endpoint_arg_summary(features, True)['same_arg_in_both_fragments'], 'undetermined')
        features[1]['arg_subject_id'] = 'refa'
        self.assertEqual(hgt.endpoint_arg_summary(features, True)['same_arg_in_both_fragments'], 'Yes')
        self.assertEqual(hgt.endpoint_arg_summary(features[:1], True)['same_arg_in_both_fragments'], 'No')

    def test_complete_arg_flag_excludes_partial_only_genes(self):
        self.hits.write_text("argPartial\tref\t90\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\t100\n")
        out = self.run_ok(extra=["--arg-hits", str(self.hits), "--orf-tsv", str(self.orfs)])
        summary = rows(out / "05_summary/HGT_summary.tsv")[0]
        self.assertEqual(summary["has_arg_overlap"], "Yes")
        self.assertEqual(summary["contains_complete_arg"], "No")
        self.assertEqual(summary["contains_complete_arg_a"], "No")

    def test_overlap_requires_positional_containment(self):
        self.assertEqual(hgt.orf_overlap(100, 199, 200, 299), (0, 0.0, "outside"))
        self.assertEqual(hgt.orf_overlap(100, 199, 150, 249), (50, 0.5, "partial_overlap"))
        self.assertEqual(hgt.orf_overlap(100, 199, 110, 189), (80, 1.0, "contained"))
        self.assertEqual(hgt.orf_overlap(100, 199, 100, 199), (100, 1.0, "contained"))
        self.assertLess(hgt.orf_overlap(100, 199, 99, 199)[1], 1)

    def test_redundancy_requires_containment_on_both_endpoints(self):
        base = dict(contig_a="A", contig_b="B", strand_b="+", start_a=1, end_a=1000,
                    start_b=1, end_b=1000, alignment_length=1000, bitscore=100)
        parent = dict(base, fragment_id="HGT1")
        inner = dict(base, fragment_id="HGT2", start_a=4, start_b=4, alignment_length=1200)
        different = dict(base, fragment_id="HGT3", start_a=500, end_a=1500)
        retained, redundant = hgt.remove_contained_fragments([inner, different, parent])
        self.assertEqual({r["fragment_id"] for r in retained}, {"HGT1", "HGT3"})
        self.assertEqual(redundant, [dict(discarded_fragment_id="HGT2", retained_fragment_id="HGT1")])

    def test_species_or_ani_rule_and_unknown_species(self):
        self.ani.write_text("A\tB\t95\n")
        tax = self.base / "tax.tsv"
        tax.write_text("mag_id\tspecies\nA\ts__Species A\nB\ts__Species B\n")
        out = self.run_ok(extra=["--taxonomy", str(tax)])
        self.assertEqual(rows(out / "01_blastn/selected_fragment_pairs.tsv")[0]["evidence"], "different_species")
        tax.write_text("mag_id\tspecies\nA\tunknown\nB\ts__Species B\n")
        out = self.run_ok("unknown", ["--taxonomy", str(tax)])
        self.assertEqual(rows(out / "01_blastn/selected_fragment_pairs.tsv"), [])

    def test_conservative_directional_ani(self):
        self.ani.write_text("A\tB\t94\nB\tA\t96\n")
        out = self.run_ok()
        self.assertEqual(rows(out / "01_blastn/selected_fragment_pairs.tsv"), [])

    def test_length_500_included_identity_99_and_ani_95_excluded(self):
        self.blast.write_text(blast("ctgA", "ctgB", 1, 500, 1, 500, length=500))
        out = self.run_ok()
        self.assertEqual(len(rows(out / "01_blastn/selected_fragment_pairs.tsv")), 1)
        self.ani.write_text("A\tB\t95\n")
        out = self.run_ok("ani_boundary")
        self.assertEqual(rows(out / "01_blastn/selected_fragment_pairs.tsv"), [])

    def test_gff_import_and_arg_best_hit(self):
        gff = self.base / "genes.gff"
        gff.write_text("##gff-version 3\nctgA\tProdigal\tCDS\t120\t419\t.\t+\t0\tID=argA;partial=00\n")
        self.hits.write_text("argA\tref1\t80\t75\t0\t0\t1\t75\t1\t75\t1e-10\t200\t100\n"
                             "argA\tref2\t85\t100\t0\t0\t1\t100\t1\t100\t1e-30\t190\t100\n")
        out = self.run_ok(extra=["--arg-hits", str(self.hits), "--orf-gff", str(gff)])
        self.assertEqual(rows(out / "04_arg_annotation/fragment_arg_details.tsv")[0]["subject_id"], "ref2")

    def test_v2_defaults_and_custom_arg_thresholds(self):
        defaults = hgt.parser().parse_args(self.args())
        self.assertEqual(defaults.annotator, "blastp")
        self.assertEqual((defaults.arg_identity, defaults.arg_query_cover, defaults.arg_evalue),
                         (80, 75, 1e-10))
        self.assertEqual(defaults.min_fragment_length, 500)
        # A boundary hit passes the defaults, then each user override independently rejects it.
        self.hits.write_text("argA\tref1\t80\t75\t0\t0\t1\t75\t1\t75\t1e-10\t200\t100\n")
        options = ["--arg-hits", str(self.hits), "--orf-tsv", str(self.orfs)]
        out = self.run_ok("default_boundary", options)
        self.assertEqual(len(rows(out / "04_arg_annotation/fragment_arg_details.tsv")), 1)
        for flag, value in (("--arg-identity", "81"), ("--arg-query-cover", "76"),
                            ("--arg-evalue", "1e-11")):
            out = self.run_ok(flag.lstrip('-'), options + [flag, value])
            self.assertEqual(rows(out / "04_arg_annotation/fragment_arg_details.tsv"), [])

    def test_gtdb_species_with_missing_ani(self):
        self.ani.write_text("")
        tax = self.base / "gtdb.tsv"
        tax.write_text("user_genome\tclassification\nA\td__Bacteria;p__X;s__Species A\n"
                       "B\td__Bacteria;p__X;s__Species B\n")
        out = self.run_ok(extra=["--taxonomy", str(tax)])
        row = rows(out / "01_blastn/selected_fragment_pairs.tsv")[0]
        self.assertEqual(row["ani"], "NA")
        self.assertEqual(row["evidence"], "different_species")

    def test_generated_blast_and_ani_commands_with_mocked_tools(self):
        calls = []
        def fake(command, work):
            calls.append(command)
            if command[0] == "blastn":
                dest = Path(command[command.index("-out") + 1])
                query = Path(command[command.index("-query") + 1]).stem
                dest.write_text(blast("ctgA", "ctgB", 101, 700, 800, 201) if query == "A"
                                else blast("ctgB", "ctgA", 800, 201, 101, 700))
            elif command[0] == "fastANI":
                Path(command[command.index("-o") + 1]).write_text("A\tB\t94\n")
        with patch.object(hgt, "run_command", side_effect=fake):
            self.assertEqual(hgt.main(["-i", str(self.mags), "-o", str(self.base / "out"),
                                       "--taxonomy", str(self.taxonomy)]), 0)
        self.assertEqual(len(rows(self.base / "out" / "01_blastn/selected_fragment_pairs.tsv")), 1)
        self.assertEqual(sum(c[0] == "blastn" for c in calls), 2)

    def test_stage_outputs_preserve_imports_and_filtered_best_hits(self):
        self.hits.write_text("argA\tref1\t90\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\t100\n"
                             "argA\tref2\t90\t100\t0\t0\t1\t100\t1\t100\t1e-30\t210\t100\n")
        final = self.run_ok(extra=["--arg-hits", str(self.hits), "--orf-tsv", str(self.orfs)])
        out = final
        self.assertEqual((out / "01_blastn/blastn.tsv").read_bytes(), self.blast.read_bytes())
        self.assertEqual((out / "02_fastani/fastani.tsv").read_bytes(), self.ani.read_bytes())
        self.assertEqual((out / "04_arg_annotation/arg_hits.tsv").read_bytes(), self.hits.read_bytes())
        self.assertTrue((out / "03_prodigal/orf_coordinates.tsv").is_file())
        self.assertEqual(len(rows(out / "04_arg_annotation/arg_hits.filtered.tsv")), 2)
        best = rows(out / "04_arg_annotation/arg_best_hits.tsv")
        self.assertEqual(len(best), 1)
        self.assertEqual(best[0]["subject_id"], "ref2")
        log = (out / "pipeline.log").read_text()
        self.assertIn("Reusing provided BLASTn results", log)
        self.assertIn("[05] DONE", log)
        self.assertTrue((out / "output_files.tsv").is_file())
        guide = (final / "05_summary/README.md").read_text()
        self.assertIn("Column groups", guide)
        self.assertIn("HGT_summary.tsv", guide)

    def test_tool_progress_records_completion_and_failure(self):
        work = self.base / "01_blastn"
        work.mkdir()
        with patch.object(hgt.shutil, "which", return_value="/mock/blastn"), \
             patch.object(hgt.subprocess, "run"), \
             patch("sys.stdout", new_callable=io.StringIO) as terminal:
            hgt.run_command(["blastn", "-version"], work)
        self.assertEqual(terminal.getvalue(), "")
        self.assertTrue((work / "logs/commands.jsonl").is_file())
        self.assertTrue((work / "logs/tools.log").is_file())
        self.assertFalse((work / "tools.log").exists())
        self.assertIn("DONE blastn", (self.base / "pipeline.log").read_text())
        with patch.object(hgt.shutil, "which", return_value="/mock/blastn"), \
             patch.object(hgt.subprocess, "run", side_effect=subprocess.CalledProcessError(2, ["blastn"])):
            with self.assertRaises(subprocess.CalledProcessError):
                hgt.run_command(["blastn", "-version"], work)
        self.assertIn("FAILED blastn", (self.base / "pipeline.log").read_text())

    def test_standard_twelve_column_hits_need_coverage_or_explicit_trust(self):
        self.hits.write_text("argA\tref1\t90\t100\t0\t0\t1\t100\t1\t100\t1e-30\t200\n")
        extra = ["--arg-hits", str(self.hits), "--orf-tsv", str(self.orfs)]
        self.assertEqual(hgt.main(self.args(extra=extra)), 1)
        out = self.run_ok("trusted", extra + ["--arg-hits-filtered"])
        self.assertEqual(len(rows(out / "04_arg_annotation/fragment_arg_details.tsv")), 1)
        self.assertEqual(rows(out / "04_arg_annotation/fragment_arg_details.tsv")[0]["arg_query_cover"], "NA")

    def test_coverage_from_protein_fasta_and_reject_low_coverage(self):
        self.hits.write_text("argA\tref1\t90\t50\t0\t0\t1\t50\t1\t50\t1e-30\t200\n")
        proteins = self.base / "genes.faa"
        proteins.write_text(">argA\n" + "M" * 100 + "\n")
        out = self.run_ok(extra=["--arg-hits", str(self.hits), "--orf-tsv", str(self.orfs),
                                "--orf-proteins", str(proteins)])
        self.assertEqual(rows(out / "04_arg_annotation/fragment_arg_details.tsv"), [])

    def test_mismatched_orf_ids_and_invalid_coordinates_rejected(self):
        self.hits.write_text("wrong\tref1\t90\t100\t0\t0\t1\t100\t1\t100\t1e-30\t200\t100\n")
        self.assertEqual(hgt.main(self.args(extra=["--arg-hits", str(self.hits),
                                                   "--orf-tsv", str(self.orfs)])), 1)
        self.blast.write_text(blast("ctgA", "ctgB", 1, 2500, 1, 600))
        self.assertEqual(hgt.main(self.args("badcoord")), 1)

    def test_duplicate_contig_ids_rejected(self):
        (self.mags / "B.fa").write_text(">ctgA\n" + self.seq_b + "\n")
        self.assertEqual(hgt.main(self.args()), 1)

    def test_output_reuse_rejected_without_overwrite(self):
        out = self.run_ok()
        before = (out / "01_blastn/selected_fragment_pairs.tsv").read_bytes()
        self.assertEqual(hgt.main(self.args()), 1)
        self.assertEqual(before, (out / "01_blastn/selected_fragment_pairs.tsv").read_bytes())

    def test_prodigal_gff_id_mapping(self):
        work = self.base / "work"
        work.mkdir()
        def fake(command, _):
            self.assertEqual(command[0], "prodigal")
            (work / "genes.gff").write_text("ctgA\tProdigal\tCDS\t120\t419\t.\t+\t0\tID=1_1;partial=00\n")
            (work / "genes.faa").write_text(">ctgA_1 # 120 # 419 # 1 # ID=1_1;partial=00;\n" + "M" * 99 + "*\n")
        with patch.object(hgt, "run_command", side_effect=fake):
            orfs, proteins = hgt.predict_orfs({"ctgA": self.seq_a}, work)
        self.assertEqual(orfs["ctgA_1"]["start"], 120)

    def test_imported_prodigal_gff_and_protein_header_mapping(self):
        gff = self.base / "prodigal.gff"
        gff.write_text("ctgA\tProdigal\tCDS\t120\t419\t.\t+\t0\tID=1_1;partial=00\n")
        proteins = self.base / "prodigal.faa"
        proteins.write_text(">ctgA_1 # 120 # 419 # 1 # ID=1_1;partial=00;\n" + "M" * 100 + "\n")
        self.hits.write_text("ctgA_1\tref1\t90\t100\t0\t0\t1\t100\t1\t100\t1e-30\t200\n")
        out = self.run_ok(extra=["--arg-hits", str(self.hits), "--orf-gff", str(gff),
                                "--orf-proteins", str(proteins)])
        self.assertEqual(rows(out / "04_arg_annotation/fragment_arg_details.tsv")[0]["orf_id"], "ctgA_1")

    def test_auto_annotation_branches_with_mocked_external_tools(self):
        db = self.base / "args.faa"
        db.write_text(">ref1\n" + "M" * 100 + "\n")
        for annotator in (None, "diamond", "blastp"):
            calls = []
            def fake(command, work):
                calls.append(command)
                if command[0] == "prodigal":
                    (work / "genes.gff").write_text("ctgA\tProdigal\tCDS\t120\t419\t.\t+\t0\tID=1_1;partial=00\n")
                    (work / "genes.faa").write_text(">ctgA_1 # 120 # 419 # 1 # ID=1_1;partial=00;\n" + "M" * 100 + "\n")
                elif command[:2] == ["diamond", "blastp"] or command[0] == "blastp":
                    (work / "arg_hits.tsv").write_text("ctgA_1\tref1\t90\t100\t0\t0\t1\t100\t1\t100\t1e-30\t200\t100\n")
            with patch.object(hgt, "run_command", side_effect=fake):
                options = ["--arg-db", str(db)]
                if annotator:
                    options += ["--annotator", annotator]
                out = self.run_ok(annotator or "default_blastp", options)
            self.assertEqual(rows(out / "04_arg_annotation/fragment_arg_details.tsv")[0]["orf_id"], "ctgA_1")
            arg_dir = out / "04_arg_annotation"
            self.assertTrue((arg_dir / "database").is_dir())
            expected_prefix = str(arg_dir / "database/arg_db")
            index_calls = [c for c in calls if c[0] == "makeblastdb" or c[:2] == ["diamond", "makedb"]]
            self.assertTrue(any(expected_prefix in c for c in index_calls))
            if annotator == "diamond":
                self.assertTrue(any(c[:2] == ["diamond", "blastp"] for c in calls))
            else:
                self.assertTrue(any(c[0] == "blastp" for c in calls))
                self.assertFalse(any(c[0] == "diamond" for c in calls))

    def test_external_failure_stops_run(self):
        with patch.object(hgt, "run_command", side_effect=subprocess.CalledProcessError(1, ["makeblastdb"])):
            args = ["-i", str(self.mags), "-o", str(self.base / "out"),
                    "--taxonomy", str(self.taxonomy)]
            self.assertEqual(hgt.main(args), 1)
        self.assertFalse((self.base / "out" / "run_summary.json").exists())

    def test_cli_help_and_imported_orf_requirement(self):
        proc = subprocess.run([sys.executable, str(ROOT / "hgt_demo.py"), "--help"], capture_output=True)
        self.assertEqual(proc.returncode, 0)
        proc = subprocess.run([sys.executable, str(ROOT / "hgt_demo.py")]
                              + self.args(extra=["--arg-hits", str(self.hits)]), capture_output=True)
        self.assertEqual(proc.returncode, 2)

    def test_taxonomy_is_required_and_must_cover_all_mags(self):
        proc = subprocess.run([sys.executable, str(ROOT / "hgt_demo.py"), "-i", str(self.mags),
                               "-o", str(self.base / "out")], capture_output=True)
        self.assertEqual(proc.returncode, 2)
        self.assertIn(b"--taxonomy", proc.stderr)
        self.taxonomy.write_text("mag_id\tspecies\nA\ts__Species A\n")
        self.assertEqual(hgt.main(self.args()), 1)
        self.assertFalse((self.base / "out").exists())

    def test_ani_is_calculated_when_not_provided(self):
        calls = []
        def fake(command, work):
            calls.append(command)
            self.assertEqual(command[0], "fastANI")
            Path(command[command.index("-o") + 1]).write_text("A\tB\t94\n")
        args = ["-i", str(self.mags), "-o", str(self.base / "out"),
                "--taxonomy", str(self.taxonomy), "--blast-results", str(self.blast)]
        with patch.object(hgt, "run_command", side_effect=fake):
            self.assertEqual(hgt.main(args), 0)
        report = json.loads((self.base / "out" / "run_summary.json").read_text())
        self.assertEqual(report["ani_source"], "computed_with_fastani")
        self.assertEqual(len(calls), 1)

    def test_supplied_ani_is_used_without_fastani_call(self):
        with patch.object(hgt, "generate_ani", side_effect=AssertionError("Unexpected ANI recalculation")):
            out = self.run_ok()
        report = json.loads((out / "run_summary.json").read_text())
        self.assertEqual(report["ani_source"], "provided")


if __name__ == "__main__":
    unittest.main()
