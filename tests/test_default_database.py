"""Bundled SARG selection, metadata and cached automatic annotation."""
import csv
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import hgt_demo as hgt


class DefaultDatabaseTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.base = Path(temporary.name)
        self.mags = self.base / "mags"
        self.mags.mkdir()
        for name in ("A", "B"):
            (self.mags / f"{name}.fa").write_text(f">ctg{name}\n" + "ACGT" * 250 + "\n")
        self.taxonomy = self.base / "taxonomy.tsv"
        self.taxonomy.write_text("mag_id\tspecies\nA\tSpecies A\nB\tSpecies B\n")

    def args(self, name="out", extra=()):
        return ["-i", str(self.mags), "-o", str(self.base / name), "--taxonomy", str(self.taxonomy),
                "--arg-db-cache", str(self.base / "cache")] + list(extra)

    def test_bundled_fasta_matches_manifest_and_annotation_map(self):
        record = json.loads((ROOT / "database/SARG3.2.1-S.source.json").read_text())
        self.assertEqual(hashlib.sha256(hgt.DEFAULT_ARG_DB.read_bytes()).hexdigest(), record["sha256"])
        names = [name for name, _, _ in hgt.fasta(hgt.DEFAULT_ARG_DB)]
        self.assertEqual(len(names), record["sequence_count"])
        self.assertEqual(len(set(names)), len(names))
        with hgt.DEFAULT_ARG_MAP.open(encoding="utf-8") as handle:
            mapping = {row["subject_id"] for row in csv.DictReader(handle, delimiter="\t")}
        self.assertTrue(set(names) <= mapping)

    def test_default_uses_bundled_sarg(self):
        with patch.object(hgt, "analyze") as analyze:
            self.assertEqual(hgt.main(self.args()), 0)
        args = analyze.call_args.args[0]
        self.assertEqual(args.arg_db, hgt.DEFAULT_ARG_DB.resolve())
        self.assertEqual(args.arg_db_source, "bundled_SARG")
        self.assertEqual(args.arg_map, hgt.DEFAULT_ARG_MAP.resolve())

    def test_custom_database_overrides_default(self):
        custom = self.base / "custom.faa"
        custom.write_text(">custom\nMAAA\n")
        with patch.object(hgt, "analyze") as analyze:
            self.assertEqual(hgt.main(self.args(extra=["--arg-db", str(custom)])), 0)
        self.assertEqual(analyze.call_args.args[0].arg_db, custom.resolve())
        self.assertEqual(analyze.call_args.args[0].arg_db_source, "user")

    def test_imported_hits_skip_default_database(self):
        hits, orfs = self.base / "hits.tsv", self.base / "orfs.tsv"
        hits.write_text("")
        orfs.write_text("orf_id\tcontig_id\tstart\tend\tstrand\n")
        with patch.object(hgt, "analyze") as analyze:
            self.assertEqual(hgt.main(self.args(extra=["--arg-hits", str(hits), "--orf-tsv", str(orfs)])), 0)
        self.assertIsNone(analyze.call_args.args[0].arg_db)
        self.assertEqual(analyze.call_args.args[0].arg_db_source, "not_used")

    def test_hgt_only_works_without_bundled_database(self):
        with patch.object(hgt, "DEFAULT_ARG_DB", self.base / "missing.faa"), patch.object(hgt, "analyze") as analyze:
            self.assertEqual(hgt.main(self.args(extra=["--hgt-only"])), 0)
        self.assertIsNone(analyze.call_args.args[0].arg_db)
        self.assertIsNone(analyze.call_args.args[0].arg_map)

    def test_missing_default_database_reports_error(self):
        with patch.object(hgt, "DEFAULT_ARG_DB", self.base / "missing.faa"), patch.object(hgt, "analyze") as analyze:
            self.assertEqual(hgt.main(self.args()), 1)
        analyze.assert_not_called()

    def test_hgt_only_and_custom_annotation_are_exclusive(self):
        for extra in (["--hgt-only", "--arg-db", "custom.faa"], ["--hgt-only", "--arg-hits", "hits.tsv"]):
            with self.subTest(extra=extra), self.assertRaises(SystemExit) as error:
                hgt.parser().parse_args(self.args(extra=extra))
            self.assertEqual(error.exception.code, 2)

    def test_default_automatic_annotation_reuses_sarg_index(self):
        blast, ani = self.base / "blast.tsv", self.base / "ani.tsv"
        blast.write_text("ctgA\tctgB\t99.9\t600\t0\t0\t101\t700\t201\t800\t1e-100\t900\n")
        ani.write_text("A\tB\t94\t10\t12\n")
        calls = []
        def fake(command, work):
            calls.append(command)
            if command[0] == "prodigal":
                (work / "genes.gff").write_text(
                    "ctgA\tProdigal\tCDS\t120\t419\t.\t+\t0\tID=1_1;partial=00\n"
                    "ctgB\tProdigal\tCDS\t220\t519\t.\t+\t0\tID=2_1;partial=00\n")
                (work / "genes.faa").write_text(
                    ">ctgA_1 # 120 # 419 # 1 # ID=1_1;partial=00;\n" + "M" * 100 + "\n"
                    ">ctgB_1 # 220 # 519 # 1 # ID=2_1;partial=00;\n" + "M" * 100 + "\n")
            elif command[0] == "makeblastdb":
                self.assertEqual(Path(command[command.index("-in") + 1]), hgt.DEFAULT_ARG_DB.resolve())
                prefix = command[command.index("-out") + 1]
                for suffix in (".phr", ".pin", ".psq"):
                    Path(prefix + suffix).write_bytes(b"index")
            elif command[0] == "blastp":
                (work / "arg_hits.tsv").write_text("".join(
                    f"{oid}\tABA08387\t100\t100\t0\t0\t1\t100\t1\t100\t1e-40\t500\t100\n"
                    for oid in ("ctgA_1", "ctgB_1")))
        options = ["--blast-results", str(blast), "--ani-results", str(ani)]
        with patch.object(hgt, "run_command", side_effect=fake):
            for name in ("first", "second"):
                self.assertEqual(hgt.main(self.args(name, options)), 0)
        self.assertEqual(sum(call[0] == "makeblastdb" for call in calls), 1)
        self.assertEqual(sum(call[0] == "blastp" for call in calls), 2)
        first = json.loads((self.base / "first/04_arg_annotation/database.json").read_text())
        second = json.loads((self.base / "second/04_arg_annotation/database.json").read_text())
        self.assertFalse(first["reused"])
        self.assertTrue(second["reused"])
        with (self.base / "second/05_summary/HGT_summary.tsv").open() as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual((row["arg_a"], row["arg_b"], row["shared_complete_args"]), ("dfrA1", "dfrA1", "dfrA1"))


if __name__ == "__main__":
    unittest.main()
