"""Reusable BLASTp and DIAMOND ARG database indexes."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import hgt_demo as hgt


class ArgDatabaseCacheTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.base = Path(temporary.name)
        self.database = self.base / "SARG.faa"
        self.database.write_text(">ref1\nMAAA\n")
        self.proteins = self.base / "genes.faa"
        self.proteins.write_text(">orf1\nMAAA\n")
        self.calls = []
        self.arguments = SimpleNamespace(arg_db=self.database, arg_db_cache=self.base / "cache",
                                         rebuild_arg_db=False, annotator="blastp", threads=2, arg_evalue=1e-10)

    def fake_tool(self, command, work):
        self.calls.append(command)
        if command[0] == "makeblastdb":
            prefix = command[command.index("-out") + 1]
            for suffix in (".phr", ".pin", ".psq"):
                Path(prefix + suffix).write_bytes(b"index")
        elif command[:2] == ["diamond", "makedb"]:
            Path(command[command.index("--db") + 1] + ".dmnd").write_bytes(b"index")
        else:
            (work / "arg_hits.tsv").write_text("orf1\tref1\t100\t4\t0\t0\t1\t4\t1\t4\t1e-20\t50\t4\n")

    def run_annotation(self, name):
        work = self.base / name
        work.mkdir()
        with patch.object(hgt, "run_command", side_effect=self.fake_tool):
            result = hgt.annotate(self.proteins, self.arguments, work)
        self.assertTrue(result.is_file())
        return json.loads((work / "database.json").read_text())

    def builders(self):
        return [call for call in self.calls if call[0] == "makeblastdb" or call[:2] == ["diamond", "makedb"]]

    def test_reuses_index_across_output_directories_for_both_aligners(self):
        for annotator in ("blastp", "diamond"):
            with self.subTest(annotator=annotator):
                self.calls.clear()
                self.arguments.annotator = annotator
                first = self.run_annotation(annotator + "_first")
                second = self.run_annotation(annotator + "_second")
                self.assertFalse(first["reused"])
                self.assertTrue(second["reused"])
                self.assertEqual(first["index_prefix"], second["index_prefix"])
                self.assertEqual(len(self.builders()), 1)
                self.assertEqual(len(self.calls), 3)

    def test_identical_content_at_another_path_reuses_index(self):
        first = self.run_annotation("first")
        relocated = self.base / "renamed.faa"
        relocated.write_bytes(self.database.read_bytes())
        self.arguments.arg_db = relocated
        second = self.run_annotation("second")
        self.assertTrue(second["reused"])
        self.assertEqual(first["index_prefix"], second["index_prefix"])
        self.assertEqual(len(self.builders()), 1)
        self.assertEqual(second["database_fasta"], str(relocated))

    def test_fasta_content_change_builds_a_new_index(self):
        first = self.run_annotation("first")
        self.database.write_text(">ref1\nMAAAA\n")
        second = self.run_annotation("changed")
        self.assertFalse(second["reused"])
        self.assertNotEqual(first["source_sha256"], second["source_sha256"])
        self.assertEqual(len(self.builders()), 2)

    def test_aligner_indexes_are_separate(self):
        first = self.run_annotation("blastp")
        self.arguments.annotator = "diamond"
        second = self.run_annotation("diamond")
        self.assertFalse(second["reused"])
        self.assertNotEqual(first["index_prefix"], second["index_prefix"])
        self.assertEqual(len(self.builders()), 2)

    def test_alignment_settings_do_not_rebuild_database(self):
        self.run_annotation("first")
        self.arguments.arg_evalue = 1e-20
        self.arguments.threads = 8
        second = self.run_annotation("stricter")
        self.assertTrue(second["reused"])
        self.assertEqual(len(self.builders()), 1)
        self.assertIn("1e-20", self.calls[-1])
        self.assertIn("8", self.calls[-1])

    def test_missing_index_file_triggers_rebuild(self):
        first = self.run_annotation("first")
        Path(first["index_prefix"] + ".psq").unlink()
        second = self.run_annotation("repaired")
        self.assertFalse(second["reused"])
        self.assertNotEqual(first["index_prefix"], second["index_prefix"])
        self.assertEqual(len(self.builders()), 2)

    def test_failed_build_is_not_reused(self):
        work = self.base / "failed"
        work.mkdir()
        def fail(command, _):
            Path(command[command.index("-out") + 1] + ".pin").write_bytes(b"partial")
            raise subprocess.CalledProcessError(1, command)
        with patch.object(hgt, "run_command", side_effect=fail):
            with self.assertRaises(subprocess.CalledProcessError):
                hgt.annotate(self.proteins, self.arguments, work)
        self.assertEqual(list(self.arguments.arg_db_cache.rglob("manifest.json")), [])
        second = self.run_annotation("retry")
        self.assertFalse(second["reused"])

    def test_force_rebuild_and_incomplete_builder_output(self):
        first = self.run_annotation("first")
        self.arguments.rebuild_arg_db = True
        second = self.run_annotation("forced")
        self.assertFalse(second["reused"])
        self.assertNotEqual(first["index_prefix"], second["index_prefix"])
        work = self.base / "empty"
        work.mkdir()
        with patch.object(hgt, "run_command"):
            with self.assertRaisesRegex(ValueError, "complete index"):
                hgt.annotate(self.proteins, self.arguments, work)
        self.arguments.rebuild_arg_db = False
        third = self.run_annotation("still_valid")
        self.assertTrue(third["reused"])
        self.assertEqual(second["index_prefix"], third["index_prefix"])


if __name__ == "__main__":
    unittest.main()
