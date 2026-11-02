"""Verify documented optional-input modes without real external tools."""
import csv
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASH = os.environ.get('HGT_BASH') or shutil.which('bash')


def table(path):
    with path.open(encoding='utf-8-sig') as f:
        return list(csv.DictReader(f, delimiter='\t'))


@unittest.skipUnless(BASH and (ROOT/'testdata').is_dir(), 'Bash and testdata are needed for scenario checks')
class UsageScenarioTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='HGT usage ')
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.config = self.base / 'case_config.sh'
        example = ROOT / 'testdata'
        self.settings = dict(MAG_DIR=example/'mags', TAXONOMY_TSV=example/'taxonomy.tsv',
                             SARG_FASTA=example/'genes.faa', ORF_GFF=example/'genes.gff',
                             ORF_TSV=example/'orfs.tsv', ORF_PROTEINS=example/'genes.faa',
                             ARG_HITS=example/'arg_hits.tsv', FILTERED_ARG_HITS=example/'arg_hits.prefiltered.tsv',
                             ANI_RESULTS=example/'fastani.tsv', BLASTN_RESULTS=example/'blastn.tsv',
                             ARG_MAP=example/'arg_map.tsv', OUTPUT_BASE=self.base/'results',
                             THREADS='2', RUNNER='python', PYTHON_BIN=Path(sys.executable))
        self.write_config()

    def write_config(self):
        self.config.write_text(''.join(f'{key}={shlex.quote(value.as_posix() if isinstance(value, Path) else value)}\n'
                                      for key, value in self.settings.items()), encoding='utf-8')

    def launch(self, case, dry=False):
        command = [BASH, str(ROOT/'usage_scenarios/run_case.sh'), case, '--config', str(self.config)]
        if dry:
            command.append('--dry-run')
        return subprocess.run(command, capture_output=True, text=True, encoding='utf-8',
                              env=dict(os.environ, PYTHONIOENCODING='utf-8'), check=True)

    def test_all_scenario_commands(self):
        flags = {}
        for case in ('auto', 'orfs', 'arg_hits', 'ani', 'blastn', 'reuse_all', 'reuse_tsv', 'filtered_hits', 'hgt_only'):
            flags[case] = set(shlex.split(self.launch(case, True).stdout))
        self.assertIn('--arg-db', flags['auto'])
        self.assertNotIn('--orf-gff', flags['auto'])
        self.assertIn('--orf-proteins', flags['orfs'])
        self.assertIn('--arg-db', flags['orfs'])
        self.assertIn('--arg-hits', flags['arg_hits'])
        self.assertNotIn('--arg-db', flags['arg_hits'])
        self.assertNotIn('--ani-results', flags['arg_hits'])
        self.assertNotIn('--blast-results', flags['arg_hits'])
        self.assertIn('--ani-results', flags['ani'])
        self.assertIn('--blast-results', flags['blastn'])
        self.assertTrue({'--ani-results', '--blast-results', '--orf-gff', '--arg-hits'} <= flags['reuse_all'])
        self.assertIn('--orf-tsv', flags['reuse_tsv'])
        self.assertNotIn('--orf-gff', flags['reuse_tsv'])
        self.assertIn('--arg-hits-filtered', flags['filtered_hits'])
        self.assertFalse({'--arg-db', '--arg-hits', '--arg-map'} & flags['hgt_only'])

    def test_import_modes_export_and_roundtrip(self):
        self.launch('reuse_all')
        source = self.base/'results/reuse_all'
        expected = table(source/'05_summary/HGT_summary.tsv')
        self.assertEqual(len(expected), 2)
        prepared = self.base/'prepared'
        subprocess.run([sys.executable, str(ROOT/'usage_scenarios/prepare_existing_inputs.py'),
                        '--run-dir', str(source), '-o', str(prepared)], check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual((prepared/'arg/arg_hits.tsv').read_bytes(), (ROOT/'testdata/arg_hits.tsv').read_bytes())
        self.assertEqual((prepared/'prodigal/genes.gff').read_bytes(), (ROOT/'testdata/genes.gff').read_bytes())
        self.assertEqual(len(table(prepared/'prodigal/orf_coordinates.tsv')), 4)
        record = json.loads((prepared/'preparation.json').read_text())
        self.assertEqual(record['prefiltered_hit_rows'], 3)
        self.settings.update(ORF_GFF=prepared/'prodigal/genes.gff', ORF_TSV=prepared/'prodigal/orf_coordinates.tsv',
                             ORF_PROTEINS=prepared/'prodigal/genes.faa', ARG_HITS=prepared/'arg/arg_hits.tsv',
                             FILTERED_ARG_HITS=prepared/'arg/arg_hits.prefiltered.tsv', ANI_RESULTS=prepared/'ani/fastani.tsv',
                             BLASTN_RESULTS=prepared/'blastn/blastn.tsv')
        self.write_config()
        for case in ('reuse_tsv', 'filtered_hits'):
            self.launch(case)
            self.assertEqual(table(self.base/f'results/{case}/05_summary/HGT_summary.tsv'), expected)
        self.settings.update(ARG_HITS=ROOT/'testdata/arg_hits.12col.tsv', OUTPUT_BASE=self.base/'results12')
        self.write_config()
        self.launch('reuse_all')
        self.assertEqual(table(self.base/'results12/reuse_all/05_summary/HGT_summary.tsv'), expected)

    def test_help_and_invalid_case(self):
        result = subprocess.run([BASH, str(ROOT/'usage_scenarios/run_case.sh'), '-h'], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0)
        self.assertIn('reuse_all', result.stdout)
        invalid = subprocess.run([BASH, str(ROOT/'usage_scenarios/run_case.sh'), 'unknown', '--config', str(self.config)],
                                 capture_output=True, text=True)
        self.assertEqual(invalid.returncode, 2)


if __name__ == '__main__':
    unittest.main()
