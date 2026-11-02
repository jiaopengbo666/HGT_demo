#!/usr/bin/env python3
"""Export reusable inputs from a completed HGT_demo run."""
import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path


def protein_lengths(path):
    lengths, current = {}, None
    with path.open(encoding='utf-8-sig') as f:
        for line in f:
            if line.startswith('>'):
                current = line[1:].split()[0]
                if current in lengths:
                    raise ValueError('Duplicate protein ID: ' + current)
                lengths[current] = 0
            elif line.strip():
                if current is None:
                    raise ValueError('Protein FASTA sequence before header')
                lengths[current] += len(line.strip())
    return lengths


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''):
            value.update(chunk)
    return value.hexdigest()


def prepare(source, output):
    source = source.resolve()
    report = json.loads((source / 'run_summary.json').read_text(encoding='utf-8'))
    if output.exists() and any(output.iterdir()):
        raise ValueError('Output must be new or empty; files are not overwritten.')
    output.mkdir(parents=True, exist_ok=True)
    copied = []
    def copy_first(candidates, target):
        original = next((source / name for name in candidates if (source / name).is_file()), None)
        if original is None:
            return False
        destination = output / target
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(original, destination)
        checksum = digest(destination)
        if digest(original) != checksum:
            raise ValueError('Copied file hash mismatch: ' + str(destination))
        copied.append(dict(source=str(original), output=target,
                           sha256=checksum))
        return True
    copy_first(['00_inputs/taxonomy.tsv'], 'taxonomy.tsv')
    copy_first(['00_inputs/MAG_files.tsv'], 'MAG_files.tsv')
    copy_first(['01_blastn/blastn.tsv', 'work/blastn.tsv'], 'blastn/blastn.tsv')
    copy_first(['02_fastani/fastani.tsv', 'work/fastani.tsv'], 'ani/fastani.tsv')
    copy_first(['03_prodigal/genes.gff', '03_prodigal/provided_coordinates.gff',
                '03_prodigal/provided_coordinates.gff3', 'work/genes.gff'], 'prodigal/genes.gff')
    copy_first(['03_prodigal/genes.faa', 'work/genes.faa'], 'prodigal/genes.faa')
    copy_first(['03_prodigal/genes.fna', 'work/genes.fna'], 'prodigal/genes.fna')
    copy_first(['03_prodigal/orf_coordinates.tsv', '03_prodigal/provided_coordinates.tsv'], 'prodigal/orf_coordinates.tsv')
    copy_first(['04_arg_annotation/arg_hits.tsv', 'work/arg_hits.tsv'], 'arg/arg_hits.tsv')
    copy_first(['00_inputs/arg_map.tsv'], 'arg/arg_map.tsv')
    params = report.get('parameters', {})
    thresholds = {key: params.get(key, value) for key, value in
                  (('arg_identity', 80), ('arg_query_cover', 75), ('arg_evalue', 1e-10))}
    raw = output / 'arg/arg_hits.tsv'
    filtered_count = None
    if raw.is_file():
        proteins = output / 'prodigal/genes.faa'
        lengths = protein_lengths(proteins) if proteins.is_file() else {}
        selected_lines = []
        with raw.open(encoding='utf-8-sig') as f:
            for number, line in enumerate(f, 1):
                if not line.strip() or line.startswith('#'):
                    continue
                columns = line.rstrip('\r\n').split('\t')
                if len(columns) not in (12, 13):
                    raise ValueError(f'Raw ARG file row {number} is not standard outfmt 6')
                if not params.get('arg_hits_filtered', False):
                    qlen = int(columns[12]) if len(columns) == 13 else lengths.get(columns[0])
                    if not qlen or qlen <= 0:
                        raise ValueError('Coverage filtering needs qlen or matching protein FASTA')
                    cover = 100 * (abs(int(columns[7])-int(columns[6]))+1) / qlen
                    if (float(columns[2]) < thresholds['arg_identity'] or
                        float(columns[10]) > thresholds['arg_evalue'] or cover < thresholds['arg_query_cover']):
                        continue
                selected_lines.append(line.rstrip('\r\n') + '\n')
        (output / 'arg/arg_hits.prefiltered.tsv').write_text(''.join(selected_lines), encoding='utf-8')
        filtered_count = len(selected_lines)
    record = dict(source_run=str(source), original_parameters=params, copied_files=copied,
                  prefilter_thresholds=thresholds, prefiltered_hit_rows=filtered_count,
                  mag_directory=params.get('input'), note='MAG FASTA paths are recorded in preparation.json.')
    (output / 'preparation.json').write_text(json.dumps(record, indent=2)+'\n', encoding='utf-8')
    (output / 'README.md').write_text('''# Reusable inputs

Source paths, checksums and thresholds are recorded in preparation.json. MAGs remain at their original paths.

Use prodigal/genes.gff with genes.faa, or the normalized orf_coordinates.tsv. arg/arg_hits.tsv contains raw outfmt6 hits; arg_hits.prefiltered.tsv applies the source run thresholds for --arg-hits-filtered. ARG labels are in arg/arg_map.tsv when available.
''', encoding='utf-8')
    print(f'Prepared {len(copied)} source files in {output.resolve()}')
    print(f'Original MAG directory: {params.get("input", "unknown")}')
    for item in copied:
        print('  ' + item['output'])
    if filtered_count is not None:
        print(f'  arg/arg_hits.prefiltered.tsv: {filtered_count} passing hits')
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, required=True, help='Complete existing output directory, including stage folders')
    parser.add_argument('-o', '--output', type=Path, required=True, help='New/empty prepared-input directory')
    args = parser.parse_args()
    try:
        prepare(args.run_dir, args.output)
    except (OSError, ValueError) as error:
        parser.exit(1, f'ERROR: {error}\n')


if __name__ == '__main__':
    main()
