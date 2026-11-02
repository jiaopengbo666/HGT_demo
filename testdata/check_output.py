#!/usr/bin/env python3
"""Compare example main table and FASTA exports with the supplied expected outputs."""
import argparse
import csv
from pathlib import Path

def table(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="Example run output root")
    args = parser.parse_args()
    expected = Path(__file__).resolve().parent / "expected"
    actual = args.output / "05_summary"
    try:
        if table(actual / "HGT_summary.tsv") != table(expected / "HGT_summary.tsv"):
            parser.exit(1, "FAIL: HGT_summary.tsv differs from the expected example.\n")
        for name in ("hgt_fragments.fna", "hgt_arg_orfs.fna"):
            if (actual / name).read_bytes() != (expected / name).read_bytes():
                parser.exit(1, f"FAIL: {name} differs from the expected example.\n")
    except OSError as error:
        parser.exit(1, f"FAIL: {error}\n")
    print("PASS: 2 fragment pairs / 2 summary rows; main table and FASTA exports match.")

if __name__ == "__main__":
    main()
