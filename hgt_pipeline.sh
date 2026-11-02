#!/usr/bin/env bash
# Command-line entry point; see README.md for inputs and examples.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec "${PYTHON:-python3}" "$script_dir/hgt_demo.py" "$@"
