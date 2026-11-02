#!/usr/bin/env bash
# Run tools from existing Conda environments.
set -euo pipefail

project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
conda_root="${HGT_CONDA_ROOT:-}"
if [[ -z "$conda_root" ]]; then
  command -v conda >/dev/null 2>&1 || { echo 'Set HGT_CONDA_ROOT or put conda on PATH.' >&2; exit 1; }
  conda_root="$(conda info --base)"
fi
conda_exe="${HGT_CONDA_EXE:-$conda_root/bin/conda}"
python_env="${HGT_PYTHON_ENV:-$conda_root/envs/arg_hgt_python}"
blast_env="${HGT_BLAST_ENV:-$conda_root/envs/blast}"
diamond_env="${HGT_DIAMOND_ENV:-$conda_root/envs/diamond}"
prodigal_env="${HGT_PRODIGAL_ENV:-$conda_root/envs/prodigal}"
fastani_env="${HGT_FASTANI_ENV:-$conda_root/envs/dRep}"
tools_dir="$project_dir/.server_tools"

[[ -x "$conda_exe" ]] || { echo "Conda executable not found: $conda_exe" >&2; exit 1; }
[[ -d "$python_env" ]] || { echo "Python environment not found: $python_env" >&2; exit 1; }
if [[ $# -eq 1 && ( "$1" == -h || "$1" == --help || "$1" == --version ) ]]; then
  exec "$conda_exe" run --no-capture-output -p "$python_env" \
    python "$project_dir/hgt_demo.py" "$1"
fi
mkdir -p "$tools_dir"

make_wrapper() {
  local tool="$1" prefix="$2"
  [[ -d "$prefix" ]] || { echo "Environment not found: $prefix" >&2; exit 1; }
  local executable="$prefix/bin/$tool"
  [[ -x "$executable" ]] || { echo "Tool executable not found: $executable" >&2; exit 1; }
  # Use absolute tool paths to avoid resolving the wrapper itself.
  printf '#!/usr/bin/env bash\nexec %q run --no-capture-output -p %q %q "$@"\n' \
    "$conda_exe" "$prefix" "$executable" > "$tools_dir/$tool"
  chmod u+x "$tools_dir/$tool"
}

make_wrapper makeblastdb "$blast_env"
make_wrapper blastn "$blast_env"
make_wrapper blastp "$blast_env"
make_wrapper diamond "$diamond_env"
make_wrapper prodigal "$prodigal_env"
make_wrapper fastANI "$fastani_env"

if [[ "${1:-}" == --check ]]; then
  "$conda_exe" run --no-capture-output -p "$python_env" python -c \
    'import sys; print("Python:", sys.version); sys.exit(0 if sys.version_info >= (3, 10) else 1)'
  "$tools_dir/blastn" -version
  "$tools_dir/makeblastdb" -version
  "$tools_dir/diamond" version
  "$tools_dir/prodigal" -v
  "$tools_dir/fastANI" --version
  printf '\nEnvironment check passed.\n'
  exit 0
fi

# Put tool wrappers first after entering the Python environment.
exec "$conda_exe" run --no-capture-output -p "$python_env" bash -c \
  'export PATH="$1:$PATH"; shift; exec python "$@"' \
  hgt-server "$tools_dir" "$project_dir/hgt_demo.py" "$@"
