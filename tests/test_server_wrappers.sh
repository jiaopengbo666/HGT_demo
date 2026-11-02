#!/usr/bin/env bash
# No real Conda/tool installation required; bounded mocks reproduce inherited-PATH recursion.
set -euo pipefail
root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
fixture="$(mktemp -d "$root/tests/.wrapper_test.XXXXXX")"
printf 'Wrapper regression fixture: %s\n' "$fixture"
# Temporary fixtures are intentionally retained for inspection; no recursive cleanup.
mkdir -p "$fixture/project" "$fixture/conda/bin"
cp "$root/run_server_existing_envs.sh" "$fixture/project/launcher.sh"
export HGT_CONDA_ROOT="$fixture/conda"
export HGT_MOCK_LOG="$fixture/calls.txt"
for env_name in arg_hgt_python blast diamond prodigal dRep; do
  mkdir -p "$HGT_CONDA_ROOT/envs/$env_name/bin"
done
cat > "$HGT_CONDA_ROOT/bin/conda" <<'MOCK_CONDA'
#!/usr/bin/env bash
set -euo pipefail
depth="${HGT_MOCK_DEPTH:-0}"
if (( depth >= 3 )); then echo 'BOUNDED_RECURSION_DETECTED' >&2; exit 86; fi
export HGT_MOCK_DEPTH=$((depth + 1))
[[ "$1" == run ]]; shift
[[ "$1" == --no-capture-output ]]; shift
[[ "$1" == -p ]]; shift 2
if [[ "$1" == python ]]; then printf 'Mock Python 3.12\n'; exit 0; fi
# Preserve PATH deliberately: this is the condition that exposed the real bug.
exec "$@"
MOCK_CONDA
chmod u+x "$HGT_CONDA_ROOT/bin/conda"
for pair in blast:blastn blast:blastp blast:makeblastdb diamond:diamond prodigal:prodigal dRep:fastANI; do
  env_name="${pair%%:*}"; tool="${pair##*:}"
  cat > "$HGT_CONDA_ROOT/envs/$env_name/bin/$tool" <<'MOCK_TOOL'
#!/usr/bin/env bash
printf '%s\n' "$(basename -- "$0")" >> "$HGT_MOCK_LOG"
printf '%s\n' "$@" >> "$HGT_MOCK_LOG"
printf 'Mock executable reached: %s\n' "$0"
MOCK_TOOL
  chmod u+x "$HGT_CONDA_ROOT/envs/$env_name/bin/$tool"
done

# Fixed launcher must reach real absolute mock executables even with wrappers first on PATH.
PATH="$fixture/project/.server_tools:$PATH" bash "$fixture/project/launcher.sh" --check
PATH="$fixture/project/.server_tools:$PATH" \
  "$fixture/project/.server_tools/makeblastdb" -in 'genome with spaces.fa'
grep -Fx 'genome with spaces.fa' "$HGT_MOCK_LOG" >/dev/null

# Reproduce the previous implementation, with a depth bound to avoid runaway processes.
sed 's/"$conda_exe" "$prefix" "$executable"/"$conda_exe" "$prefix" "$tool"/' \
  "$root/run_server_existing_envs.sh" > "$fixture/project/old_launcher.sh"
set +e
PATH="$fixture/project/.server_tools:$PATH" bash "$fixture/project/old_launcher.sh" --check \
  > "$fixture/old_stdout.txt" 2> "$fixture/old_stderr.txt"
status=$?
set -e
[[ "$status" == 86 ]] || { echo "Expected bounded recursion (86), got $status" >&2; exit 1; }
grep -F 'BOUNDED_RECURSION_DETECTED' "$fixture/old_stderr.txt" >/dev/null
printf 'PASS: former recursion reproduced safely; absolute-path fix and argument forwarding verified.\n'
