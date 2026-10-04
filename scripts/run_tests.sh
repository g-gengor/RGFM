#!/usr/bin/env bash
# Regression tests for RGFM: run the bundled test cases and compare the
# results against the stored reference solutions.
#
#   scripts/run_tests.sh [all|grid|h]     (default: all)
#
#   grid  test_p1grid: displacement field on GRIDS/GRID_128 from the stored
#         H.txt, compared with US/U_truth.            (~10 s)
#   h     test_p1h: solves the dense 1116x1116 system for the H vector,
#         compared with H_truth.txt.                   (~12 min, 4 threads)
#
# run_RGFM reads every input relative to its working directory, so each case
# is copied to a scratch directory and run there; the repository is never
# modified. Environment overrides:
#   RGFM_BIN       path to the executable (default: <repo>/bin/run_RGFM)
#   RGFM_WORKDIR   where to put the scratch copies (default: mktemp -d)
#   GRID_RTOL / GRID_ATOL, H_RTOL / H_ATOL   comparison tolerances
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${RGFM_BIN:-$ROOT/bin/run_RGFM}"
COMPARE="$ROOT/scripts/compare_outputs.py"
WHICH="${1:-all}"

# U_truth is stored with 8 significant digits; H_truth with 16.
GRID_RTOL="${GRID_RTOL:-1e-6}"
GRID_ATOL="${GRID_ATOL:-1e-8}"
H_RTOL="${H_RTOL:-1e-6}"
H_ATOL="${H_ATOL:-1e-9}"

case "$WHICH" in
  all|grid|h) ;;
  *) echo "usage: $0 [all|grid|h]" >&2; exit 2 ;;
esac

if [[ ! -x "$BIN" ]]; then
  echo "error: $BIN not found; build first (cmake -S . -B build && cmake --build build)" >&2
  exit 2
fi
BIN="$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")"

WORKDIR="${RGFM_WORKDIR:-$(mktemp -d "${TMPDIR:-/tmp}/rgfm-tests.XXXXXX")}"
mkdir -p "$WORKDIR"
echo "RGFM binary: $BIN"
echo "Scratch dir: $WORKDIR"

failures=0

# run_case NAME OUTPUT ARGS... : copy test case NAME to the scratch dir and run
# it there. OUTPUT (the file to be checked) is deleted first, so a failed run
# cannot leave the copy committed in the repo behind to be compared.
run_case() {
  local name="$1" output="$2"; shift 2
  local dir="$WORKDIR/$name"
  rm -rf "$dir"
  cp -R "$ROOT/$name" "$dir"
  rm -f "$dir/$output"
  echo
  echo "== $name: running run_RGFM $* (log: $dir/run.log)"
  local start=$SECONDS
  if ! (cd "$dir" && "$BIN" "$@" > run.log 2>&1); then
    echo "FAIL: run_RGFM exited with an error; last lines of the log:"
    tail -n 20 "$dir/run.log"
    return 1
  fi
  echo "   finished in $((SECONDS - start)) s"
}

# check NAME RESULT REFERENCE RTOL ATOL
check() {
  local dir="$WORKDIR/$1"
  printf '   %s vs %s: ' "$2" "$3"
  python3 "$COMPARE" "$dir/$2" "$dir/$3" --rtol "$4" --atol "$5"
}

if [[ "$WHICH" == all || "$WHICH" == grid ]]; then
  if run_case test_p1grid US/U_128 128 \
     && check test_p1grid US/U_128 US/U_truth "$GRID_RTOL" "$GRID_ATOL"; then :; else
    failures=$((failures + 1))
  fi
fi

if [[ "$WHICH" == all || "$WHICH" == h ]]; then
  if run_case test_p1h H.txt \
     && check test_p1h H.txt H_truth.txt "$H_RTOL" "$H_ATOL"; then :; else
    failures=$((failures + 1))
  fi
fi

echo
if (( failures > 0 )); then
  echo "$failures test case(s) FAILED"
  exit 1
fi
echo "All requested test cases passed."
