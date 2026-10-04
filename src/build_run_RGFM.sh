#!/usr/bin/env bash
# Thin convenience wrapper around CMake: configures a Release build in
# <repo>/build and puts the executable in <repo>/bin/run_RGFM.
# Extra arguments are passed to the configure step, e.g.
#   src/build_run_RGFM.sh -DCMAKE_PREFIX_PATH=/path/to/deps
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cmake -S "$ROOT" -B "$ROOT/build" -DCMAKE_BUILD_TYPE=Release "$@"
cmake --build "$ROOT/build" --parallel
