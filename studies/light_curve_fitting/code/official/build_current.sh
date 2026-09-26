#!/usr/bin/env bash
# Rebuild the archived SNANA source with output-only diagnostics; no global install.
set -euo pipefail
task_root="$(cd "$(dirname "$0")/../../.." && pwd)"
task_out="$task_root/phase2/official"
export SNANA_DIR="$task_out/build/SNANA-current"
export GSL_DIR="$task_out/build/sysroot/usr"
export PATH="$GSL_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$GSL_DIR/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export FC="$GSL_DIR/bin/gfortran -B$GSL_DIR/lib/gcc/x86_64-pc-linux-gnu/16/ -B/usr/lib/gcc/x86_64-pc-linux-gnu/16/"
if [[ ! -d "$SNANA_DIR" ]]; then
    cp -a "$task_root/sources/repos/RickKessler__SNANA" "$SNANA_DIR"
    patch -d "$SNANA_DIR" -p1 < "$task_out/inputs/snana-output-only.patch"
fi
cd "$SNANA_DIR"
autoreconf -fi
./configure
mkdir -p obj bin
make -C src set_Cpreproc_flags
make -C src -j4 snlc_fit
