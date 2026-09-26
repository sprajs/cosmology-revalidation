#!/usr/bin/env bash
# Isolated historical NIR fitter build. Never installs into global paths.
set -euo pipefail
task_root="$(cd "$(dirname "$0")/../.." && pwd)"
build_root="$task_root/runs/research_2026_09_26/raisin_nir_timing_sensitivity/SNANA-v11_03c-source"
expected_commit=06f2ccfdbf99d62c23dec66d0bf7b7e9452604d2
test "$(git -C "$build_root" rev-parse HEAD)" = "$expected_commit"
export SNANA_DIR="$build_root"
export GSL_DIR="$task_root/phase2/official/build/sysroot/usr"
export CFITSIO_DIR=/usr
export PATH="$GSL_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$GSL_DIR/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$SNANA_DIR/src"
make set_Cpreproc_flags makeDirs
task_compiler="gfortran -B$GSL_DIR/lib/gcc/x86_64-pc-linux-gnu/16/ -B/usr/lib/gcc/x86_64-pc-linux-gnu/16/"
make ../bin/fcasplit FFC="$task_compiler"
make -j2 ../bin/snlc_fit.exe FFC="$task_compiler" \
  LGSL="-L$GSL_DIR/lib -lgsl -lgslcblas" \
  EXTRA_FLAGS_C="-O1 -fcommon -std=gnu99 -Wno-error=implicit-function-declaration"
