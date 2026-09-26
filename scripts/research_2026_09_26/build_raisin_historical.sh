#!/usr/bin/env bash
# Isolated source-version-matched RAISIN fitter. No global installation.
set -euo pipefail
task_root="$(cd "$(dirname "$0")/../.." && pwd)"
export SNANA_DIR="$task_root/phase2/official/build/SNANA-v11_04k"
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
