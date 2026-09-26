#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "$0")/../../.." && pwd)"
export SNANA_DIR="$task_root/phase2/official/build/SNANA-2fe0f56"
export GSL_DIR="$task_root/phase2/official/build/sysroot/usr"
export CFITSIO_DIR=/usr
export PATH="$GSL_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$GSL_DIR/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd "$SNANA_DIR/src"
make set_Cpreproc_flags makeDirs
compiler="gfortran -B$GSL_DIR/lib/gcc/x86_64-pc-linux-gnu/16/ -B/usr/lib/gcc/x86_64-pc-linux-gnu/16/"
make ../bin/fcasplit FFC="$compiler"
make -j4 ../bin/snlc_fit.exe ../bin/snlc_sim.exe FFC="$compiler" \
  LGSL="-L$GSL_DIR/lib -lgsl -lgslcblas" \
  EXTRA_FLAGS_C="-O1 -fcommon -std=gnu99 -Wno-error=implicit-function-declaration" \
  GIT_VERSION=2fe0f564
