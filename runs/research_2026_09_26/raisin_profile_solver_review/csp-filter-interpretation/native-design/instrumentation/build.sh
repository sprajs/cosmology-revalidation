#!/usr/bin/env bash
set -euo pipefail
task_root="/home/szymon/Documents/ChatGPT/supernova"
task_build="/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_profile_solver_review/csp-filter-interpretation/native-design/instrumentation/SNANA-v11_04k-output"
export SNANA_DIR="$task_build"
export GSL_DIR="$task_root/phase2/official/build/sysroot/usr"
export CFITSIO_DIR=/usr
export PATH="$GSL_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$GSL_DIR/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$SNANA_DIR/src"
task_compiler="gfortran -B$GSL_DIR/lib/gcc/x86_64-pc-linux-gnu/16/ -B/usr/lib/gcc/x86_64-pc-linux-gnu/16/"
make -j2 ../bin/snlc_fit.exe FFC="$task_compiler" \
  LGSL="-L$GSL_DIR/lib -lgsl -lgslcblas" \
  EXTRA_FLAGS_C="-O1 -fcommon -std=gnu99 -Wno-error=implicit-function-declaration"
