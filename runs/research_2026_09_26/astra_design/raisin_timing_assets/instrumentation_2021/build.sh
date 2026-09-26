#!/usr/bin/env bash
set -euo pipefail
task_root=/home/szymon/Documents/ChatGPT/supernova
task_build="$task_root/runs/research_2026_09_26/astra_design/raisin_timing_assets/instrumentation_2021/build"
export SNANA_DIR="$task_build"
export GSL_DIR="$task_root/phase2/official/build/sysroot/usr"
export CFITSIO_DIR=/usr
export PATH="$GSL_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$GSL_DIR/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset CERN_DIR ROOT_DIR SNANA_PYTHON_DIR
cd "$task_build/src"
make set_Cpreproc_flags makeDirs
task_compiler="gfortran -B$GSL_DIR/lib/gcc/x86_64-pc-linux-gnu/16/ -B/usr/lib/gcc/x86_64-pc-linux-gnu/16/"
make ../bin/fcasplit FFC="$task_compiler"
make -j2 ../bin/snlc_fit.exe FFC="$task_compiler" \
  LGSL="-L$GSL_DIR/lib -lgsl -lgslcblas" \
  EXTRA_FLAGS_C="-O1 -fcommon -std=gnu99 -Wno-error=implicit-function-declaration"
