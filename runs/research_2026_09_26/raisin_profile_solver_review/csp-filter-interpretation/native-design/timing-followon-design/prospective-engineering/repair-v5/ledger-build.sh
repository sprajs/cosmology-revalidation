#!/usr/bin/env bash
set -euo pipefail
task_root=/home/szymon/Documents/ChatGPT/supernova
task_build="$task_root/phase2/pte/repair-v5/ledger-build"
export SNANA_DIR="$task_build"
export GSL_DIR="$task_root/phase2/official/build/sysroot/usr"
export CFITSIO_DIR=/usr
export PATH="$GSL_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$GSL_DIR/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset CERN_DIR ROOT_DIR SNANA_PYTHON_DIR
cd "$task_build/src"
task_compiler="gfortran -B$GSL_DIR/lib/gcc/x86_64-pc-linux-gnu/16/ -B/usr/lib/gcc/x86_64-pc-linux-gnu/16/"
make -j1 ../bin/snlc_sim.exe FFC="$task_compiler" LGSL="-L$GSL_DIR/lib -lgsl -lgslcblas" EXTRA_FLAGS_C="-O1 -fcommon -std=gnu99 -Wno-error=implicit-function-declaration"
