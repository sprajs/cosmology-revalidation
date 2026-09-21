#!/usr/bin/env bash
set -euo pipefail
task_root="$(cd "$(dirname "$0")/../../.." && pwd)"
export SNANA_DIR="$task_root/phase2/official/build/SNANA-current"
export SNDATA_ROOT="$task_root/phase2/official/inputs/SNDATA_ROOT"
export LD_LIBRARY_PATH="$task_root/phase2/official/build/sysroot/usr/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
label="${1:-pilot12}"
"$SNANA_DIR/bin/snlc_fit.exe" "$task_root/phase2/official/inputs/snana_${label}.nml"
