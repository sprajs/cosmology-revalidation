#!/usr/bin/env bash
# Separate output-only executable: preserve the baseline executable and provenance.
set -euo pipefail
task_root="$(cd "$(dirname "$0")/../../.." && pwd)"
task_label="${1:-audit-v3}"
case "$task_label" in audit|audit-v2|audit-v3) ;; *) exit 2 ;; esac
export SNANA_DIR="$task_root/phase2/official/build/SNANA-$task_label"
export GSL_DIR="$task_root/phase2/official/build/sysroot/usr"
export LD_LIBRARY_PATH="$GSL_DIR/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
if [[ ! -d "$SNANA_DIR" ]]; then
  cp -a "$task_root/phase2/official/build/SNANA-current" "$SNANA_DIR"
  patch -d "$SNANA_DIR" -p1 < "$task_root/phase2/official/inputs/snana-$task_label-output-only.patch"
fi
make -C "$SNANA_DIR/src" -j4 snlc_fit
