#!/usr/bin/env python3
"""Freeze the lightweight scientific record without copying bulk source data."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess

ROOT=Path(__file__).resolve().parents[1]
files=subprocess.check_output(['git','ls-files','--cached','--others','--exclude-standard','-z'],cwd=ROOT).decode().split('\0')
omit={'runs/final-manifest.json','runs/verification/record-check.json'}
records={}
for name in sorted(set(files)-omit-{''}):
    path=ROOT/name
    if path.is_file():records[name]=hashlib.sha256(path.read_bytes()).hexdigest()
verification=json.loads((ROOT/'runs/verification/record-check.json').read_text())
assert not verification['failures'] and not verification['broken_local_links']
replay=json.loads((ROOT/'runs/cosmology/replay-comparison.json').read_text())
assert len(replay)==19 and all(x['scientific_summary_identical'] for x in replay.values())
report={
    'created_utc':datetime.now(timezone.utc).isoformat(),
    'purpose':'Final local scientific record; no external publication',
    'code_revision_parent_at_freeze':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
    'branch':subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip(),
    'python':platform.python_version(),
    'configuration':{'source_freeze_date':'2026-09-20','scientific_report':'docs/final-report.md','reproducibility':'docs/reproducibility.md','goal_scope':'Completed bounded investigation with documented exact-reproduction barriers'},
    'verification_before_packaging':{'active_digest_records':verification['checked_digest_records'],'local_links':verification['checked_local_links'],'failed_hashes':0,'broken_links':0,'identical_scientific_cosmology_replays':19},
    'record_files_sha256':records,
    'bulk_retention':'Original data/papers/sources, large numerical chains and covariance arrays are local and Git-ignored. Exact dependencies are identified by active branch manifests. Copy the full workspace or reacquire exact hashed inputs; Git alone is incomplete.',
    'mutable_verifier_output':'runs/verification/record-check.json is excluded from packaging hashes to avoid a self-referential digest cycle; rerun scripts/verify_record.py to verify this package and its upstream manifests.'}
(ROOT/'runs/final-manifest.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'record_files':len(records),'verification':report['verification_before_packaging']}))
