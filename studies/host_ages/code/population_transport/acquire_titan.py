#!/usr/bin/env python3
"""Restore pinned public author products; retain upstream code only as input."""
from pathlib import Path
import hashlib
import json
import urllib.request

ROOT = Path(__file__).resolve().parents[4]
manifest = ROOT / 'studies/host_ages/results/population_transport/author-inputs.json'
for record in json.loads(manifest.read_text())['files']:
    destination = ROOT / record['path']
    if destination.exists():
        raw = destination.read_bytes()
    else:
        raw = urllib.request.urlopen(record['url'], timeout=60).read()
    if hashlib.sha256(raw).hexdigest() != record['sha256']:
        raise RuntimeError(f'Input hash mismatch: {destination}')
    if len(raw) != record['bytes']:
        raise RuntimeError(f'Input size mismatch: {destination}')
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(raw)
    print(f'Verified {destination.relative_to(ROOT)}')
