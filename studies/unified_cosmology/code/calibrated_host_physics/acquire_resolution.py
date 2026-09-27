#!/usr/bin/env python3
"""Recover the pinned author definition of the native DESI DIA response."""
import hashlib
import json
from pathlib import Path
import urllib.request
from datetime import datetime,timezone

ROOT=Path(__file__).resolve().parents[4]
WORK=ROOT/'.work/unified-cosmology/calibrated-host-physics'
RESULT=ROOT/'studies/unified_cosmology/results/calibrated_host_physics'
URL='https://raw.githubusercontent.com/desihub/desispec/661d87dee6d374ce343b754089d41580f468459a/py/desispec/resolution.py'
EXPECTED='a3562748ec82253121d3e08008a99e4bac7482be5d3e9481ed050ba675152284'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    WORK.mkdir(parents=True,exist_ok=True)
    path=WORK/'desispec-resolution.py'
    if not path.exists():
        with urllib.request.urlopen(URL) as response:path.write_bytes(response.read())
    assert sha(path)==EXPECTED
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),url=URL,
        commit='661d87dee6d374ce343b754089d41580f468459a',source_sha256=EXPECTED,
        code_sha256={str(Path(__file__).relative_to(ROOT)):sha(__file__)},
        output_sha256={str(path.relative_to(ROOT)):sha(path)},
        scope='Reference source retained outside Git; response implementation independently checked against direct FITS DIA column scattering.')
    (RESULT/'resolution-reference.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
