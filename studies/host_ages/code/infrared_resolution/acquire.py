#!/usr/bin/env python3
"""Recover the public empirical SPIRE beams without bundling third-party data."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import urllib.request

ROOT=Path(__file__).resolve().parents[4]
WORK=ROOT/'.work/infrared-resolution'
OUT=ROOT/'studies/host_ages/results/infrared_resolution'
BASE='https://archives.esac.esa.int/hsa/legacy/ADP/PSF/SPIRE/SPIRE-P/'

def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def fetch(name):
    dest=WORK/name
    if not dest.exists():
        temp=dest.with_suffix('.partial')
        with urllib.request.urlopen(BASE+name,timeout=60) as response,temp.open('wb') as stream:
            while block:=response.read(1024*1024):stream.write(block)
        temp.replace(dest)
    return dict(path=str(dest.relative_to(ROOT)),url=BASE+name,bytes=dest.stat().st_size,sha256=sha(dest))

def main():
    WORK.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
    names=['README.html']+[f'0x5000241aL_{b}_bgmod10_1arcsec.fits.gz' for b in ['PSW','PMW','PLW']]
    with ThreadPoolExecutor(3) as pool:records=list(pool.map(fetch,names))
    target=OUT/'acquisition.json'
    if target.exists():
        old={x['path']:x['sha256'] for x in json.loads(target.read_text())['files']}
        assert all(old[r['path']]==r['sha256'] for r in records)
    record=dict(acquired_utc=datetime.now(timezone.utc).isoformat(),files=records,code_sha256=sha(Path(__file__)),scope='Official empirical Neptune beams, not exact exposure-weighted coadd PSFs or source-specific colour corrections.')
    target.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))

if __name__=='__main__':main()
