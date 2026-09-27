"""Download pinned public host catalogues and joint posterior archive."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import urllib.request

ROOT = Path(__file__).resolve().parents[4]
WORK = ROOT/'.work/unified-cosmology/host-likelihood'
OUT = ROOT/'studies/unified_cosmology/results/host_likelihood'

def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--posterior-archive',action='store_true')
    args=parser.parse_args()
    WORK.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
    url='https://zenodo.org/api/records/16953206'
    meta=WORK/'frankenblast-record.json'
    if not meta.exists():meta.write_bytes(urllib.request.urlopen(url,timeout=60).read())
    release=json.loads(meta.read_text())
    records=[]
    for item in release['files']:
        name=item['key']
        if not name.endswith('.csv') and not(args.posterior_archive and name=='spectroscopic_SN_host_sbipp_output.zip'):continue
        path=WORK/name
        if not path.exists() or path.stat().st_size!=item['size']:
            temp=path.with_suffix(path.suffix+'.partial')
            with urllib.request.urlopen(item['links']['self'],timeout=120) as src,temp.open('wb') as dst:
                while block:=src.read(4*1024*1024):dst.write(block)
            assert temp.stat().st_size==item['size'],name
            temp.replace(path)
        algorithm,digest=item['checksum'].split(':')
        with path.open('rb') as stream:actual=hashlib.file_digest(stream,algorithm).hexdigest()
        assert actual==digest,(name,actual,digest)
        records.append(dict(path=str(path.relative_to(ROOT)),url=item['links']['self'],bytes=path.stat().st_size,
                            published_checksum=item['checksum'],sha256=sha(path)))
        print(name,path.stat().st_size,flush=True)
    result=dict(completed_utc=datetime.now(timezone.utc).isoformat(),record_url=url,record_sha256=sha(meta),
                code_sha256=sha(__file__),design_sha256=sha(Path(__file__).with_name('design.json')),files=records)
    (OUT/'acquisition.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
