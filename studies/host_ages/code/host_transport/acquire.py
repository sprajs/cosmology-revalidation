#!/usr/bin/env python3
"""Retrieve byte-pinned original host catalogues; never vendor downloaded data."""
from pathlib import Path
import argparse,hashlib,json,urllib.request
ROOT=Path(__file__).resolve().parents[4]
REGISTRY=ROOT/'studies/host_ages/results/host_transport/inputs.json'
def digest(p):
    with p.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify-only',action='store_true');args=parser.parse_args()
    rows=json.loads(REGISTRY.read_text())['files']
    for row in rows:
        p=ROOT/row['path']
        if not p.exists():
            if args.verify_only:raise FileNotFoundError(p)
            p.parent.mkdir(parents=True,exist_ok=True)
            temporary=p.with_suffix(p.suffix+'.part')
            with urllib.request.urlopen(row['url'],timeout=60) as source,temporary.open('wb') as dest:
                while block:=source.read(4*1024**2):dest.write(block)
            temporary.rename(p)
        assert digest(p)==row['sha256'],p
    print(f'Verified {len(rows)} pinned original sources.')
if __name__=='__main__':main()
