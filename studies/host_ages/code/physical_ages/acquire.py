#!/usr/bin/env python3
"""Acquire pinned stellar models and Galactic foreground maps, outside Git."""
import hashlib,json,subprocess,urllib.request
from datetime import datetime,timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).parent
WORK=ROOT/'.work/physical-ages'
OUT=ROOT/'studies/host_ages/results/physical_ages'
FSPS_COMMIT='05b5e550ddd7ffb1e71102b75cfbfdf057c211fd'
SFD_COMMIT='7a5fe7fadf086561ba4748756e59a4c51d0ec632'

def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def main():
    WORK.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
    folder=WORK/'fsps'
    if not folder.exists():
        subprocess.run(['git','clone','--depth','1','--branch','v4.0','https://github.com/cconroy20/fsps.git',str(folder)],check=True)
    commit=subprocess.check_output(['git','-C',str(folder),'rev-parse','HEAD'],text=True).strip()
    assert commit==FSPS_COMMIT,commit
    subprocess.run(['git','-C',str(folder),'diff','--exit-code','HEAD'],check=True,stdout=subprocess.DEVNULL)
    records=[]
    for name in ['SFD_dust_4096_ngp.fits','SFD_dust_4096_sgp.fits']:
        dest=WORK/'sfd'/name;url=f'https://raw.githubusercontent.com/kbarbary/sfddata/{SFD_COMMIT}/{name}'
        if not dest.exists():
            dest.parent.mkdir(parents=True,exist_ok=True)
            temp=dest.with_suffix('.partial')
            with urllib.request.urlopen(url,timeout=60) as response,temp.open('wb') as f:
                while chunk:=response.read(1048576):f.write(chunk)
            temp.replace(dest)
        records.append(dict(path=str(dest.relative_to(ROOT)),url=url,sha256=sha(dest),bytes=dest.stat().st_size))
    previous=OUT/'acquisition.json'
    if previous.exists():
        old=json.loads(previous.read_text());expected={x['path']:x['sha256'] for x in old['files']}
        assert all(expected.get(x['path'],x['sha256'])==x['sha256'] for x in records)
    tracked=subprocess.check_output(['git','-C',str(folder),'ls-files','-z']).decode().split('\0')
    identities={p:sha(folder/p) for p in tracked if p}
    # Full stellar-data byte identities stay in a compact local manifest; the
    # committed record pins that manifest and the public immutable Git tree.
    inventory=WORK/'stellar-file-hashes.json';inventory.write_text(json.dumps(identities,indent=2)+'\n')
    result=dict(acquired_utc=datetime.now(timezone.utc).isoformat(),fsps_repository='https://github.com/cconroy20/fsps',fsps_commit=commit,stellar_files=len(identities),stellar_inventory_path=str(inventory.relative_to(ROOT)),stellar_inventory_sha256=sha(inventory),files=records,code_sha256=sha(Path(__file__)),scope='Third-party data/code stay ignored; Git identity and complete tracked byte inventory retained. Stellar-model assumptions are not validated by an identity check.')
    previous.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='files'},indent=2))

if __name__=='__main__':main()
