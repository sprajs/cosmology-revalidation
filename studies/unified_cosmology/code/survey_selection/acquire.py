#!/usr/bin/env python3
"""Freeze official current release, public likelihood inputs and selection assets."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime
import hashlib
import json
from pathlib import Path
import tarfile
import urllib.request
from common import ROOT,HERE,WORK,RESULTS,sha
REPO='des-science/DES-SN5YR'
ZENODO=19503606


def get(url,path):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if not path.exists():
        temporary=path.with_suffix(path.suffix+'.partial')
        req=urllib.request.Request(url,headers={'User-Agent':'cosmology-revalidation-research'})
        with urllib.request.urlopen(req,timeout=90) as r,temporary.open('wb') as f:
            while b:=r.read(2**20):f.write(b)
        temporary.replace(path)
    return {'url':url,'path':str(path),'bytes':path.stat().st_size,'sha256':sha(path)}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--bundle',action='store_true');a=parser.parse_args()
    commit_path=WORK/'release-commit.json'
    get(f'https://api.github.com/repos/{REPO}/commits/main',commit_path)
    commit=json.loads(commit_path.read_text())['sha']
    treepath=WORK/'release-tree.json'
    get(f'https://api.github.com/repos/{REPO}/git/trees/{commit}?recursive=1',treepath)
    tree=json.loads(treepath.read_text());assert not tree.get('truncated')
    selected=[]
    for item in tree['tree']:
        p=item['path']
        if item['type']!='blob':continue
        keep=(p=='README.md' or p.startswith(('0_DATA/','3_CLASSIFICATION/','7_PIPPIN_FILES/'))
              or p in ['4_DISTANCES_COVMAT/README.md','4_DISTANCES_COVMAT/DES-Dovekie_HD.csv',
                       '4_DISTANCES_COVMAT/DES-Dovekie_Metadata.csv','4_DISTANCES_COVMAT/STAT+SYS.npz',
                       '4_DISTANCES_COVMAT/STATONLY.npz','4_DISTANCES_COVMAT/DES-Dovekie-SN_Likelihood.py']
              or p.startswith('4_DISTANCES_COVMAT/SingleSYS_CovMatrix/'))
        if keep:selected.append(item)
    def download(item):
        p=item['path'];record=get(f'https://raw.githubusercontent.com/{REPO}/{commit}/{p}',WORK/'release'/p)
        h=hashlib.sha1();h.update(f"blob {record['bytes']}\0".encode())
        with Path(record['path']).open('rb') as f:
            for b in iter(lambda:f.read(2**20),b''):h.update(b)
        assert h.hexdigest()==item['sha'],p
        record.update(release_path=p,git_blob_sha1=item['sha'],is_lfs_pointer=Path(record['path']).read_bytes()[:40].startswith(b'version https://git-lfs'))
        return record
    with ThreadPoolExecutor(max_workers=4) as pool:records=list(pool.map(download,selected))
    metadata=WORK/'sndataroot-v13.json';get(f'https://zenodo.org/api/records/{ZENODO}',metadata)
    z=json.loads(metadata.read_text());bundle=None
    if a.bundle:
        file=z['files'][0]
        archive=WORK/file['key']
        bundle=get(file['links']['self'],archive)
        alg,digest=file['checksum'].split(':');assert sha(archive,alg)==digest
        bundle['release_checksum']=file['checksum']
        extract=WORK/'SNDATA_ROOT_2026-04-10';extract.mkdir(exist_ok=True)
        members=[]
        with tarfile.open(archive,'r:gz') as tar:
            for member in tar:
                if not member.isfile():continue
                p=Path(member.name)
                if p.is_absolute() or '..' in p.parts:raise ValueError(p)
                # Only information-bearing survey/calibration/selection assets.
                if not any(s in member.name for s in ['sample_input_files/DES-SN5YR/','models/searcheff/','lcmerge/','kcor/','SIMLIB/','simlib/','HOSTLIB/','hostlib/']):continue
                target=extract/p;target.parent.mkdir(parents=True,exist_ok=True)
                if not target.exists():
                    with tar.extractfile(member) as src,target.open('wb') as dst:
                        while b:=src.read(2**20):dst.write(b)
                members.append({'member':member.name,'bytes':target.stat().st_size,'sha256':sha(target)})
        bundle['extracted_assets']=members
    result={'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'repository':REPO,'commit':commit,'tree_sha256':sha(treepath),
            'code_sha256':sha(__file__),'git_assets':records,
            'zenodo_record':ZENODO,'zenodo_metadata_sha256':sha(metadata),'new_bundle':bundle,
            'scope':'Official released photons, identifiers, corrected distances, covariance and documented inputs; no assertion that exact production realizations or full selection likelihood are reproduced.'}
    (RESULTS/'acquisition.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'commit':commit,'git_files':len(records),'lfs_pointers':sum(r['is_lfs_pointer'] for r in records),'bundle':bundle and bundle['bytes'],'extracted':bundle and len(bundle['extracted_assets'])}))
if __name__=='__main__':main()
