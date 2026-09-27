#!/usr/bin/env python3
"""Acquire the original spectroscopic DES3YR release as a separate physical route."""
import json,tarfile
from pathlib import Path
from common import ROOT,WORK,RESULTS,sha
from acquire import get

def main():
    folder=WORK/'des3yr';folder.mkdir(exist_ok=True)
    url='https://desdr-server.ncsa.illinois.edu/despublic/sn_files/y3/tar_files/PublicDataRelease.tar.gz'
    r=get(url,folder/'PublicDataRelease.tar.gz');members=[]
    with tarfile.open(r['path']) as t:
        for m in t:
            if not m.isfile():continue
            p=Path(m.name);assert not p.is_absolute() and '..' not in p.parts
            q=folder/p;q.parent.mkdir(parents=True,exist_ok=True)
            if not q.exists():
                with t.extractfile(m) as src,q.open('wb') as dst:
                    while b:=src.read(2**20):dst.write(b)
            members.append({'path':str(q.relative_to(ROOT)),'bytes':m.size,'sha256':sha(q)})
    nested=[]
    for archive in sorted(folder.glob('0*.tar.gz')):
        with tarfile.open(archive) as t:
            for m in t:
                if not m.isfile():continue
                p=Path(m.name);assert not p.is_absolute() and '..' not in p.parts
                if any(v.startswith('._') for v in p.parts):continue
                # Chains alreadypublishedare notneededforphysical-event closure.
                if 'COSMOMC_CHAINS' in p.parts:continue
                q=folder/p;q.parent.mkdir(parents=True,exist_ok=True)
                if not q.exists():
                    with t.extractfile(m) as src,q.open('wb') as dst:
                        while b:=src.read(2**20):dst.write(b)
                nested.append({'archive':archive.name,'path':str(q.relative_to(ROOT)),'bytes':m.size,'sha256':sha(q)})
    result={'code_sha256':sha(__file__),'author_archive':r,'members':members,'nested_members':nested,'scope':'OriginalDES3YR spectroscopicallytypedrelease; distinctphotometry,SALTcalibration,scattermodelandspectroscopictargetingmodel. Do notsubstituteforupdatedDovekieorapply3yrspecselectionto5yrSNe.'}
    (RESULTS/'des3yr-acquisition.json').write_text(json.dumps(result,indent=2)+'\n');print('members',len(members))
if __name__=='__main__':main()
