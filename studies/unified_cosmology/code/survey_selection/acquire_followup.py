#!/usr/bin/env python3
"""Acquire rejected-candidate and failed-redshift records, not cosmology flux substitutes."""
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from common import WORK,RESULTS,sha
from acquire import get


def main():
    target=WORK/'followup';target.mkdir(exist_ok=True)
    metadata=target/'diffimg-zenodo.json';get('https://zenodo.org/api/records/14332953',metadata)
    z=json.loads(metadata.read_text())
    def one(f):
        r=get(f['links']['self'],target/f['key'])
        alg,digest=f['checksum'].split(':');assert sha(r['path'],alg)==digest
        r['release_checksum']=f['checksum'];return r
    with ThreadPoolExecutor(max_workers=3) as pool:records=list(pool.map(one,z['files']))
    for name in ['ReadMe','ozdesdr2.dat.gz']:
        records.append(get('https://cdsarc.cds.unistra.fr/ftp/J/MNRAS/496/19/'+name,target/name))
    result={'code_sha256':sha(__file__),'assets':records,
            'DIFFIMG_scope':'31,636 detected-candidate denominator and crosscheck only. Authors explicitly prohibit using DIFFIMG flux for cosmology; use SMP for inference.',
            'OzDES_scope':'Observed follow-up targets including insecure/no-redshift records; not a complete census of unobserved potential targets or undetectedIa.'}
    (RESULTS/'followup-acquisition.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'files':len(records),'bytes':sum(x['bytes'] for x in records)}))
if __name__=='__main__':main()
