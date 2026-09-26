#!/usr/bin/env python3
"""Fetch only frozen small author files and verify exact Git blobs."""
import hashlib
import json
from pathlib import Path
import urllib.parse
import urllib.request

OUT = Path(__file__).resolve().parent
plan = json.loads((OUT/'source-plan.json').read_text())
assert len(plan['files']) == 11 and sum(x['size'] for x in plan['files']) < plan['max_bytes']
records=[]
for item in plan['files']:
    path=item['path']
    target=OUT/'author'/path
    assert not target.exists()
    url='https://raw.githubusercontent.com/djones1040/RAISIN_cosmo/'+plan['commit']+'/'+urllib.parse.quote(path)
    with urllib.request.urlopen(url, timeout=30) as response:
        data=response.read(item['size']+1)
    assert len(data)==item['size'],(path,len(data),item['size'])
    blob=hashlib.sha1(f'blob {len(data)}\0'.encode()+data).hexdigest()
    assert blob==item['git_blob'],(path,blob,item['git_blob'])
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_bytes(data)
    records.append({'path':path,'git_blob':blob,'size':len(data),'sha256':hashlib.sha256(data).hexdigest(),'url':url})
(OUT/'acquisition.json').write_text(json.dumps({'commit':plan['commit'],'files':records},indent=2)+'\n')
print('verified',len(records),'bytes',sum(x['size'] for x in records))
