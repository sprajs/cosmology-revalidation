"""Recover pinned official SH0ES compact products, including Git LFS covariance."""
import json
from urllib.request import urlopen
from common import ROOT,HERE,WORK,OUT,sha,relative,write

def main():
 pins=json.loads((HERE/'sources.json').read_text());WORK.mkdir(parents=True,exist_ok=True);rows=[]
 for item in pins['files']:
  path=WORK/item['name']
  if not path.exists():
   with urlopen(item['url'],timeout=180)as response:body=response.read()
   assert len(body)==item['bytes']
   temporary=path.with_name(path.name+'.part');temporary.write_bytes(body)
   assert sha(temporary)==item['sha256'];temporary.replace(path)
  assert path.stat().st_size==item['bytes']and sha(path)==item['sha256'],item['name']
  rows.append(dict(item,path=relative(path)))
 pointer=(WORK/'covariance-lfs-pointer.txt').read_text()
 assert 'oid sha256:'+next(x['sha256']for x in rows if x['name'].startswith('allc_')) in pointer
 result={'status':'passed_pinned_acquisition','official_commit':pins['official_commit'],'files':rows,
  'source_sha256':{relative(p):sha(p)for p in[HERE/'acquire.py',HERE/'common.py',HERE/'sources.json',HERE/'design.json']},
  'scope':'Released compressed likelihood and author reference code; downloaded code is inspected, not installed or executed. No active cosmology inputs changed.'}
 write(OUT/'acquisition.json',result);print(json.dumps({'status':result['status'],'files':len(rows),'bytes':sum(r['bytes']for r in rows)}))
if __name__=='__main__':main()
