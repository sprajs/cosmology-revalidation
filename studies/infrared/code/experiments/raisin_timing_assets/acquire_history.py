from pathlib import Path
import urllib.request,json,hashlib,concurrent.futures,datetime
P=Path(__file__).parent;O=P/'public_metadata'
reqs=[('author-20211111-tree','https://api.github.com/repos/djones1040/RAISIN_cosmo/git/trees/aaa709ead7a7339d56a2a4604d78991329d9368f?recursive=1'),('zenodo-latest','https://zenodo.org/api/records/6349657/versions/latest')]
def get(x):
 name,url=x;r={'name':name,'url':url,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:
  with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'supernova-provenance-audit'}),timeout=45) as f:b=f.read(15_000_001);assert len(b)<=15_000_000;r.update(status=f.status,headers=dict(f.headers))
  p=O/(name+'.json');p.write_bytes(b);r.update(bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),path=str(p))
 except Exception as e:r['error']=repr(e)
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:r=list(pool.map(get,reqs))
(P/'history-acquisition.json').write_text(json.dumps(r,indent=2)+'\n')
print([{k:v for k,v in x.items() if k in ['name','bytes','status','error']} for x in r])
