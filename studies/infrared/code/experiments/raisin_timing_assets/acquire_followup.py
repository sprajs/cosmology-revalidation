from pathlib import Path
import urllib.request,urllib.parse,hashlib,json,concurrent.futures,datetime
P=Path(__file__).parent;O=P/'public_metadata';O.mkdir(exist_ok=True)
base='https://api.github.com/repos/djones1040/'
reqs=[('release-v1p0-tree',base+'RAISIN_DataRelease/git/trees/318fbc4ff20b92dd94cd8fc0cef544d139e9afd3?recursive=1'),('author-commits',base+'RAISIN_cosmo/commits?per_page=100'),('zenodo-6349657','https://zenodo.org/api/records/6349657')]
for name in ['README.md','.gitignore']:
 reqs.append(('author-'+name.replace('.','_'),'https://raw.githubusercontent.com/djones1040/RAISIN_cosmo/b888214a5cbae38ac0bf886488ce7734f5b77e87/'+name))
for typ in ['FITRES','LCPLOT']:
 path='output/fit_nir/DES_RAISIN_NIR_SIM/DES_RAISIN_SIM/FITOPT000.'+typ+'.gz'
 reqs.append(('author-history-'+typ,base+'RAISIN_cosmo/commits?per_page=100&path='+urllib.parse.quote(path)))
def fetch(x):
 name,url=x;r={'name':name,'url':url,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:
  with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'supernova-provenance-audit'}),timeout=30) as f:b=f.read(12_000_001);assert len(b)<=12_000_000;r.update(status=f.status,headers=dict(f.headers))
  p=O/(name+'.txt' if name.startswith('author-README') or name.startswith('author-_gitignore') else name+'.json');p.write_bytes(b);r.update(path=str(p),bytes=len(b),sha256=hashlib.sha256(b).hexdigest())
 except Exception as e:r['error']=repr(e)
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:out=list(pool.map(fetch,reqs))
(P/'followup-acquisition.json').write_text(json.dumps(out,indent=2)+'\n')
for r in out:print({k:v for k,v in r.items() if k in ['name','status','bytes','error']})
