from pathlib import Path
import urllib.request,json,hashlib,concurrent.futures,datetime
P=Path(__file__).parent;out=P/'public_metadata';out.mkdir(exist_ok=True)
requests=[]
for short,repo in [('author','RAISIN_cosmo'),('release','RAISIN_DataRelease')]:
 for label,tail in [('repository',''),('branches','/branches?per_page=100'),('tags','/tags?per_page=100'),('releases','/releases?per_page=100')]:
  requests.append((short+'-'+label,'https://api.github.com/repos/djones1040/'+repo+tail))
def get(x):
 label,url=x;record={'name':label,'url':url,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:
  req=urllib.request.Request(url,headers={'User-Agent':'supernova-provenance-audit','Accept':'application/vnd.github+json'})
  with urllib.request.urlopen(req,timeout=30) as f:
   b=f.read(12_000_001);assert len(b)<=12_000_000;record.update(status=f.status,headers=dict(f.headers),size=len(b))
  p=out/(label+'.json');p.write_bytes(b);record.update(path=str(p),sha256=hashlib.sha256(b).hexdigest());d=json.loads(b)
  if isinstance(d,list):record['summary']=[{'name':r.get('name',r.get('tag_name')),'sha':r.get('commit',{}).get('sha'),'assets':[{k:a.get(k) for k in ['name','size','browser_download_url']} for a in r.get('assets',[])]} for r in d]
  else:record['summary']={k:d.get(k) for k in ['default_branch','pushed_at','size','archived','html_url']}
 except Exception as e:record['error']=repr(e)
 return record
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:r=list(pool.map(get,requests))
(P/'public-metadata-acquisition.json').write_text(json.dumps(r,indent=2)+'\n')
for x in r:print({k:v for k,v in x.items() if k in ['name','status','size','summary','error']})
