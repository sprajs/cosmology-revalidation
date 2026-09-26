from pathlib import Path
import urllib.request,json,datetime,hashlib
from astropy.io import fits
O=Path('phase2/assumptions/image_access');v=json.load(open(O/'stamp_manifest.json'));out=[]
for r in v:
 if 'ext1.' not in r['name']:continue
 u=r['url'].replace('extn=1&','extn=3&');p=O/r['name'].replace('ext1.','ext3.');a={'url':u,'query_basis':'Probe next FITS extension via same public cutout API; not advertised by SIA, must inspect EXTNAME before interpreting.','utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:
  q=urllib.request.urlopen(u,timeout=45);p.write_bytes(q.read());a.update(status=q.status,path=str(p),bytes=p.stat().st_size,sha256=hashlib.file_digest(p.open('rb'),'sha256').hexdigest())
  with fits.open(p) as h:a['hdus']=[{'shape':None if x.data is None else list(x.data.shape),'header':dict(x.header)} for x in h]
 except Exception as e:a['error']=str(e)
 out.append(a);(O/'weight_access_manifest.json').write_text(json.dumps(out,indent=2,default=str)+'\n');print(a.get('error',a.get('hdus')),flush=True)
