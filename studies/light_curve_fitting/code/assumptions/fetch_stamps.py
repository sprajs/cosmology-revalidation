from pathlib import Path
import urllib.request,json,datetime,hashlib
from astropy.io import fits
O=Path('phase2/assumptions/image_access');v=json.load(open(O/'sia_des_dr2_se_1246314.json'));chosen=[]
for ref in ['D00230179_','D00230180_','D00230181_']:
 for r in v:
  if ref in r['access_url']:
   url=r['access_url'].replace('SIZE=0.003,0.003','SIZE=0.02,0.02');name=url.split('siaRef=')[1].split('&')[0]+'.ext'+url.split('extn=')[1].split('&')[0]+'.cutout.fits';p=O/name
   entry={'url':url,'original_metadata':r,'name':name,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
   try:
    response=urllib.request.urlopen(url,timeout=50);p.write_bytes(response.read());entry.update(bytes=p.stat().st_size,sha256=hashlib.file_digest(p.open('rb'),'sha256').hexdigest())
    with fits.open(p) as h:
     entry['hdus']=[{'shape':None if x.data is None else list(x.data.shape),'header':dict(x.header)} for x in h];print(name,[(x.data.shape if x.data is not None else None) for x in h],flush=True)
   except Exception as e:entry['error']=str(e)
   chosen.append(entry);(O/'stamp_manifest.json').write_text(json.dumps(chosen,indent=2,default=str)+'\n')
