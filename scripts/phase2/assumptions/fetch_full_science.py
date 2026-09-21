from pathlib import Path
import urllib.request,json,hashlib,datetime
from astropy.io import fits
O=Path('phase2/assumptions/image_access');out=[]
for r in json.load(open(O/'stamp_manifest.json')):
 if 'ext1.' not in r['name']:continue
 u=r['url'].replace('SIZE=0.02,0.02','SIZE=0.4,0.4');p=O/r['name'].replace('.cutout.','.fullccd.');q=urllib.request.urlopen(u,timeout=50);p.write_bytes(q.read())
 h=fits.getheader(p);out.append({'url':u,'path':str(p),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest(),'NAXIS1':h['NAXIS1'],'NAXIS2':h['NAXIS2'],'CRPIX1':h['CRPIX1'],'CRPIX2':h['CRPIX2'],'purpose':'Recover uncut CCD WCS to align full weight plane with the small science/mask cutout.'});(O/'full_science_manifest.json').write_text(json.dumps(out,indent=2)+'\n');print(out[-1],flush=True)
