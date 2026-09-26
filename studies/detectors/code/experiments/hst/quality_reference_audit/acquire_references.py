"""Fetch exactly two header-named official CRDS calibration tables."""
import hashlib,json,time,urllib.request
from pathlib import Path
B=Path(__file__).resolve().parent
P=json.loads((B/'reference-products.json').read_text());assert P['total_size']<5000000
OUT=B/'reference_files';OUT.mkdir(exist_ok=True)
start=time.monotonic();total=0;out={'status':'partial','records':[]}
for x in P['records']:
    if time.monotonic()-start>90 or total+x['size']>5000000:break
    p=OUT/x['filename'];entry={'filename':x['filename'],'url':x['url'],'expected_size':x['size']};out['records'].append(entry)
    try:
        req=urllib.request.Request(x['url'],headers={'User-Agent':'Codex-RAISIN-reference-audit/1'})
        with urllib.request.urlopen(req,timeout=20) as res,p.open('wb') as f:
            if res.status!=200:raise RuntimeError(f'HTTP {res.status}')
            while True:
                if time.monotonic()-start>90:raise RuntimeError('time cap')
                chunk=res.read(min(65536,5000000-total))
                if not chunk:break
                total+=len(chunk);f.write(chunk)
        if p.stat().st_size!=x['size']:raise RuntimeError('size mismatch')
        entry['status']='complete';entry['sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
    except Exception as exc:
        entry['status']='failed';entry['error']=str(exc);break
else:out['status']='complete'
out['bytes_transferred']=total;out['elapsed_seconds']=time.monotonic()-start
(B/'reference-acquisition.json').write_text(json.dumps(out,indent=2)+'\n')
