"""Download only three listed small support/log products for template pilot dark."""
import hashlib,json,time,urllib.parse,urllib.request
from pathlib import Path
B=Path(__file__).resolve().parent
P=json.loads((B/'support-products.json').read_text())
assert P['total_size']<10000000
OUT=B/'support_files';OUT.mkdir(exist_ok=True)
start=time.monotonic();used=0;result={'status':'partial','records':[]}
def save():
    result['elapsed_seconds']=time.monotonic()-start;result['bytes_transferred']=used
    (B/'support-acquisition.json').write_text(json.dumps(result,indent=2)+'\n')
for p in P['records']:
    if time.monotonic()-start>90 or used+p['size']>10000000:break
    uri=p['dataURI'];url='https://mast.stsci.edu/api/v0.1/Download/file?uri='+urllib.parse.quote(uri,safe='')
    path=OUT/p['productFilename'];r={'filename':p['productFilename'],'uri':uri,'expected_size':p['size']}
    result['records'].append(r);save()
    try:
        req=urllib.request.Request(url,headers={'User-Agent':'Codex-RAISIN-quality-support/1'})
        with urllib.request.urlopen(req,timeout=25) as response, path.open('wb') as f:
            if response.status!=200:raise RuntimeError(f'HTTP {response.status}')
            r['http_status']=response.status
            while True:
                if time.monotonic()-start>90:raise RuntimeError('90s time cap')
                chunk=response.read(min(1<<20,10000000-used))
                if not chunk:break
                f.write(chunk);used+=len(chunk)
        if path.stat().st_size!=p['size']:raise RuntimeError('size mismatch')
        r['sha256']=hashlib.sha256(path.read_bytes()).hexdigest();r['status']='complete'
    except Exception as exc:
        r['status']='failed';r['error']=str(exc);save();break
    save()
else:result['status']='complete';save()
