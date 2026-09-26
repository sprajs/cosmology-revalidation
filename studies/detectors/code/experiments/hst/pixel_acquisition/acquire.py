"""Bounded acquisition only. No photometry or noise summaries."""
from pathlib import Path
import concurrent.futures, hashlib, json, time, urllib.request, urllib.parse
P=Path(__file__).resolve().parent
protocol=json.loads((P/'protocol.json').read_text())
def acquire(uri):
    start=time.monotonic(); name=uri.rsplit('/',1)[-1]
    record={'uri':uri,'filename':name}
    target=P/name; partial=P/(name+'.partial')
    try:
        req=urllib.request.Request('https://mast.stsci.edu/api/v0.1/Download/file?'+urllib.parse.urlencode({'uri':uri}),headers={'User-Agent':'SupernovaCalibrationAudit/1.0'})
        with urllib.request.urlopen(req,timeout=30) as r,partial.open('xb') as f:
            record['http_status']=r.status; record['headers']=dict(r.headers)
            if r.status!=200:raise ValueError('Full product requires HTTP200')
            total=0; digest=hashlib.sha256()
            while True:
                b=r.read(min(1048576,32000001-total))
                if not b:break
                total+=len(b)
                if total>32000000:raise ValueError('32MB cap')
                digest.update(b); f.write(b)
            if 'Content-Length' in r.headers and total!=int(r.headers['Content-Length']):raise ValueError('size mismatch')
            record.update(bytes=total,sha256=digest.hexdigest())
        partial.rename(target);record['success']=True
    except Exception as e:record.update(success=False,error=type(e).__name__+': '+str(e))
    record['elapsed_seconds']=time.monotonic()-start
    (P/(name+'.acquisition.json')).write_text(json.dumps(record,indent=2)+'\n')
    return record
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    rows=list(pool.map(acquire,protocol['uris']))
(P/'result.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps([{k:v for k,v in r.items() if k!='headers'} for r in rows],indent=2))
