"""Frozen four-file public RAW acquisition; no FITS arrays are interpreted."""
import hashlib,json,time,urllib.request
from pathlib import Path

B=Path(__file__).resolve().parent
P=json.loads((B/'verified-products.json').read_text())
DEST=B/'files';DEST.mkdir(exist_ok=True)
START=time.monotonic();CAP_SECONDS=200;CAP_BYTES=140_000_000
used=0;out={'status':'partial','records':[],'verified_products_sha256':hashlib.sha256((B/'verified-products.json').read_bytes()).hexdigest()}
def save():
    out['elapsed_seconds']=time.monotonic()-START
    out['transferred_bytes']=used
    (B/'acquisition-results.json').write_text(json.dumps(out,indent=2)+'\n')
for src in P['records']:
    path=DEST/src['filename'];record={'filename':src['filename'],'uri':src['uri'],'expected_size':src['size'],'attempts':[]}
    out['records'].append(record);save()
    for attempt in range(2):
        if time.monotonic()-START>=CAP_SECONDS or used>=CAP_BYTES:
            record['status']='resource_stop';save();raise SystemExit(0)
        old=path.stat().st_size if path.exists() else 0
        headers={'User-Agent':'Codex-RAISIN-RAW-pilot/1'}
        if old:
            headers['Range']=f'bytes={old}-'
            headers['If-Range']=src['head_headers'].get('ETag','')
        entry={'attempt':attempt+1,'resume_offset':old,'started_seconds':time.monotonic()-START}
        record['attempts'].append(entry);save()
        try:
            req=urllib.request.Request(src['url'],headers=headers)
            with urllib.request.urlopen(req,timeout=30) as res:
                expected_status=206 if old else 200
                if res.status!=expected_status:raise RuntimeError(f'HTTP {res.status}, expected {expected_status}')
                if old and not res.headers.get('Content-Range','').startswith(f'bytes {old}-'):
                    raise RuntimeError('resume Content-Range mismatch')
                entry['http_status']=res.status
                with path.open('ab' if old else 'wb') as f:
                    while True:
                        if time.monotonic()-START>=CAP_SECONDS or used>=CAP_BYTES:
                            raise RuntimeError('cumulative resource cap')
                        chunk=res.read(min(1<<20,CAP_BYTES-used))
                        if not chunk:break
                        f.write(chunk);used+=len(chunk)
                        if used%(8<<20)<len(chunk):save()
            if path.stat().st_size!=src['size']:
                raise RuntimeError(f'file size {path.stat().st_size}, expected {src["size"]}')
            record['status']='downloaded'
            record['bytes']=path.stat().st_size
            record['sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
            entry['status']='complete';save();break
        except Exception as exc:
            entry['status']='failed';entry['error']=str(exc);entry['partial_size']=path.stat().st_size if path.exists() else 0
            save()
            if 'resource cap' in str(exc) or attempt==1:
                record['status']='failed';save();raise SystemExit(0)
out['status']='complete';save()
