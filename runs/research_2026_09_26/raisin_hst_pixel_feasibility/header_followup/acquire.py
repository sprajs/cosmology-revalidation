from pathlib import Path
import json, hashlib, urllib.request, urllib.parse, time
P=Path(__file__).resolve().parent
protocol=json.loads((P/'protocol.json').read_text())
result=[]
for uri in protocol['uris']:
    url='https://mast.stsci.edu/api/v0.1/Download/file?'+urllib.parse.urlencode({'uri':uri})
    rec={'uri':uri,'url':url}; start=time.monotonic()
    try:
        request=urllib.request.Request(url,headers={'Range':'bytes=0-287999','User-Agent':'SupernovaCalibrationAudit/1.0'})
        with urllib.request.urlopen(request,timeout=20) as r:
            rec.update(http_status=r.status, headers=dict(r.headers),resolved_url=r.url)
            blocks=[]
            for index in range(100):
                block=r.read(2880)
                if len(block)!=2880:raise ValueError('incomplete FITS header block')
                if index==0 and block[:8]!=b'SIMPLE  ':raise ValueError('response is not FITS')
                blocks.append(block)
                cards=[block[i:i+80] for i in range(0,2880,80)]
                if any(card[:8]==b'END     ' for card in cards):break
            else:raise ValueError('header cap exceeded')
            raw=b''.join(blocks); name=uri.rsplit('/',1)[-1]+'.primary-header'
            (P/name).write_bytes(raw)
            rec.update(success=True,header_bytes=len(raw),header_file=name,sha256=hashlib.sha256(raw).hexdigest())
    except Exception as e:rec.update(success=False,error=type(e).__name__+': '+str(e))
    rec['elapsed_seconds']=time.monotonic()-start;result.append(rec)
    (P/'result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps([{k:v for k,v in r.items() if k not in ['headers','url','resolved_url']} for r in result],indent=2))
