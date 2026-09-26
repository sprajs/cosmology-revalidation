"""Verify four exact MAST RAW product names and sizes before acquisition."""
import hashlib,json,urllib.parse,urllib.request
from pathlib import Path

B=Path(__file__).resolve().parent
M=B.parent/'dark_ramp_metadata'
ROOTS=('idbx43p7q','idp247tnq','icxoi1bcq','icxoi4hgq')
dark={}
for name in ('search_remaining_products.json','template_remaining_products.json'):
    p=M/'queries'/name
    d=json.loads(p.read_text())
    for x in d['data']:
        if x.get('productFilename') in {r+'_raw.fits' for r in ROOTS[:2]}:
            dark[x['productFilename']]=x
out={'purpose':'exact four RAW product URI/size verification before frozen acquisition','records':[],'source_product_json_sha256':{}}
for name in ('search_remaining_products.json','template_remaining_products.json'):
    p=M/'queries'/name;out['source_product_json_sha256'][name]=hashlib.sha256(p.read_bytes()).hexdigest()
for root in ROOTS:
    filename=root+'_raw.fits';uri='mast:HST/product/'+filename
    if root in ROOTS[:2]:
        record=dark[filename]
        assert record['dataURI']==uri and record['dataRights']=='PUBLIC'
    else:record=None
    url='https://mast.stsci.edu/api/v0.1/Download/file?uri='+urllib.parse.quote(uri,safe='')
    req=urllib.request.Request(url,method='HEAD',headers={'User-Agent':'Codex-RAISIN-RAW-pilot/1'})
    with urllib.request.urlopen(req,timeout=25) as res:
        assert res.status==200
        h={k:v for k,v in res.headers.items() if k.lower() in ('content-length','etag','last-modified','content-type','accept-ranges','content-disposition')}
    size=int(h['Content-Length'])
    if record is not None:assert size==record['size']
    out['records'].append({'root':root,'filename':filename,'uri':uri,'url':url,'size':size,'product_obsid':None if record is None else record['obsID'],'product_record_available':record is not None,'head_status':200,'head_headers':h})
out['total_size']=sum(x['size'] for x in out['records'])
(B/'verified-products.json').write_text(json.dumps(out,indent=2)+'\n')
