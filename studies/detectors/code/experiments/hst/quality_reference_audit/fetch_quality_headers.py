"""Frozen metadata-only primary-header Range audit of 14 dark candidates."""
import csv,hashlib,json,time,urllib.parse,urllib.request
from pathlib import Path

B=Path(__file__).resolve().parent
ROOT=B.parent
C=list(csv.DictReader(open(ROOT/'dark_ramp_metadata/candidate-cohort.csv')))
assert len(C)==14 and len({r['root'] for r in C})==14
LOCAL=ROOT/'dark_ramp_raw_pilot/files'
HEADERS=B/'primary_headers';HEADERS.mkdir(exist_ok=True)
START=time.monotonic();used=0;out={'status':'partial','records':[]}

def parse_cards(raw):
    vals={}
    for i in range(0,len(raw),80):
        c=raw[i:i+80].decode('ascii','replace');key=c[:8].strip()
        if key=='END':break
        if c[8:10]=='= ' and key:
            field=c[10:].split(' / ')[0].strip()
            if field.startswith("'"):value=field.split("'",2)[1].strip()
            elif field in ('T','F'):value=field=='T'
            else:
                try:value=float(field.replace('D','E')) if any(x in field for x in '.EeDd') else int(field)
                except ValueError:value=field
            vals[key]=value
    return vals

def read_local(path):
    with path.open('rb') as f:
        blocks=[]
        while len(blocks)<24:
            b=f.read(2880)
            if len(b)!=2880:raise RuntimeError('short local header block')
            blocks.append(b)
            if any(b[i:i+80].startswith(b'END ') for i in range(0,2880,80)):
                return b''.join(blocks)
    raise RuntimeError('local primary header >24 blocks')

def read_range(uri):
    global used
    url='https://mast.stsci.edu/api/v0.1/Download/file?uri='+urllib.parse.quote(uri,safe='')
    blocks=[]
    for j in range(24):
        if time.monotonic()-START>90 or used+2880>2000000:raise RuntimeError('phase resource cap')
        lo=j*2880;hi=lo+2879
        req=urllib.request.Request(url,headers={'Range':f'bytes={lo}-{hi}','User-Agent':'Codex-RAISIN-quality-header/1'})
        with urllib.request.urlopen(req,timeout=20) as r:
            cr=r.headers.get('Content-Range','')
            if r.status!=206 or not cr.startswith(f'bytes {lo}-{hi}/'):raise RuntimeError(f'unsafe HTTP range {r.status} {cr}')
            block=r.read(2881)
        if len(block)!=2880:raise RuntimeError('unexpected range length')
        used+=len(block);blocks.append(block)
        if any(block[i:i+80].startswith(b'END ') for i in range(0,2880,80)):
            return b''.join(blocks)
    raise RuntimeError('remote primary header >24 blocks')

def save():
    out['elapsed_seconds']=time.monotonic()-START;out['network_header_bytes']=used
    (B/'quality-header-results.json').write_text(json.dumps(out,indent=2)+'\n')

for row in C:
    root=row['root'];filename=root+'_raw.fits';uri='mast:HST/product/'+filename
    result={'root':root,'visit':row['visit'],'uri':uri,'source':'local pilot' if (LOCAL/filename).exists() else 'MAST HTTP206 primary header'}
    out['records'].append(result);save()
    try:
        raw=read_local(LOCAL/filename) if (LOCAL/filename).exists() else read_range(uri)
        (HEADERS/(root+'.hdr')).write_bytes(raw)
        result.update({'status':'complete','header_bytes':len(raw),'header_sha256':hashlib.sha256(raw).hexdigest(),'keys':parse_cards(raw)})
    except Exception as exc:
        result.update({'status':'failed','error':str(exc)})
        save();break
    save()
else:out['status']='complete';save()
