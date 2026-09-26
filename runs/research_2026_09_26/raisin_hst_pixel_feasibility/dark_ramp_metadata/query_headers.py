#!/usr/bin/env python3
"""Metadata-only FITS primary-header Range query for frozen CAOM dark candidates."""
import hashlib
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

BASE = Path(__file__).resolve().parent
PRODUCTS = json.loads((BASE / 'product-results.json').read_text())
OUT = BASE / 'header-results.json'
START = time.monotonic()
TOTAL = int(PRODUCTS['total_metadata_bytes'])
KEYS = ('ROOTNAME','DETECTOR','SUBARRAY','SUBTYPE','SAMP_SEQ','NSAMP','EXPSTART','EXPTIME','ATODGNA','ATODGNB','ATODGNC','ATODGND','READNSEA','READNSEB','READNSEC','READNSED','CAL_VER','DARKFILE','NLINFILE','ZSIGFILE','PFLTFILE','DFLTFILE','LFLTFILE','UNITCORR','FLATCORR','CRCORR','ZOFFCORR','BLEVCORR','DQICORR','FLUXCORR','SAMPZERO')


def header(uri):
    global TOTAL
    url = 'https://mast.stsci.edu/api/v0.1/Download/file?uri=' + urllib.parse.quote(uri, safe='')
    blocks = []
    for ib in range(24):
        if time.monotonic() - START > 300 or TOTAL + 2880 > 20000000:
            raise RuntimeError('cumulative resource cap')
        lo, hi = ib*2880, (ib+1)*2880-1
        req = urllib.request.Request(url, headers={'Range': f'bytes={lo}-{hi}', 'User-Agent':'Codex-MAST-metadata-audit/1.0'})
        with urllib.request.urlopen(req, timeout=20) as res:
            status = res.status
            cr = res.headers.get('Content-Range','')
            if status != 206 or not cr.startswith(f'bytes {lo}-{hi}/'):
                raise RuntimeError(f'unsafe range response status={status} content-range={cr}')
            body = res.read(2881)
            if len(body) != 2880:
                raise RuntimeError(f'unexpected range body size {len(body)}')
        TOTAL += len(body)
        blocks.append(body)
        cards = [body[i:i+80] for i in range(0,2880,80)]
        if any(c.startswith(b'END ') for c in cards):
            raw = b''.join(blocks)
            vals={}
            for off in range(0,len(raw),80):
                line=raw[off:off+80].decode('ascii',errors='replace')
                key=line[:8].strip()
                if key=='END':break
                if line[8:10]=='= ' and key:
                    field=line[10:].split(' / ')[0].strip()
                    if field.startswith("'"):
                        value=field.split("'",2)[1]
                    elif field in ('T','F'):value=field=='T'
                    else:
                        try:value=float(field.replace('D','E')) if any(x in field for x in '.EeDd') else int(field)
                        except ValueError:value=field
                    vals[key]=value
            return {'keys':{k:vals.get(k) for k in KEYS},'all_key_names':list(vals),'header_bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest(),'uri':uri,'url':url}
    raise RuntimeError('primary header exceeds 24 blocks')


def save(result):
    result['total_metadata_bytes']=TOTAL
    result['elapsed_seconds']=time.monotonic()-START
    OUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')


def main():
    result={'status':'partial','visits':{},'source_products_sha256':hashlib.sha256((BASE/'product-results.json').read_bytes()).hexdigest()}
    for visit,info in PRODUCTS['visits'].items():
        products={p['productFilename']:p for p in info['products']}
        rows=[]
        eligible=0
        for obs in info['candidate_roots']:
            root=obs['obs_id']
            p=products.get(root+'_raw.fits')
            row={'root':root,'obsid':obs['obsid'],'mjd':obs['t_min'],'exptime_caom':obs['t_exptime'],'target':obs['target_name'],'proposal_id':obs['proposal_id'],'raw_available':p is not None,'ima_available':root+'_ima.fits' in products}
            try:
                if p is None:raise RuntimeError('RAW product absent')
                row.update(header(p['dataURI']))
                k=row['keys']
                row['exact_eligible']=(k['DETECTOR']=='IR' and k['SUBARRAY'] is False and k['SAMP_SEQ']=='SPARS50' and isinstance(k['NSAMP'],int) and k['NSAMP']>=8)
                eligible+=bool(row['exact_eligible'])
            except Exception as exc:
                row['error']=str(exc)
            rows.append(row)
            result['visits'][visit]={'records':rows,'eligible_count':eligible}
            save(result)
            if eligible>=8:break
            if time.monotonic()-START>280 or TOTAL>19500000:break
    result['status']='complete' if all(len(v['records'])>0 for v in result['visits'].values()) else 'partial'
    save(result)

if __name__=='__main__':main()
