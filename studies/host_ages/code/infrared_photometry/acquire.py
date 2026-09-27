#!/usr/bin/env python3
"""Acquire public IRSA catalogue rows around the fixed DES/SPIRE host cohort."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import hashlib, io, json, time
import pandas as pd

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).parent
WORK=ROOT/'.work/infrared-photometry'
OUT=ROOT/'studies/host_ages/results/infrared_photometry'
TABLES={'C3':['chandra_cat_f05','chandra_24_cat_f05','servscdfsi12'],
        'X3':['xmm_cat_s05','xmm_24_cat_s05','servsxmmi12']}
DOCS={
 'servs-survey-paper.pdf':'https://arxiv.org/pdf/1206.4060',
 'swire-columns.html':'https://irsa.ipac.caltech.edu/data/SPITZER/SWIRE/SWIRE_EN1_columns.html',
 'swire24-columns.html':'https://irsa.ipac.caltech.edu/data/SPITZER/SWIRE/SWIRE_24only_columns.html',
 'swire-delivery.pdf':'https://irsa.ipac.caltech.edu/data/SPITZER/SWIRE/docs/delivery_doc_r2_v2.pdf',
 'servs-delivery.pdf':'https://irsa.ipac.caltech.edu/data/SPITZER/SERVS/docs/SERVS_DR1_v1.4.pdf',
 'servs-overview.html':'https://irsa.ipac.caltech.edu/data/SPITZER/SERVS/overview.html',
}
PINNED={}
if (OUT/'acquisition.json').exists():
    previous=json.loads((OUT/'acquisition.json').read_text())
    for record in previous['catalogues']:
        for item in record['requests']:
            PINNED[item['path']]=(item['sha256'],item['url'])
    for item in previous['documentation']:
        PINNED[item['path']]=(item['sha256'],item['url'])

def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def sample():
    p=ROOT/'.work/host-transport'
    q=pd.read_csv(p/'des-deep-crosswalk.csv',dtype={'SNID':str,'deep_ID':'Int64'})
    q=q[q.primary_match & q.selected_Dovekie & q.deep_quality & q.host_dlr_lt4]
    ir=pd.read_csv(p/'des-spire-photometry.csv',dtype={'sn_id':str})
    ir=ir[ir.eightband_331 & ir.valid_map_measurement & (ir.common_sky_offsets>=64)]
    ids=ir.groupby('sn_id').band_um.nunique();q=q[q.SNID.isin(ids[ids==3].index)]
    assert len(q)==265 and q.deep_ID.nunique()==265
    return q.sort_values(['field','SNID'])

def fetch(name,url):
    path=WORK/name
    key=str(path.relative_to(ROOT));pin=PINNED.get(key)
    if pin is not None:
        assert url==pin[1],f'Changed pinned query or URL: {key}'
    # Never silently overwrite a pinned response. Existing errors are retried.
    if path.exists() and b'QUERY_STATUS" value="ERROR' not in path.read_bytes():
        if pin is not None:assert sha(path)==pin[0],f'Changed pinned response: {key}'
        return path
    for attempt in range(4):
        try:
            raw=urlopen(url,timeout=60).read()
            if b'QUERY_STATUS" value="ERROR' in raw:raise RuntimeError(raw.decode()[:1000])
            if pin is not None:assert hashlib.sha256(raw).hexdigest()==pin[0],f'Upstream response changed: {key}'
            path.write_bytes(raw);return path
        except Exception:
            if attempt==3:raise
            time.sleep(1+attempt)

def tap(name,query):
    url='https://irsa.ipac.caltech.edu/TAP/sync?'+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',MAXREC=200000,QUERY=query))
    path=fetch(name,url)
    return {'path':str(path.relative_to(ROOT)),'sha256':sha(path),'url':url,'query':query,'bytes':path.stat().st_size}

def acquire_table(job):
    field,table,hosts=job;records=[];frames=[]
    records.append(tap(table+'-columns.csv',f"SELECT column_name,datatype,unit,description FROM TAP_SCHEMA.columns WHERE table_name='{table}'"))
    for i in range(0,len(hosts)):
        cones=[f"CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{r.deep_RA:.12f},{r.deep_DEC:.12f},{75/3600:.12f}))=1" for r in hosts.iloc[i:i+1].itertuples()]
        query=f'SELECT * FROM {table} WHERE '+ ' OR '.join(cones)
        rec=tap(f'{table}-{i:03d}.csv',query);records.append(rec)
        raw=pd.read_csv(ROOT/rec['path'],dtype=str)
        assert 'cntr' in raw and 'ra' in raw
        frames.append(raw)
    frame=pd.concat(frames,ignore_index=True)
    # Same ID returned through multiple overlapping cones must have same values.
    assert frame.groupby('cntr',dropna=False).nunique(dropna=False).max().max()<=1
    frame=frame.drop_duplicates('cntr').sort_values('cntr')
    frame.to_csv(WORK/f'{table}.csv',index=False)
    return {'field':field,'catalogue':table,'rows':len(frame),'path':str((WORK/f'{table}.csv').relative_to(ROOT)),
            'sha256':sha(WORK/f'{table}.csv'),'requests':records}

def main():
    WORK.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
    hosts=sample();hosts[['SNID','deep_ID','field','deep_RA','deep_DEC','redshift']].to_csv(WORK/'hosts.csv',index=False)
    inputs={str((ROOT/'.work/host-transport'/n).relative_to(ROOT)):sha(ROOT/'.work/host-transport'/n) for n in ['des-deep-crosswalk.csv','des-spire-photometry.csv']}
    jobs=[(field,t,hosts[hosts.field==field]) for field,ts in TABLES.items() for t in ts]
    with ThreadPoolExecutor(6) as pool:catalogues=list(pool.map(acquire_table,jobs))
    docs=[]
    for name,url in DOCS.items():
        p=fetch(name,url);docs.append(dict(path=str(p.relative_to(ROOT)),url=url,sha256=sha(p)))
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),hosts=265,host_fields=hosts.field.value_counts().to_dict(),
        acquisition_radius_arcsec=75,inputs=inputs,catalogues=catalogues,documentation=docs,
        code_sha256=sha(Path(__file__)),design_sha256=sha(HERE/'design.json'),hosts_sha256=sha(WORK/'hosts.csv'))
    result['acquisition_amendment_sha256']=sha(HERE/'acquisition-amendment.json')
    (OUT/'acquisition.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({c['catalogue']:c['rows'] for c in catalogues}))

if __name__=='__main__':main()
