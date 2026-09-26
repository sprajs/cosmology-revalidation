"""Metadata-only mode-diagnostic Range headers for nearest 352.939514 s darks."""
import hashlib,json,time
from pathlib import Path
from query_caom import request
import query_headers as qh

BASE=Path(__file__).resolve().parent
P=json.loads((BASE/'protocol.json').read_text())
C=json.loads((BASE/'caom-results.json').read_text())
H=json.loads((BASE/'header-results.json').read_text())
assert all(x['eligible_count']==0 for x in H['visits'].values())
started=time.monotonic();total_time=C['elapsed_seconds']+json.loads((BASE/'product-results.json').read_text())['product_seconds']+H['elapsed_seconds']+json.loads((BASE/'expanded-caom-results.json').read_text())['elapsed_seconds']
qh.START=started-(total_time); qh.TOTAL=json.loads((BASE/'expanded-caom-results.json').read_text())['cumulative_metadata_bytes']
result={'status':'partial','visits':{},'prior_elapsed_seconds':total_time,'prior_metadata_bytes':qh.TOTAL}
for visit,mid in P['science_visit_midpoints_mjd'].items():
    candidates=sorted((r for r in C['visits'][visit][0]['rows'] if abs(r['t_exptime']-352.939514)<.001),key=lambda r:(abs(r['t_min']-mid),r['obs_id']))[:8]
    payload={'service':'Mast.Caom.Products','format':'json','params':{'obsid':','.join(str(r['obsid']) for r in candidates)}}
    prod,rec=request(payload,visit+'_352s_products',P['limits'],started,qh.TOTAL)
    qh.TOTAL+=rec['bytes']; products={x['productFilename']:x for x in prod['data']}
    records=[]
    for obs in candidates:
        root=obs['obs_id'];p=products.get(root+'_raw.fits')
        row={'root':root,'obsid':obs['obsid'],'mjd':obs['t_min'],'exptime_caom':obs['t_exptime'],'proposal_id':obs['proposal_id'],'raw_available':p is not None,'ima_available':root+'_ima.fits' in products}
        try:
            if p is None:raise RuntimeError('RAW product absent')
            row.update(qh.header(p['dataURI']))
            k=row['keys'];row['exact_eligible']=(str(k['DETECTOR']).strip()=='IR' and k['SUBARRAY'] is False and str(k['SAMP_SEQ']).strip()=='SPARS50' and isinstance(k['NSAMP'],int) and k['NSAMP']>=8)
        except Exception as exc:row['error']=str(exc)
        records.append(row)
        result['visits'][visit]={'records':records,'query':rec}
        result['cumulative_metadata_bytes']=qh.TOTAL;result['cumulative_seconds']=total_time+time.monotonic()-started
        (BASE/'targeted-header-results.json').write_text(json.dumps(result,indent=2)+'\n')
        if result['cumulative_seconds']>290:break
result['status']='complete' if all(len(x['records'])==8 for x in result['visits'].values()) else 'partial'
(BASE/'targeted-header-results.json').write_text(json.dumps(result,indent=2)+'\n')
