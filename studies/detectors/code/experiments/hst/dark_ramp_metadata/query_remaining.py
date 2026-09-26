"""Additive <=120 s, <=10 MB header-only continuation over remaining ±30 d roots."""
import hashlib,json,time
from pathlib import Path
from query_caom import request
import query_headers as qh

B=Path(__file__).resolve().parent
P=json.loads((B/'protocol.json').read_text())
C=json.loads((B/'caom-results.json').read_text())
H=json.loads((B/'header-results.json').read_text())
T=json.loads((B/'targeted-header-results.json').read_text())
start=time.monotonic();qh.START=start-180;qh.TOTAL=0
prior={v:{r['root'] for r in H['visits'][v]['records']}|{r['root'] for r in T['visits'][v]['records']} for v in ('search','template')}
by_visit={}
for v in ('search','template'):
    mid=P['science_visit_midpoints_mjd'][v]
    by_visit[v]=[r for r in sorted((x for x in C['visits'][v][0]['rows'] if x['t_exptime']>=300),key=lambda x:(abs(x['t_min']-mid),x['obs_id'])) if r['obs_id'] not in prior[v]]
out={'status':'partial','visits':{'search':{'records':[]},'template':{'records':[]}},'prior_header_sha256':hashlib.sha256((B/'header-results.json').read_bytes()).hexdigest(),'prior_targeted_sha256':hashlib.sha256((B/'targeted-header-results.json').read_bytes()).hexdigest()}
def save():
    out['elapsed_seconds']=time.monotonic()-start;out['metadata_bytes']=qh.TOTAL
    (B/'remaining-header-results.json').write_text(json.dumps(out,indent=2)+'\n')

# Product metadata for each remaining visit, no guessed filenames.
product_maps={}
for v,rows in by_visit.items():
    payload={'service':'Mast.Caom.Products','format':'json','params':{'obsid':','.join(str(r['obsid']) for r in rows)}}
    answer,rec=request(payload,v+'_remaining_products',{'wall_seconds':120,'metadata_total_bytes':10000000},start,qh.TOTAL)
    qh.TOTAL+=rec['bytes'];product_maps[v]={p['productFilename']:p for p in answer['data']}
    out['visits'][v]['product_query']=rec;save()

for rank in range(max(map(len,by_visit.values()))):
    for v in ('search','template'):
        if rank>=len(by_visit[v]):continue
        if time.monotonic()-start>115 or qh.TOTAL>9500000:
            save();raise SystemExit(0)
        obs=by_visit[v][rank];root=obs['obs_id'];products=product_maps[v]
        product=products.get(root+'_raw.fits')
        row={'root':root,'rank_remaining':rank+1,'obsid':obs['obsid'],'mjd':obs['t_min'],'exptime_caom':obs['t_exptime'],'proposal_id':obs['proposal_id'],'raw_available':product is not None,'ima_available':root+'_ima.fits' in products}
        try:
            if product is None:raise RuntimeError('RAW product absent')
            row.update(qh.header(product['dataURI']))
            k=row['keys'];row['exact_eligible']=(str(k['DETECTOR']).strip()=='IR' and k['SUBARRAY'] is False and str(k['SAMP_SEQ']).strip()=='SPARS50' and isinstance(k['NSAMP'],int) and k['NSAMP']>=8)
        except Exception as exc:row['error']=str(exc)
        out['visits'][v]['records'].append(row);save()
out['status']='complete';save()
