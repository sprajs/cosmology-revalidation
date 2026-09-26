"""Offline additive summary; never changes the first-pass result or manifest."""
import csv,hashlib,json
from pathlib import Path

B=Path(__file__).resolve().parent
P=json.loads((B/'protocol.json').read_text())
C=json.loads((B/'caom-results.json').read_text())
H=json.loads((B/'header-results.json').read_text())
T=json.loads((B/'targeted-header-results.json').read_text())
R=json.loads((B/'remaining-header-results.json').read_text())
lookup={v:{x['obs_id']:x for x in C['visits'][v][0]['rows']} for v in ('search','template')}
rows=[];cohorts={};checked={};left={}
for v in ('search','template'):
    allrows=H['visits'][v]['records']+T['visits'][v]['records']+R['visits'][v]['records']
    assert len({x['root'] for x in allrows})==len(allrows)
    for row in allrows:
        o=lookup[v][row['root']];k=row.get('keys',{})
        rows.append({'visit':v,'root':row['root'],'obsid':row['obsid'],'target':o['target_name'],'filter':o['filters'],'proposal_id':o['proposal_id'],'mjd':row['mjd'],'days_from_visit':row['mjd']-P['science_visit_midpoints_mjd'][v],'raw_uri':row.get('uri'),'ima_available':row.get('ima_available'),'samp_seq':str(k.get('SAMP_SEQ','')).strip(),'nsamp':k.get('NSAMP'),'subarray':k.get('SUBARRAY'),'subtype':str(k.get('SUBTYPE','')).strip(),'expstart':k.get('EXPSTART'),'exptime':k.get('EXPTIME'),'sampzero':k.get('SAMPZERO'),'header_sha256':row.get('sha256'),'header_error':row.get('error','')})
    eligible=[r for r in rows if r['visit']==v and r['target']=='DARK' and r['filter']=='BLANK' and r['samp_seq']=='SPARS50' and r['subarray'] is False and isinstance(r['nsamp'],int) and r['nsamp']>=8]
    eligible.sort(key=lambda r:(abs(r['days_from_visit']),r['root']))
    cohorts[v]=eligible[:8]
    checked[v]=len(allrows)
    left[v]=sum(x['t_exptime']>=300 for x in C['visits'][v][0]['rows'])-len(allrows)
with (B/'all-checked-header-ledger.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
with (B/'candidate-cohort.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(cohorts['search']+cohorts['template'])
summary={'status':'metadata feasible conditional on first-eight-prefix timing and RAW processing closure','checked_by_visit':checked,'uninspected_30d_long_exposures_by_visit':left,'candidate_count':{v:len(x) for v,x in cohorts.items()},'candidate_roots':{v:[x['root'] for x in xs] for v,xs in cohorts.items()},'candidate_mode':'WFC3/IR BLANK-filter DARK, full-frame SPARS50 NSAMP16, EXPTIME702.938171 s; science FLT SPARS50 NSAMP8, so first-eight prefix not yet verified','cumulative_metadata_elapsed_seconds':json.loads((B/'targeted-header-results.json').read_text())['cumulative_seconds']+R['elapsed_seconds'],'cumulative_metadata_bytes':json.loads((B/'targeted-header-results.json').read_text())['cumulative_metadata_bytes']+R['metadata_bytes'],'additional_120s_elapsed_seconds':R['elapsed_seconds'],'remaining_stage_complete':R['status']=='complete','no_image_pixels_or_noise_outcomes':True}
(B/'remaining-result.json').write_text(json.dumps(summary,indent=2)+'\n')
manifest={'status':'additive continuation preserved first-pass files','sha256':{n:hashlib.sha256((B/n).read_bytes()).hexdigest() for n in ('remaining-protocol.json','query_remaining.py','remaining-header-results.json','all-checked-header-ledger.csv','candidate-cohort.csv','remaining-result.json')}}
(B/'remaining-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
