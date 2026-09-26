"""Independent local metadata aggregation only; no queries, downloads or pixels."""
from pathlib import Path
from collections import Counter
import json,hashlib
O=Path(__file__).resolve().parent
M=O.parents[2]/'dark_ramp_metadata'
# M is corrected by locating the shared feasibility root, without changing it.
M=next(p for p in O.parents if p.name=='raisin_hst_pixel_feasibility')/'dark_ramp_metadata'
names=['protocol.json','products-protocol.json','headers-protocol.json','expand-protocol.json','caom-results.json','product-results.json','header-results.json','expanded-caom-results.json']
a={n:json.loads((M/n).read_text()) for n in names}
r={'scope':'Local metadata only. No native fits, pixels, duplicate network requests or download authorization.','visits':{},'inputs':[{'path':str(M/n),'sha256':hashlib.sha256((M/n).read_bytes()).hexdigest()} for n in names]}
for visit,v in a['header-results.json']['visits'].items():
 rows=v['records'];pr=a['product-results.json']['visits'][visit]
 c=Counter((x['keys']['SAMP_SEQ'].strip(),str(x['keys']['SUBARRAY']),int(x['keys']['NSAMP'])) for x in rows)
 step=[]
 for x in rows:
  if x['keys']['SAMP_SEQ'].strip()=='STEP50' and not x['keys']['SUBARRAY']:
   pp=[{'URI':p['dataURI'],'bytes':p['size'],'subgroup':p['productSubGroupDescription']} for p in pr['products'] if p['productFilename'].startswith(x['root']) and p['productSubGroupDescription'] in ['RAW','IMA','FLT']]
   step.append({'root':x['root'],'EXPSTART':x['keys']['EXPSTART'],'distance_to_visit_days':abs(x['keys']['EXPSTART']-a['protocol.json']['science_visit_midpoints_mjd'][visit]),'NSAMP':x['keys']['NSAMP'],'EXPTIME':x['keys']['EXPTIME'],'proposal_id':x['proposal_id'],'products':pp})
 expanded=a['expanded-caom-results.json']['visits'][visit]['rows']
 r['visits'][visit]={'headers_checked':len(rows),'exact_eligible':sum(x['exact_eligible'] for x in rows),'sequence_layout_counts':[{'SAMP_SEQ':s,'SUBARRAY':l,'NSAMP':n,'count':v} for (s,l,n),v in sorted(c.items())],'nearer_long_exposure_CAOM_roots':pr['eligible_caom_exptime_count'],'nearer_unqueried_headers':pr['eligible_caom_exptime_count']-len(rows),'full_STEP50':step,'expanded_CAOM_count':len(expanded),'expanded_exact_DARK_count':sum(x['target_name']=='DARK' for x in expanded),'expanded_targets':dict(Counter(x['target_name'] for x in expanded)),'expanded_dark_702_second_candidates_not_sampling_verified':sum(x['target_name']=='DARK' and 702<x['t_exptime']<704 for x in expanded)}
(O/'fallback-metadata-result.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps({k:{q:v[q] for q in ['headers_checked','exact_eligible','nearer_unqueried_headers','expanded_CAOM_count','expanded_dark_702_second_candidates_not_sampling_verified']} for k,v in r['visits'].items()},indent=2))
