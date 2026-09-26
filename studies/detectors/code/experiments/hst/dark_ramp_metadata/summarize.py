"""Offline summary of immutable metadata-only MAST query ledgers."""
import csv,hashlib,json
from collections import Counter
from pathlib import Path

B=Path(__file__).resolve().parent
caom=json.loads((B/'caom-results.json').read_text())
head=json.loads((B/'header-results.json').read_text())
target=json.loads((B/'targeted-header-results.json').read_text())
expanded=json.loads((B/'expanded-caom-results.json').read_text())
rows=[]
for kind,source in [('date_ranked',head),('352s_diagnostic',target)]:
    for visit,entry in source['visits'].items():
        for rank,r in enumerate(entry['records'],1):
            k=r.get('keys',{})
            rows.append({'stage':kind,'visit':visit,'rank':rank,'root':r['root'],'obsid':r['obsid'],'mjd':r['mjd'],'days_from_science':r['mjd']-({'search':57715.6503,'template':58077.6242}[visit]),'proposal_id':r.get('proposal_id'),'target':r.get('target'),'caom_exptime':r['exptime_caom'],'raw_uri':r.get('uri'),'raw_available':r.get('raw_available'),'ima_available':r.get('ima_available'),'detector':str(k.get('DETECTOR','')).strip(),'subarray':k.get('SUBARRAY'),'subtype':str(k.get('SUBTYPE','')).strip(),'samp_seq':str(k.get('SAMP_SEQ','')).strip(),'nsamp':k.get('NSAMP'),'header_exptime':k.get('EXPTIME'),'expstart':k.get('EXPSTART'),'header_sha256':r.get('sha256'),'header_bytes':r.get('header_bytes'),'error':r.get('error','')})
with (B/'candidate-header-ledger.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
summary={'scope':'public metadata and FITS primary headers only, no image pixels or FITS data arrays','caom_30d_counts':{v:len(z[0]['rows']) for v,z in caom['visits'].items()},'caom_180d_counts':{v:len(z['rows']) for v,z in expanded['visits'].items()},'header_rows':len(rows),'header_by_stage_visit':{},'exact_eligible_count':0,'cumulative_seconds':target['cumulative_seconds'],'cumulative_metadata_bytes':target['cumulative_metadata_bytes'],'constraint':'nearest 40 exptime>=300 per visit plus eight per visit in repeated 352.939514-s diagnostic; remaining 30d/180d headers uninspected under cap'}
for kind in ('date_ranked','352s_diagnostic'):
    for visit in ('search','template'):
        subset=[r for r in rows if r['stage']==kind and r['visit']==visit]
        modes=Counter((r['samp_seq'],r['nsamp'],r['subarray']) for r in subset)
        summary['header_by_stage_visit'][kind+'_'+visit]={'count':len(subset),'modes':[{'samp_seq':a,'nsamp':b,'subarray':c,'count':n} for (a,b,c),n in sorted(modes.items(),key=lambda x:str(x[0]))],'raw_product_count':sum(bool(r['raw_available']) for r in subset),'ima_product_count':sum(bool(r['ima_available']) for r in subset),'range_errors':sum(bool(r['error']) for r in subset)}
        summary['exact_eligible_count']+=sum(r['detector']=='IR' and r['subarray'] is False and r['samp_seq']=='SPARS50' and isinstance(r['nsamp'],int) and r['nsamp']>=8 for r in subset)
(B/'result.json').write_text(json.dumps(summary,indent=2)+'\n')
manifest={'status':'offline summary complete','sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(B.iterdir()) if p.is_file() and p.name not in ('manifest.json',)}}
(B/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
