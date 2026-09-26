from pathlib import Path
import csv,json,hashlib,collections
R=Path.cwd();p=R/'runs/research_2026_09_26/raisin_hst_pixel_feasibility/dark_ramp_metadata';o=Path(__file__).resolve().parent
sha=lambda x:hashlib.sha256(x.read_bytes()).hexdigest()
manifest=json.loads((p/'remaining-manifest.json').read_text())
assert all(sha(p/n)==h for n,h in manifest['sha256'].items())
records={}
files=['header-results.json','targeted-header-results.json','remaining-header-results.json']
for n in files:
 d=json.loads((p/n).read_text())
 for visit,v in d['visits'].items():
  for row in v['records']:
   key=(visit,row['root'])
   if key in records: assert records[key]['keys']==row['keys']
   records[key]=row
candidate={k:r for k,r in records.items() if r.get('keys',{}).get('SAMP_SEQ','').strip()=='SPARS50' and r['keys'].get('NSAMP',0)>=8 and r['keys'].get('SUBARRAY') is False}
ledger=list(csv.DictReader((p/'candidate-cohort.csv').open()))
assert set(candidate)=={(r['visit'],r['root']) for r in ledger}
for r in ledger:
 x=candidate[r['visit'],r['root']];k=x['keys']
 assert r['header_sha256']==x['sha256'] and r['raw_uri']==x['uri']
 assert int(r['nsamp'])==k['NSAMP']==16
 assert float(r['exptime'])==k['EXPTIME']==702.938171
 assert float(r['sampzero'])==k['SAMPZERO']==2.911756
 assert r['target']=='DARK' and r['filter']=='BLANK'
result={'status':'PASS_SAVED_METADATA_COHORT_REVIEW','independent_queries':0,'pixels_read':0,'checked_roots_by_visit':dict(collections.Counter(v for v,r in records)), 'matched_candidate_count':dict(collections.Counter(v for v,r in candidate)), 'candidate_roots':{v:sorted(r for vv,r in candidate if vv==v) for v in ['search','template']},'matching_science_sequence':'SPARS50 full-frame, but candidate NSAMP16 versus science NSAMP8; no exact per-read prefix established by primary header','source_hashes':{n:sha(p/n) for n in files+['candidate-cohort.csv','remaining-result.json','remaining-protocol.json','remaining-manifest.json']}}
(o/'matched-candidate-result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['status','checked_roots_by_visit','matched_candidate_count']}))
