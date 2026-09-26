from pathlib import Path
import csv,gzip,json,hashlib,collections
import numpy as np
R=Path.cwd();P=Path(__file__).parent;S=R/'runs/research_2026_09_26/raisin_simulation_execution_review/nominal_fitres'
proto=P/'metadata-match-protocol.json'
fit=[]
with gzip.open(S/'nir.FITRES.gz','rt') as f:
 for line in f:
  a=line.split()
  if not a:continue
  if a[0]=='VARNAMES:':keys=a[1:]
  if a[0]=='SN:':fit.append(dict(zip(keys,a[1:])))
fids=np.array([int(r['CID']) for r in fit]);peaks=np.array([float(r['PKMJD']) for r in fit]);libs=np.array([int(r['SIM_LIBID']) for r in fit]);rows=[]
for r in csv.DictReader((P/'cadence-500.csv').open()):
 ids=[int(x) for x in r['cadence_candidates_0p0056'].split(',') if x]
 known=len(ids)==1; candidate=(abs(peaks-float(r['model_clock']))<=.0056)&(libs==ids[0]) if known else np.zeros(len(fit),bool)
 ix=fids[candidate]; rows.append({'plot_CID':int(r['CID']),'model_clock':float(r['model_clock']),'identified_SIMLIB':ids[0] if known else None,'cadence_candidate_count':len(ids),'FITRES_candidate_count':len(ix),'original_CID_retained':int(int(r['CID']) in ix),'candidate_CIDs':'|'.join(map(str,ix))})
with (P/'metadata-candidates-500.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
res={'protocol_sha256':hashlib.sha256(proto.read_bytes()).hexdigest(),'plot_objects':len(rows),'FITRES_objects':len(fit),'unknown_cadence':sum(r['identified_SIMLIB'] is None for r in rows),'candidate_counts':{'zero':sum(r['FITRES_candidate_count']==0 for r in rows),'unique':sum(r['FITRES_candidate_count']==1 for r in rows),'multiple':sum(r['FITRES_candidate_count']>1 for r in rows)},'candidate_count_range':[min(r['FITRES_candidate_count'] for r in rows),max(r['FITRES_candidate_count'] for r in rows)],'candidate_count_median':float(np.median([r['FITRES_candidate_count'] for r in rows])),'original_CID_retained':sum(r['original_CID_retained'] for r in rows),'interpretation':'No event identity or execution linkage inferred. No distance or objective field entered matching.'}
(P/'metadata-match-result.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
