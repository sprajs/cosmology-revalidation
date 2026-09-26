#!/usr/bin/env python3
"""Documented 2021 plot-cadence amendment; preserves v1 timing summaries."""
import csv,gzip,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from audit import read_fitres,stats
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
A=json.loads((HERE/'plot-cadence-amendment.json').read_text())
P=json.loads((HERE/'protocol.json').read_text())
V1=json.loads((HERE/'result.json').read_text())
plot=ROOT/A['plot_source'];assert hashlib.sha256(plot.read_bytes()).hexdigest()==A['plot_sha256']
prior=ROOT/A['existing_exact_join'];assert hashlib.sha256((ROOT/'runs/research_2026_09_26/astra_design/raisin_timing_assets/historical-coherence-result.json').read_bytes()).hexdigest()==A['existing_join_result_sha256']
fitpath=ROOT/P['inputs']['nir'];assert hashlib.sha256(fitpath.read_bytes()).hexdigest()==P['inputs']['nir_sha256']
fit,_,_=read_fitres(fitpath)
with (ROOT/P['inputs']['plotted500']).open(newline='') as f:fixed=list(csv.DictReader(f))
ids=[r['CID'] for r in fixed];assert len(ids)==len(set(ids))==500
old={r['CID']:r for r in fixed}
with prior.open(newline='') as f:joined={r['CID']:r for r in csv.DictReader(f)}
assert set(joined)==set(ids)
counts=defaultdict(Counter);all_plot_ids=set();varnames=None
with gzip.open(plot,'rt') as f:
 for line in f:
  if line.startswith('VARNAMES:'):varnames=line.split()[1:]
  elif line.startswith('OBS:'):
   vals=line.split()[1:];assert len(vals)==len(varnames)
   d=dict(zip(varnames,vals));cid=d['CID'];all_plot_ids.add(cid)
   if d['DATAFLAG']=='1':
    assert d['BAND'] in {'J','H'},(cid,d['BAND'])
    counts[cid][d['BAND']]+=1
assert set(ids)==all_plot_ids==set(counts)
ledger=[]
for cid in ids:
 nJ=counts[cid]['J'];nH=counts[cid]['H'];r=old[cid];j=joined[cid]
 assert nJ+nH==int(fit[cid]['NDOF'])+1==int(j['accepted_rows'])
 assert j['NIR_present']=='True' and j['clock_matches']=='True' and j['cadence_matches']=='True'
 assert j['FITRES_SIM_LIBID']==fit[cid]['SIM_LIBID']
 ledger.append({'CID':cid,'n_J_2021plot':nJ,'n_H_2021plot':nH,'zHEL_2021fit':fit[cid]['zHEL'],
                'SIM_LIBID_2021fit':fit[cid]['SIM_LIBID'],'n_J_2022ledger':r['n_J'],
                'n_H_2022ledger':r['n_H'],'zHEL_2022ledger':r['zHEL'],'SIM_LIBID_2022ledger':r['SIM_LIBID']})
with (HERE/'plot-cadence-v2.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(ledger[0]));w.writeheader();w.writerows(ledger)
cadence={'n_J':stats([r['n_J_2021plot'] for r in ledger]),
         'n_H':stats([r['n_H_2021plot'] for r in ledger]),
         'total_J':sum(r['n_J_2021plot'] for r in ledger),
         'total_H':sum(r['n_H_2021plot'] for r in ledger),
         'plot_count_matches_NIR_NDOF_plus1':500,
         'inherited_2022_metadata_mismatch_counts':{
          'n_J':sum(int(r['n_J_2021plot'])!=int(r['n_J_2022ledger']) for r in ledger),
          'n_H':sum(int(r['n_H_2021plot'])!=int(r['n_H_2022ledger']) for r in ledger),
          'zHEL':sum(float(r['zHEL_2021fit'])!=float(r['zHEL_2022ledger']) for r in ledger),
          'SIM_LIBID':sum(r['SIM_LIBID_2021fit']!=r['SIM_LIBID_2022ledger'] for r in ledger)}}
R=dict(V1)
R['fixed_plot_cadence']=cadence
R['v2_amendment_sha256']=hashlib.sha256((HERE/'plot-cadence-amendment.json').read_bytes()).hexdigest()
R['v1_result_sha256']=hashlib.sha256((HERE/'result.json').read_bytes()).hexdigest()
(HERE/'result-v2.json').write_text(json.dumps(R,indent=2,sort_keys=True)+'\n')
with (HERE/'first8-v2.csv').open('w',newline='') as f:
 extra={r['CID']:r for r in ledger}
 with (HERE/'first8.csv').open(newline='') as src:old8=list(csv.DictReader(src))
 fields=list(old8[0]);fields.remove('n_J');fields.remove('n_H');fields=fields[:3]+['n_J_2021plot','n_H_2021plot']+fields[3:]
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
 for r in old8:
  out={k:v for k,v in r.items() if k not in {'n_J','n_H'}}
  out['n_J_2021plot']=extra[r['CID']]['n_J_2021plot'];out['n_H_2021plot']=extra[r['CID']]['n_H_2021plot'];w.writerow(out)
print(json.dumps(cadence,indent=2))
