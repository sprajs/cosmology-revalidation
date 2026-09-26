#!/usr/bin/env python3
"""Nontrivial validation: finite normalization and exact trigger-state enumeration."""
import itertools,json
import numpy as np
from selection_assets import PopulationPDF,HostEfficiency,DetectionEfficiency,ROOT
p=PopulationPDF();d=DetectionEfficiency();h=HostEfficiency();checks=[]
for name,(names,axes,_,_) in p.maps.items():
 for mass in [7.0,9.0,10.0,12.7]:
  context={'LOGMASS':mass,'ZTRUE':.6}
  values=p.density(name,axes[0],**context)
  integral=np.trapezoid(values,axes[0]);assert abs(integral-1)<1e-12
 checks.append({'check':'density_normalization','name':name,'pass':True})
# Exposure groups [0,.1], [.7,.8], [2.] by historical first-epoch grouping.
mjd=np.array([0,.1,.7,.8,2.]);bands=np.array(list('grizi'));snr=np.array([3.,6.,5.,2.,8.])
probs=np.array([d.probability(b,s) for b,s in zip(bands,snr)])
answer=0.
for state in itertools.product([0,1],repeat=5):
 if sum([any(state[:2]),any(state[2:4]),state[4]])>=2:
  answer+=np.prod([pr if bit else 1-pr for pr,bit in zip(probs,state)])
result=d.trigger_probability(mjd,bands,snr);assert abs(result-answer)<1e-14
checks.append({'check':'trigger_against_complete_32_state_enumeration','pass':True,'probability':result})
assert d.probability('g',0)==0 and d.probability('g',-1)==0
assert d.trigger_probability([1.],['g'],[100.])==0
for field in ['C1','C2','C3','E1','E2','S1','S2','X1','X2','X3']:
 for mjd in [57000,57700,58000]:
  for r in [14.,21.,23.,26.,35.]:
   for gr in [-4.,1.,2.,80.]:assert 0<=h.probability(r,gr,field,mjd)<=1
checks.append({'check':'host_maps_all_fields_seasons_boundary_probabilities','pass':True})
report={'checks':checks,'limits':['Tests validate this reader and the conditional trigger combinatorics, not full historical SNANA parity, generated density source equality, or final fitted selection.']}
(ROOT/'phase2/literature/selection-validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
