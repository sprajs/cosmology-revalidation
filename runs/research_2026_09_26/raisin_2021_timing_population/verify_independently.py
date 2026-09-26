#!/usr/bin/env python3
"""Independent FITRES parser and standard-library checks; no cohort changes."""
import csv,gzip,hashlib,json,statistics
from pathlib import Path
here=Path(__file__).resolve().parent
root=here.parents[2]
p=json.loads((here/'protocol.json').read_text())
r=json.loads((here/'result-v2.json').read_text())
def load(path):
 data={};keys=None
 with gzip.open(root/path,'rt') as f:
  for line in f:
   if line.startswith('VARNAMES:'):keys=line.split()[1:]
   elif line.startswith('SN:'):
    values=line.split()[1:];assert len(values)==len(keys)
    d=dict(zip(keys,values));assert d['CID'] not in data;data[d['CID']]=d
 return data
n=load(p['inputs']['nir']);j=load(p['inputs']['optnir']);c=sorted(set(n)&set(j),key=int)
assert len(n)==30000 and len(j)==len(c)==29995
x=[float(n[i]['PKMJDINI'])-float(n[i]['SIM_PKMJD']) for i in c]
y=[float(j[i]['PKMJD'])-float(j[i]['SIM_PKMJD']) for i in c]
z=[float(j[i]['PKMJD'])-float(n[i]['PKMJD']) for i in c]
expected=r['summaries']['all_common'];checks={}
for key,a in [('initializer_minus_sim_peak_day',x),('joint_peak_minus_sim_peak_day',y),('joint_peak_minus_nir_peak_day',z)]:
 checks[key]={'n':len(a),'mean':statistics.mean(a),'sample_sd':statistics.stdev(a),'median':statistics.median(a)}
 for metric in ['mean','sample_sd','median']:
  assert abs(checks[key][metric]-expected[key][metric])<1e-10,(key,metric)
with (root/p['inputs']['plotted500']).open(newline='') as f: ids=[a['CID'] for a in csv.DictReader(f)]
assert len(ids)==len(set(ids))==500 and all(i in c for i in ids)
out={'scope':'Independent parser and standard-library arithmetic on same frozen 2021 cohort; no additional selection','checks':checks,'fixed500_unique_and_common':True,'result_v2_sha256':hashlib.sha256((here/'result-v2.json').read_bytes()).hexdigest()}
(here/'independent-verification.json').write_text(json.dumps(out,indent=2)+'\n')
print('independent arithmetic/join PASS')
