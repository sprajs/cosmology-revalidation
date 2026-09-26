#!/usr/bin/env python3
"""Post-verification comparison; never used to derive independent native states."""
import json,hashlib
from pathlib import Path
here=Path(__file__).resolve().parent
root=here.parents[2]
base=root/'runs/research_2026_09_26/astra_design/raisin_timing_assets/rounding_2021'
ours=json.loads((here/'case-results.json').read_text())
out={'scope':'Post-independent-result comparison against Astra summaries; no source adjustment','stages':{}}
for stage,name in [('endpoints','endpoint-result.json'),('corners','corner-result.json')]:
 a=json.loads((base/name).read_text()); b={x['CID']:x for x in ours[stage]}
 assert set(b)=={x['new_CID'] for x in a['cases']}
 maxdiff={k:0. for k in ['D','dataQ','deltaD','deltaQ','dataQ_cap']}
 gate_mismatch=[]
 for x in a['cases']:
  y=b[x['new_CID']]
  pairs={'D':'native_D','dataQ':'native_data_Q','deltaD':'delta_D','deltaQ':'delta_data_Q','dataQ_cap':'data_Q_cap'}
  for k,v in pairs.items(): maxdiff[k]=max(maxdiff[k],abs(y[k]-x[v]))
  if ('distance_cap' not in y['issues'])!=x['distance_cap_pass'] or ('dataQ_cap' not in y['issues'])!=x['data_Q_cap_pass']:
   gate_mismatch.append(x['new_CID'])
 assert not gate_mismatch
 out['stages'][stage]={'cases':len(b),'astra_result_sha256':hashlib.sha256((base/name).read_bytes()).hexdigest(),
  'max_abs_numeric_difference':maxdiff,'cap_gate_mismatch_ids':gate_mismatch,'astra_all_gates_pass':a['all_gates_pass']}
(here/'author-result-crosscheck.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
