#!/usr/bin/env python3
import collections,hashlib,json
from pathlib import Path
from audit import checked,read_fitres
here=Path(__file__).resolve().parent
p=json.loads((here/'protocol.json').read_text())
q=json.loads((here/'nuisance-protocol.json').read_text())
out={'nuisance_protocol_sha256':hashlib.sha256((here/'nuisance-protocol.json').read_bytes()).hexdigest(),'tables':{}}
for label,prefix in [('nir','nir'),('optnir','optnir')]:
 path=checked(p['inputs'][prefix],p['inputs'][prefix+'_sha256'])
 rows,_,_=read_fitres(path)
 table={}
 for k in q['fields']:
  a=[float(r[k]) for r in rows.values()]
  c=collections.Counter(r[k] for r in rows.values())
  table[k]={'distinct_printed':len(c),'min':min(a),'max':max(a),'zero_count':sum(x==0 for x in a),'negative_count':sum(x<0 for x in a),'printed_mode':c.most_common(1)[0][0],'mode_count':c.most_common(1)[0][1]}
 out['tables'][label]=table
(here/'nuisance-fields.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2))
