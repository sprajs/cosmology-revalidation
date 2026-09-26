#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
here=Path(__file__).resolve().parent
root=here.parents[2]
files=[p for p in here.rglob('*') if p.is_file() and p.name not in {'manifest.json','write_manifest.py'}]
base=root/'runs/research_2026_09_26/astra_design/raisin_timing_assets'
files += [base/'rounding_2021/precision-ledger.json',
 base/'rounding_2021/fits/endpoints/case-ledger.json',base/'rounding_2021/fits/corners/case-ledger.json',
 base/'rounding_2021/fits/endpoints/fit.log',base/'rounding_2021/fits/corners/fit.log',
 base/'rounding_2021/fits/endpoints/fit.FITRES.TEXT',base/'rounding_2021/fits/corners/fit.FITRES.TEXT',
 base/'rounding_2021/endpoint-result.json',base/'rounding_2021/corner-result.json',
 base/'instrumentation_2021/fits/absent/fit.log',base/'instrumentation_2021/fits/absent/fit.FITRES.TEXT',
 base/'author_20211111/nir.FITRES.gz']
files += list((base/'instrumentation_2021/fits/absent/data/DES_RAISIN_SIM').glob('*.DAT'))
for stage in ['endpoints','corners']:
 files += list((base/'rounding_2021/fits'/stage/'data/DES_RAISIN_SIM').glob('*.DAT'))
rows=[]
for p in sorted(set(files)):
 b=p.read_bytes();rows.append({'path':str(p.relative_to(root)),'sha256':hashlib.sha256(b).hexdigest(),'bytes':len(b)})
(here/'manifest.json').write_text(json.dumps({'scope':'Independent source and native-output verification of frozen 2021 precision fits','files':rows},indent=2)+'\n')
