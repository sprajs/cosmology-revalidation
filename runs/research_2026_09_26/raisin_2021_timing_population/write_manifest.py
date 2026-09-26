#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
here=Path(__file__).resolve().parent
root=here.parents[2]
files=[p for p in here.rglob('*') if p.is_file() and p.name not in {'manifest.json','write_manifest.py'}]
files += [root/'docs/research-2026-09-26/raisin-2021-timing-population.md',
 root/'runs/research_2026_09_26/astra_design/raisin_timing_assets/author_20211111/nir.FITRES.gz',
 root/'runs/research_2026_09_26/astra_design/raisin_timing_assets/author_20211111/optnir.FITRES.gz',
 root/'runs/research_2026_09_26/raisin_profile_solver_review/native-recovery-design/all500-membership.csv',
 root/'runs/research_2026_09_26/raisin_simulation_execution_review/nominal_fitres/lcplot-feasibility/FITOPT000.LCPLOT.gz',
 root/'runs/research_2026_09_26/astra_design/raisin_timing_assets/historical-coherence-result.json',
 root/'runs/research_2026_09_26/astra_design/raisin_timing_assets/historical-coherence-500.csv']
rows=[]
for f in sorted(set(files)):
 b=f.read_bytes();rows.append({'path':str(f.relative_to(root)),'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()})
(here/'manifest.json').write_text(json.dumps({'scope':'Coherent 2021 RAISIN timing metadata only, amended plot cadence','files':rows},indent=2)+'\n')
