#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
here=Path(__file__).resolve().parent
root=here.parents[3]
paths=[p for p in here.rglob('*') if p.is_file() and p.name not in {'manifest.json','write_manifest.py'}]
paths.append(root/'docs/research-2026-09-26/csp-wirc-operator-vega-bridge.md')
files=[]
for p in sorted(paths):
 b=p.read_bytes();files.append({'path':str(p.relative_to(root)),'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()})
(here/'manifest.json').write_text(json.dumps({'purpose':'Pinned SNooPy VegaB reference bridge; earlier CSP branch unmodified','files':files},indent=2)+'\n')
