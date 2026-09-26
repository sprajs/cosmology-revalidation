#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
r=Path(__file__).resolve().parents[3];o=r/'phase2/data_audit'
files=[p for p in o.iterdir() if p.is_file() and p.name!='outputs_manifest.json']
files+=list((r/'scripts/phase2/data_audit').glob('*.py'))
p=r/'docs/phase2/data-audit.md'
if p.exists():files.append(p)
x=[{'path':str(p.relative_to(r)),'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest()} for p in sorted(files)]
(o/'outputs_manifest.json').write_text(json.dumps({'files':x},indent=2)+'\n')
print('Hashed',len(x),'outputs/code/report files')
