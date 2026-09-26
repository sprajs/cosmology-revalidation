#!/usr/bin/env python3
"""Validate document structure/links and save an indexed, hashed audit record."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re
import subprocess

ROOT=Path(__file__).resolve().parents[2]
DOC=ROOT/'docs/physics-audit'
OUT=ROOT/'runs/physics_audit'
chapters=['luminosity.md','extinction.md','observation-chain.md','cosmology-populations.md']
entries=[]
for name in chapters:
    p=DOC/name;text=p.read_text();section=''
    assert text.count(r'\[')==text.count(r'\]'), name
    inline=re.sub(r'\\\[(.*?)\\\]', '', text, flags=re.S)
    assert inline.count(r'\(')==inline.count(r'\)'), name
    for equation in re.findall(r'\\\[(.*?)\\\]',text,re.S):
        depth=0
        for token in re.finditer(r'(?<!\\)[{}]',equation):
            depth+=1 if token.group()=='{' else -1
            assert depth>=0,(name,equation)
        assert depth==0,(name,equation)
    for n,line in enumerate(text.splitlines(),1):
        if line.startswith('#'):section=line.lstrip('# ').strip()
        for eid in re.findall(r'\\tag\{([A-Z]\d+)\}',line):
            entries.append({'id':eid,'chapter':str(p.relative_to(ROOT)),'line':n,'section':section})
ids=[e['id'] for e in entries];assert len(ids)==len(set(ids))
findings=json.loads((DOC/'findings.json').read_text())
for f in findings['findings']:
    assert set(f.get('equations',[]))<=set(ids), f['id']
(DOC/'equation-index.json').write_text(json.dumps({'format':'LaTeX in Markdown','numbered_equation_groups':len(entries),'scope':'Grouped equations may contain multiple physical equalities; unnumbered supporting coefficients and derivation steps are in the chapters.','equations':entries},indent=2)+'\n')
# Check every Markdown local file link, including the generated links in README.
for p in DOC.glob('*.md'):
    for target in re.findall(r'\]\(([^)]+)\)',p.read_text()):
        if '://' in target or target.startswith('#'):continue
        target=target.split('#')[0]
        linked=(p.parent/target).resolve()
        if linked==OUT/'manifest.json':continue  # generated below after validation
        assert linked.exists(),(str(p),target)

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
# Confirm saved evidence describes the present executable bytes.
for filename in ['cosmology-check.json','extinction-check.json','observation-check.json','luminosity-check.json']:
    data=json.loads((OUT/filename).read_text())
    for key in ['source_sha256','input_sha256','inputs_sha256']:
        for path,digest in data.get(key,{}).items():
            assert sha(ROOT/path)==digest,(filename,path,'changed since verification')
    if 'script_sha256' in data:
        script=ROOT/'scripts/physics_audit'/filename.replace('-check.json','_check.py')
        assert sha(script)==data['script_sha256']

changed=['README.md','docs/first-principles.md','scripts/cosmology/core.py','scripts/standardization/dust_identifiability.py','scripts/phase2/hierarchy/core.py','scripts/phase2/hierarchy/run.py','scripts/phase2/independent_flux/engine.py']
paths=[ROOT/p for p in changed]+list(DOC.glob('*'))+list((ROOT/'scripts/physics_audit').glob('*.py'))+list(OUT.glob('*-check.json'))+list(OUT.glob('*-check.log'))+[OUT/'flux-before.json']+list((OUT/'before').glob('*'))
commands=[f'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 {exe} scripts/physics_audit/{name}_check.py' for name,exe in [('luminosity','phase2/env-official/bin/python'),('extinction','phase2/env-official/bin/python'),('observation','phase2/env-official/bin/python'),('cosmology','.venv/bin/python')]]
manifest={'created_utc':datetime.now(timezone.utc).isoformat(),'purpose':'Physical foundations and bounded implementation audit','git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'numbered_equation_groups':len(entries),'validation':'All four check scripts passed; finalizer verifies source hashes, equation identifiers, display braces and local links. This does not certify empirical assumptions.','commands':commands,'sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(set(paths)) if p.is_file()},'handoff':'docs/physics-audit/HANDOFF.md','historical_results':'Preserved; no new cosmology posterior, original-image photometry reconstruction or production BBC rerun.'}
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps({'numbered_equation_groups':len(entries),'chapters':len(chapters),'findings':len(findings['findings']),'structure_links_source_hashes':'passed','manifest':str((OUT/'manifest.json').relative_to(ROOT))},indent=2))
