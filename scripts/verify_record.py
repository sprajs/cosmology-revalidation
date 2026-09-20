#!/usr/bin/env python3
"""Verify active scientific manifests and local documentation links, read-only.

Historical failed/preliminary manifests are intentionally not statements about
the current output paths. Their exclusion is explicit in the saved result.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
HISTORICAL = {'preliminary', 'unsuccessful', 'initial-results', 'validation-initial-16'}
HASH = re.compile(r'^[0-9a-f]{64}$')
CACHE = {}

def digest(p):
    if p not in CACHE:
        h = hashlib.sha256()
        with p.open('rb') as stream:
            while block := stream.read(1024 * 1024):
                h.update(block)
        CACHE[p] = h.hexdigest()
    return CACHE[p]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='runs/verification/record-check.json')
    args = parser.parse_args()
    archive = {}
    for base in [ROOT/'scripts', ROOT/'docs/experiments', ROOT/'runs/cosmology/provenance', ROOT/'runs/age_signal/provenance', ROOT/'runs/directional/provenance']:
        if base.exists():
            for p in base.rglob('*'):
                if p.is_file() and p.suffix in {'.py', '.md'}:
                    archive.setdefault(digest(p), []).append(str(p.relative_to(ROOT)))
    checked = []; historical = []; failures = []; snapshots = []; listed = []
    def check(path, expected, manifest, kind):
        p = ROOT/path
        if p.is_file() and digest(p) == expected:
            checked.append({'manifest':manifest,'kind':kind,'path':path,'sha256':expected})
        elif (path.startswith('scripts/') or path.startswith('docs/experiments/') or kind=='code') and expected in archive:
            snapshots.append({'manifest':manifest,'requested_path':path,'sha256':expected,'exact_snapshot':archive[expected]})
        else:
            failures.append({'manifest':manifest,'kind':kind,'path':path,'expected':expected,'actual':digest(p) if p.is_file() else None})
    manifests=set((ROOT/'runs').rglob('*manifest*.json')) | set((ROOT/'runs/directional').glob('experiment-*.json'))
    for p in sorted(manifests):
        rel = str(p.relative_to(ROOT))
        if HISTORICAL.intersection(p.parts):
            historical.append(rel);continue
        d = json.loads(p.read_text())
        if not isinstance(d,dict):
            historical.append(rel+' (source acquisition list, not an experiment manifest)');continue
        for key in ['inputs','inputs_sha256','input_sha256','outputs','outputs_sha256','output_sha256','code_sha256','record_files_sha256']:
            value=d.get(key)
            if isinstance(value,dict):
                for path,h in value.items():
                    if isinstance(h,str) and HASH.fullmatch(h):check(path,h,rel,key)
            elif key=='outputs' and isinstance(value,list):
                for path in value:
                    if not isinstance(path,str):continue
                    f=ROOT/path
                    if not f.is_file():failures.append({'manifest':rel,'kind':'listed_output_missing','path':path})
                    else:listed.append({'manifest':rel,'path':path,'current_sha256':digest(f),'notice':'Original manifest listed path only; final inventory records current digest.'})
        for key in ['code_sha256','script_sha256','plan_sha256','lock_sha256']:
            h=d.get(key)
            if isinstance(h,str) and HASH.fullmatch(h):
                if key=='lock_sha256':check('uv.lock',h,rel,key)
                elif h in archive:snapshots.append({'manifest':rel,'kind':key,'sha256':h,'exact_snapshot':archive[h]})
                else:failures.append({'manifest':rel,'kind':key,'unresolved_sha256':h})
        def structured_records(value):
            if isinstance(value,dict):
                if isinstance(value.get('path'),str) and isinstance(value.get('sha256'),str) and HASH.fullmatch(value['sha256']):
                    check(value['path'],value['sha256'],rel,'structured_record')
                for nested in value.values():structured_records(nested)
            elif isinstance(value,list):
                for nested in value:structured_records(nested)
        structured_records(d)
    bad_links=[];links=0
    for p in [ROOT/'README.md',*(ROOT/'docs').rglob('*.md')]:
        prose=re.sub(r'```.*?```|`[^`\n]*`','',p.read_text(),flags=re.S)
        for raw in re.findall(r'\]\(([^)]+)\)',prose):
            target=raw.split('#')[0].strip('<>')
            if not target or '://' in target or target.startswith('mailto:'):continue
            if ' ' in target and not raw.startswith('<'):continue
            dest=(p.parent/target).resolve() if not target.startswith('/') else Path(target)
            links+=1
            if not dest.exists():bad_links.append({'document':str(p.relative_to(ROOT)),'target':target})
    report={'created_utc':datetime.now(timezone.utc).isoformat(),'checked_digest_records':len(checked),'resolved_code_plan_snapshots':snapshots,'listed_outputs_now_hashed':listed,'historical_manifests_excluded':historical,'checked_local_links':links,'failures':failures,'broken_local_links':bad_links,'checked':checked}
    out=ROOT/args.output;out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['checked_digest_records','checked_local_links','failures','broken_local_links']},indent=2))
    raise SystemExit(bool(failures or bad_links))

if __name__=='__main__':main()
