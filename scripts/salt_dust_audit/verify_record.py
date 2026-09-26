#!/usr/bin/env python3
"""Check audit-record integrity, not physical adequacy or global bias closure."""
from pathlib import Path
import hashlib
import json
import re
from datetime import datetime, timezone

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/salt_dust_audit'


def sha(path):
    return hashlib.file_digest(path.open('rb'),'sha256').hexdigest()


def main():
    records=[]; errors=[]
    for manifest in sorted(OUT.rglob('manifest.json')):
        content=json.loads(manifest.read_text())
        for kind in ['files','inputs','outputs']:
            entries=content.get(kind,[])
            if isinstance(entries,dict):
                entries=[{'path':key,'sha256':value} for key,value in entries.items()]
            for entry in entries:
                if not isinstance(entry,dict) or not {'path','sha256'}<=entry.keys():
                    errors.append({'manifest':str(manifest.relative_to(ROOT)),'unrecognized_entry':entry})
                    continue
                target=ROOT/entry['path']
                observed=sha(target) if target.is_file() else None
                record={'manifest':str(manifest.relative_to(ROOT)), 'kind':kind,
                        'path':entry['path'],'expected':entry['sha256'],'observed':observed,
                        'matches':observed==entry['sha256']}
                records.append(record)
                if not record['matches']:errors.append(record)
        if manifest.parent.name=='matrix_audit':
            target=ROOT/'scripts/salt_dust_audit/matrix_audit.py'
            if content['script_sha256']!=sha(target):errors.append({'path':str(target),'reason':'matrix script hash mismatch'})
    links=[]
    for path in sorted((ROOT/'docs/salt-dust-audit').glob('*.md')):
        for target in re.findall(r'\]\(([^)]+)\)',path.read_text()):
            if '://' in target or target.startswith('#'):continue
            local=path.parent/target.split('#')[0]
            link={'document':str(path.relative_to(ROOT)),'target':target,'exists':local.exists()}
            links.append(link)
            if not link['exists']:errors.append(link)
    summary={'utc':datetime.now(timezone.utc).isoformat(),'scope':'Recorded-file SHA256 integrity and local documentation links. Does not test missing physical assumptions or certify bias-free analysis.',
             'verifier_sha256':sha(Path(__file__)),'manifests_checked':len(list(OUT.rglob('manifest.json'))),
             'file_hashes_checked':len(records),'local_links_checked':len(links),'errors':errors,
             'file_checks':records,'link_checks':links}
    (OUT/'record-verification.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({key:summary[key] for key in ['manifests_checked','file_hashes_checked','local_links_checked','errors']},indent=2))
    if errors:raise SystemExit(1)


if __name__=='__main__':main()
