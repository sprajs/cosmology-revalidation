#!/usr/bin/env python3
"""Resolve public Git LFS pointers and verify their upstream SHA-256 OIDs."""
import concurrent.futures
import hashlib
import json
from pathlib import Path
import re
import shutil
from urllib.parse import quote
from collect import ROOT, fetch, now

def hydrate(item):
    slug, revision, relative, oid, size = item
    repo = slug.split('@')[0].replace('__', '/')
    url = f'https://media.githubusercontent.com/media/{repo}/{revision}/{quote(relative)}'
    object_path = f'data/lfs/{oid}'
    p = fetch(url, object_path, 'git_lfs_object', max_bytes=1_000_000_000)
    result = dict(repository=slug, revision=revision, path=relative, sha256=oid, bytes=size, url=url)
    if p and p.stat().st_size == size and hashlib.file_digest(p.open('rb'), 'sha256').hexdigest() == oid:
        shutil.copyfile(p, ROOT/'sources/repos'/slug/relative)
        result['status'] = 'hydrated_verified'
    else:
        result['status'] = 'unresolved'
    return result

if __name__ == '__main__':
    tasks = []
    for r in json.loads((ROOT/'catalog/repository_inventory.json').read_text()):
        for relative in r['lfs_pointers']:
            p = ROOT/r['directory']/relative
            b = p.read_bytes()
            if not b.startswith(b'version https://git-lfs.github.com/spec/v1'): continue
            oid = re.search(rb'oid sha256:(\w+)', b)[1].decode()
            size = int(re.search(rb'size (\d+)', b)[1])
            tasks.append((r['repository'], r['revision'], relative, oid, size))
    # Repositories can reference the same object: avoid concurrent writes to it.
    unique = {}; aliases = []
    for t in tasks:
        if t[3] in unique: aliases.append(t)
        else: unique[t[3]] = t
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
        results = list(pool.map(hydrate, unique.values()))
    results.extend(hydrate(t) for t in aliases)
    dest = ROOT/'catalog/lfs_hydration.json'
    old = json.loads(dest.read_text())['objects'] if dest.exists() else []
    by_path = {(r['repository'], r['path']): r for r in old + results}
    dest.write_text(json.dumps(dict(checked_at_utc=now(), objects=list(by_path.values())), indent=2)+'\n')
    print('Hydrated', sum(r['status']=='hydrated_verified' for r in results), 'of', len(results))
