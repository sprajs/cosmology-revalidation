#!/usr/bin/env python3
"""Acquire public research artifacts; record provenance without executing them."""
import argparse
import concurrent.futures
import datetime as dt
import hashlib
import html
import json
from pathlib import Path
import re
import subprocess
import threading
import urllib.request
import gzip
import tarfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
LOCK = threading.Lock()
MANIFEST = ROOT / 'catalog' / 'acquisition.jsonl'

def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def record(entry):
    with LOCK:
        MANIFEST.parent.mkdir(parents=True, exist_ok=True)
        with MANIFEST.open('a') as out:
            out.write(json.dumps(entry, ensure_ascii=False) + '\n')

def fetch(url, relative, kind='source', max_bytes=700_000_000):
    target = ROOT / relative
    if not target.resolve().is_relative_to(ROOT):
        raise ValueError('Output must be inside workspace')
    entry = dict(url=url, path=str(target.relative_to(ROOT)), kind=kind, retrieved_at_utc=now())
    if target.exists():
        print('EXISTS', relative, flush=True)
        return target
    target.parent.mkdir(parents=True, exist_ok=True)
    part = target.with_name(target.name + '.part')
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Supernova-literature-inventory/1.0 (public research collection)'})
        with urllib.request.urlopen(req, timeout=60) as response:
            entry.update(final_url=response.url, http_status=response.status,
                         content_type=response.headers.get('Content-Type'),
                         etag=response.headers.get('ETag'), last_modified=response.headers.get('Last-Modified'))
            if response.status == 206:
                raise ValueError('Partial HTTP response: ' + str(response.headers.get('Content-Range')))
            if int(response.headers.get('Content-Length', 0)) > max_bytes:
                raise ValueError(f'Exceeds per-file acquisition limit {max_bytes}')
            digest, size = hashlib.sha256(), 0
            with part.open('wb') as out:
                while chunk := response.read(1024 * 1024):
                    size += len(chunk)
                    if size > max_bytes:
                        raise ValueError(f'Exceeds per-file acquisition limit {max_bytes}')
                    digest.update(chunk)
                    out.write(chunk)
            expected_size = response.headers.get('Content-Length')
            if expected_size is not None and size != int(expected_size):
                raise ValueError(f'Truncated response: {size} of {expected_size} bytes')
        if target.suffix == '.pdf' and not part.read_bytes()[:1024].lstrip().startswith(b'%PDF-'):
            raise ValueError('Response is not a PDF')
        if target.suffix == '.pdf':
            subprocess.run(['pdfinfo', str(part)], capture_output=True, check=True)
        if target.suffix == '.json':
            json.loads(part.read_bytes())
        if target.name.endswith(('.gz', '.tgz')):
            with gzip.open(part, 'rb') as stream:
                while stream.read(1024 * 1024): pass
        if target.name.endswith(('.tar.gz', '.tgz')):
            with tarfile.open(part) as archive:
                for member in archive: pass
        if target.suffix == '.zip':
            with zipfile.ZipFile(part) as archive:
                if archive.testzip() is not None: raise ValueError('Bad ZIP member CRC')
        part.rename(target)
        entry.update(status='downloaded', bytes=size, sha256=digest.hexdigest())
        record(entry)
        print('OK', relative, size, flush=True)
        return target
    except Exception as exc:
        part.unlink(missing_ok=True)
        entry.update(status='failed', error=str(exc))
        record(entry)
        print('FAILED', relative, str(exc), flush=True)
        return None

def paper(arxiv_id):
    landing = fetch(f'https://arxiv.org/abs/{arxiv_id}', f'papers/metadata/{arxiv_id}.html', 'arxiv_metadata')
    if not landing:
        return
    content = landing.read_text()
    versions = re.findall(r'\[v(\d+)\]', content)
    version = max(map(int, versions), default=1)
    version_id = f'{arxiv_id}v{version}'
    pdf = fetch(f'https://arxiv.org/pdf/{version_id}', f'papers/pdf/{version_id}.pdf', 'paper_pdf')
    if not pdf:
        pdf = fetch(f'https://arxiv.org/pdf/{version_id}?download=1', f'papers/pdf/{version_id}.pdf', 'paper_pdf')
    if not pdf:
        pdf = fetch(f'https://arxiv.org/pdf/{arxiv_id}', f'papers/pdf/{version_id}.pdf', 'paper_pdf_latest_resolved_from_saved_metadata')
    if pdf:
        out = ROOT / 'papers' / 'text' / f'{version_id}.txt'
        out.parent.mkdir(exist_ok=True)
        subprocess.run(['pdftotext', '-layout', str(pdf), str(out)], check=True)
    print('PAPER', version_id, flush=True)

def repository(full_name):
    full_name, _, requested_ref = full_name.partition('@')
    slug = full_name.replace('/', '__') + ('@' + requested_ref if requested_ref else '')
    info = fetch(f'https://api.github.com/repos/{full_name}', f'catalog/repositories/{slug}.json', 'github_repository_metadata')
    if not info:
        return
    metadata = json.loads(info.read_text())
    branch = requested_ref or metadata['default_branch']
    tip = fetch(f'https://api.github.com/repos/{full_name}/commits/{branch}', f'catalog/repositories/{slug}.commit.json', 'github_commit_metadata')
    if not tip:
        return
    sha = json.loads(tip.read_text())['sha']
    fetch(f'https://api.github.com/repos/{full_name}/git/trees/{sha}?recursive=1', f'catalog/repositories/{slug}.tree.json', 'github_tree')
    fetch(f'https://codeload.github.com/{full_name}/tar.gz/{sha}', f'sources/archives/{slug}--{sha}.tar.gz', 'repository_snapshot')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    p = sub.add_parser('papers'); p.add_argument('ids', nargs='+')
    p = sub.add_parser('repos'); p.add_argument('names', nargs='+')
    p = sub.add_parser('url'); p.add_argument('url'); p.add_argument('path'); p.add_argument('--kind', default='source')
    args = parser.parse_args()
    if args.action == 'url':
        fetch(args.url, args.path, args.kind)
    else:
        items = args.ids if args.action == 'papers' else args.names
        operation = paper if args.action == 'papers' else repository
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
            list(pool.map(operation, items))
