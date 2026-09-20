#!/usr/bin/env python3
"""Index and check acquired files. Does not run scientific analysis."""
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import subprocess
import tarfile
import gzip
import zipfile
from collect import ROOT, record, now

class Meta(HTMLParser):
    def __init__(self):
        super().__init__(); self.items = {}
    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == 'meta' and a.get('name', '').startswith('citation_'):
            self.items.setdefault(a['name'][9:], []).append(a.get('content', ''))

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        while chunk := f.read(1024*1024): h.update(chunk)
    return h.hexdigest()

def extract_repositories():
    inventory = []
    for archive in sorted((ROOT/'sources/archives').glob('*.tar.gz')):
        slug, revision = archive.name[:-7].split('--')
        target = ROOT/'sources/repos'/slug
        target.mkdir(parents=True, exist_ok=True)
        count, skipped = 0, []
        with tarfile.open(archive) as tf:
            for item in tf:
                parts = item.name.split('/', 1)
                if len(parts) != 2 or not parts[1]: continue
                item.name = parts[1]
                if item.issym() or item.islnk():
                    skipped.append({'path':item.name,'target':item.linkname,'reason':'Nonportable link; archived unchanged, not materialized'})
                    continue
                dest = target/item.name
                if not dest.resolve().is_relative_to(target.resolve()): raise ValueError('Unsafe archive member')
                if not dest.exists(): tf.extract(item, target, filter='data')
                if item.isfile(): count += 1
        lfs = []
        for p in target.rglob('*'):
            if p.is_file() and p.stat().st_size < 1024:
                if p.read_bytes().startswith(b'version https://git-lfs.github.com/spec/v1'):
                    lfs.append(str(p.relative_to(target)))
        inventory.append(dict(repository=slug,revision=revision,archive=str(archive.relative_to(ROOT)),directory=str(target.relative_to(ROOT)),regular_files=count,skipped_links=skipped,lfs_pointers=lfs))
        print(slug, count, 'files;', len(skipped), 'links;',len(lfs),'LFS pointers',flush=True)
    (ROOT/'catalog/repository_inventory.json').write_text(json.dumps(inventory, indent=2)+'\n')

def index_papers():
    records=[]
    for path in sorted((ROOT/'papers/metadata').rglob('*.html')):
        p=Meta(); p.feed(path.read_text())
        aid=str(path.relative_to(ROOT/'papers/metadata').with_suffix(''))
        if not p.items.get('arxiv_id'): continue
        pdfs=sorted((ROOT/'papers/pdf').glob(aid+'v*.pdf'))
        records.append(dict(arxiv_id=aid,metadata=p.items,landing_path=str(path.relative_to(ROOT)),pdf_paths=[str(p.relative_to(ROOT)) for p in pdfs]))
    (ROOT/'catalog/papers.json').write_text(json.dumps(records,ensure_ascii=False,indent=2)+'\n')
    bib=[]
    for p in records:
        m=p['metadata']; aid=p['arxiv_id']
        title=m.get('title',[''])[0].replace('{','').replace('}','')
        authors=' and '.join(m.get('author',[]))
        year=m.get('date',[''])[0][:4]
        bib.append(f'@misc{{arxiv{aid.replace(".", "")},\n  title = {{{title}}},\n  author = {{{authors}}},\n  year = {{{year}}},\n  eprint = {{{aid}}},\n  archivePrefix = {{arXiv}},\n  url = {{https://arxiv.org/abs/{aid}}}\n}}')
    (ROOT/'catalog/references.bib').write_text('\n\n'.join(bib)+'\n')

def check():
    latest={}
    for line in (ROOT/'catalog/acquisition.jsonl').read_text().splitlines():
        r=json.loads(line)
        latest[r['path']]=r
    errors=[]; checked=0
    for relative,r in latest.items():
        if r['status'] != 'downloaded': continue
        p=ROOT/relative
        if not p.is_file(): errors.append(f'Missing: {relative}'); continue
        if p.stat().st_size!=r['bytes'] or digest(p)!=r['sha256']:
            errors.append(f'Checksum or size mismatch: {relative}')
        checked+=1
    pdf_pages={}
    for p in (ROOT/'papers/pdf').rglob('*.pdf'):
        r=subprocess.run(['pdfinfo',str(p)],capture_output=True,text=True)
        if r.returncode: errors.append(f'Unreadable PDF: {p.name}')
        else:
            pdf_pages[p.name]=next((x.split(':',1)[1].strip() for x in r.stdout.splitlines() if x.startswith('Pages:')),None)
    containers=[]
    for relative,r in latest.items():
        if r['status'] != 'downloaded': continue
        p=ROOT/relative
        try:
            if p.name.endswith(('.gz','.tgz')):
                with gzip.open(p,'rb') as f:
                    while f.read(1024*1024): pass
                containers.append(relative)
            elif p.suffix=='.zip':
                with zipfile.ZipFile(p) as f:
                    if f.testzip() is not None: raise ValueError('Bad CRC')
                containers.append(relative)
        except Exception as exc: errors.append(f'Invalid compressed file: {relative}: {exc}')
    report=dict(checked_at_utc=now(),scope='File integrity and inventory only; no cosmological fits or scientific replication',verified_downloads=checked,pdf_pages=pdf_pages,verified_compressed_files=containers,errors=errors)
    (ROOT/'catalog/integrity.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    if errors: raise SystemExit(1)

if __name__=='__main__':
    extract_repositories(); index_papers(); check()
