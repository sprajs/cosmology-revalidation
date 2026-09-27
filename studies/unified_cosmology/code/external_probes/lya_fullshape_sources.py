"""Acquire versioned Ly-alpha papers and record the checked public interfaces."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parents[4]
TREE_COMMIT = 'bb0c1c9009dc76d1391300e169e8df38fd1096db'
SOURCES = {
    'paper.html': 'https://arxiv.org/html/2607.27410v3',
    'paper.pdf': 'https://arxiv.org/pdf/2607.27410v3',
    'validation.html': 'https://arxiv.org/html/2607.27411v2',
    'desi-publications.html': 'https://data.desi.lbl.gov/doc/papers/dr2/',
    'desi-y3-directory.html': 'https://data.desi.lbl.gov/public/papers/y3/',
    'bao-tree.json': f'https://api.github.com/repos/CobayaSampler/bao_data/git/trees/{TREE_COMMIT}?recursive=1',
    'zenodo-arxiv-search.json': 'https://zenodo.org/api/records?'+urllib.parse.urlencode({'q':'2607.27410','size':10}),
    'current-lya-mean.txt': f'https://raw.githubusercontent.com/CobayaSampler/bao_data/{TREE_COMMIT}/desi_bao_dr2/desi_gaussian_bao_Lya_GCcomb_mean.txt',
    'current-lya-cov.txt': f'https://raw.githubusercontent.com/CobayaSampler/bao_data/{TREE_COMMIT}/desi_bao_dr2/desi_gaussian_bao_Lya_GCcomb_cov.txt',
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class Equations(HTMLParser):
    def __init__(self):
        super().__init__(); self.values = []
    def handle_starttag(self, tag, attributes):
        attributes = dict(attributes)
        if tag == 'math' and 'alttext' in attributes:
            self.values.append(attributes['alttext'])


def acquire(cache, output):
    assert not output.exists(), 'Keep the old acquisition receipt; choose a fresh output.'
    cache = cache.resolve()
    assert cache.is_relative_to(ROOT/'.work') and not cache.exists()
    cache.mkdir(parents=True)
    def download(item):
        name, url = item
        request = urllib.request.Request(url, headers={'User-Agent':'cosmology-revalidation/1.0'})
        with urllib.request.urlopen(request, timeout=60) as response:
            content = response.read()
            final_url = response.url
        assert content
        path = cache/name
        path.write_bytes(content)
        return name, {'url':url, 'resolved_url':final_url, 'path':str(path.relative_to(ROOT)),
                      'sha256':sha(path), 'bytes':len(content)}
    with ThreadPoolExecutor(max_workers=4) as pool:
        resources = dict(pool.map(download, SOURCES.items()))
    equations = Equations(); equations.feed((cache/'paper.html').read_text())
    selected = [v for v in equations.values if all(k in v for k in ['39.32','0.33','8.600','0.066','0.225'])]
    assert len(selected) == 1, 'Published Eq26 did not match the independently transcribed summary.'
    tree = json.loads((cache/'bao-tree.json').read_text())
    assert tree['sha'] == TREE_COMMIT and not tree['truncated']
    paths = [r['path'] for r in tree['tree'] if r['type'] == 'blob']
    matches = [p for p in paths if any(s in p.lower() for s in ['lya','fullshape','full_shape','2026'])]
    search = json.loads((cache/'zenodo-arxiv-search.json').read_text())
    hits = [{'id':r['id'], 'title':r['metadata']['title'], 'url':r['links']['self_html']}
            for r in search['hits']['hits']]
    release = ROOT/'.work/unified-cosmology/external-probes/packages/data/bao_data/desi_bao_dr2'
    originals = {'current-lya-mean.txt':release/'desi_gaussian_bao_Lya_GCcomb_mean.txt',
                 'current-lya-cov.txt':release/'desi_gaussian_bao_Lya_GCcomb_cov.txt'}
    old_checks = {name:sha(cache/name) == sha(path) for name,path in originals.items()}
    assert all(old_checks.values()), 'Current Cobaya Ly-alpha files differ from the frozen BAO release; review them.'
    record = {'status':'versioned_sources_acquired_printed_Gaussian_is_explicit_approximation',
              'retrieved_utc':datetime.now(timezone.utc).isoformat(), 'resources':resources,
              'printed_equation_26':selected[0], 'mean_sigma_rho':[39.32,8.600,.33,.066,.225],
              'current_bao_tree_commit':TREE_COMMIT, 'tree_paths_examined':len(paths),
              'relevant_bao_tree_paths':matches, 'current_lya_files_equal_frozen_BAO':old_checks,
              'zenodo_arxiv_query_hits':hits,
              'availability':'No author machine likelihood was identified in these checked interfaces. The primary paper promises later DR2/figure releases; this is not proof of absence everywhere.',
              'limits':'Use the printed Gaussian only as an explicitly approximate replacement. No released non-Gaussian tail or raw-forest reconstruction is claimed. Existing Ly-alpha BAO is overlapping information.',
              'input_sha256':{str(path.relative_to(ROOT)):sha(path) for path in originals.values()},
              'source_sha256':{str(Path(__file__).resolve().relative_to(ROOT)):sha(__file__)}}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    return record


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a = p.parse_args(); result = acquire(a.cache,a.output)
    print(json.dumps({'status':result['status'],'resources':len(result['resources']),
                      'printed_equation_verified':True,'current_bao_tree_commit':TREE_COMMIT}))
