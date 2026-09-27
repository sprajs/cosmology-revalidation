"""Bounded, read-only search for the missing released covariance inputs.

Enumerates filenames in a supplied checkout/archive; reads only scoped small
text files and bounded relevant tar headers. This is not a claim to have searched
all public repositories, private production directories, or large tar interiors.
No covariance or scientific likelihood is changed.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tarfile
import tempfile
import urllib.request

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
WORK = ROOT / '.work/unified-cosmology/calibration-input-search'
SOURCE_FILE = Path(__file__).with_name('calibration-input-sources.json')
NEEDLE = 'v6_9_duplicate_cid.cov'
NAME_PATTERN = re.compile(r'v6[_\.]?9|duplicate.?cid|ceph.*cov|cov.*ceph', re.I)
CONTENT_PATTERN = re.compile(r'v6[_\.]?9_duplicate_cid|ceph[^\n]{0,60}cov|cov[^\n]{0,60}ceph', re.I)
METADATA_PATTERN = re.compile(r'(acquisition|manifest|provenance|registry|sources?|tree|checksum)', re.I)
MAX_TEXT_BYTES = 2 * 1024**2
MAX_TAR_BYTES = 64 * 1024**2
PUBLIC_QUERIES = {'github-code-exact': NEEDLE,
                  'github-code-author-duplicate': 'duplicate_cid user:djbrout',
                  'github-code-author-ceph-cov': 'cepheid cov user:djbrout',
                  'github-code-scolnic-ceph': 'CEPH user:dscolnic',
                  'github-code-cepheid-cov': 'cepheid covariance repo:PantheonPlusSH0ES/DataRelease'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def label(path, archive):
    # rg is called with absolute roots. Preserve symlink entry names and avoid
    # resolving every path through the filesystem in the filename-only pass.
    path = Path(path).absolute()
    for base, prefix in [(ROOT, '$REPO'), (archive, '$ARCHIVE')]:
        if path.is_relative_to(base):
            return prefix + '/' + str(path.relative_to(base))
    raise ValueError('Unexpected input outside specified search roots')


def source_inventory(acquire):
    records = json.loads(SOURCE_FILE.read_text())
    for name, entry in records.items():
        path = WORK / name
        if not path.exists():
            if not acquire:
                raise FileNotFoundError(f'{path}; use --acquire')
            request = urllib.request.Request(entry['url'], headers={'User-Agent': 'cosmology-revalidation-input-audit'})
            with urllib.request.urlopen(request, timeout=60) as response:
                body = response.read()
            if hashlib.sha256(body).hexdigest() != entry['sha256']:
                raise ValueError(f'Upstream bytes changed: {name}')
            path.write_bytes(body)
        assert path.stat().st_size == entry['bytes'] and sha(path) == entry['sha256']
    trees = []
    for name in records:
        if not name.endswith('-tree.json'):
            continue
        tree = json.loads((WORK / name).read_text())
        assert tree['truncated'] is False
        trees.append({'source': name, 'commit': tree['sha'], 'entries': len(tree['tree']),
                      'exact_filename_matches': [x['path'] for x in tree['tree'] if Path(x['path']).name == NEEDLE],
                      'candidate_names': [x['path'] for x in tree['tree'] if NAME_PATTERN.search(x['path'])]})
    # Check raw files against their recorded Git tree blobs independently of SHA256.
    blobs = []
    for key in ['pippin', 'snana', 'website']:
        tree = json.loads((WORK / (key + '-tree.json')).read_text())
        by_path = {x['path']: x for x in tree['tree']}
        for name in records:
            if not name.startswith(key + '-') or name.endswith('.json') or name.endswith('-historical.py'):
                continue
            upstream = name[len(key) + 1:].replace('__', '/')
            if upstream not in by_path:
                continue
            body = (WORK / name).read_bytes()
            digest = hashlib.sha1(b'blob ' + str(len(body)).encode() + b'\0' + body).hexdigest()
            assert digest == by_path[upstream]['sha']
            blobs.append({'file': name, 'git_blob_verified': digest})
    return records, trees, blobs


def historical_assignment_probe():
    """Execute only the isolated public external-covariance reader on fake rows."""
    source = WORK / 'snana-historical.py'
    tree = ast.parse(source.read_text())
    function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'get_cov_from_covfile')
    namespace = {'np': np, 'pd': pd}
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(source), 'exec'), namespace)
    data = pd.DataFrame({'CIDstr': ['A_1', 'B_2', 'C_3']})
    rows = pd.DataFrame({'CID1': ['A', 'B', 'A'], 'IDSURVEY1': [1, 2, 1],
                         'CID2': ['B', 'A', 'B'], 'IDSURVEY2': [2, 1, 2],
                         'MU_COV': [.004, .004, .007]})
    with tempfile.TemporaryDirectory(dir=WORK) as temp:
        path = Path(temp) / 'synthetic.csv'
        rows.to_csv(path, index=False)
        cov, _ = namespace['get_cov_from_covfile'](data, str(path), 1.)
    expected = np.array([[0., .007, 0.], [.004, 0., 0.], [0., 0., 0.]])
    assert np.array_equal(cov, expected)
    return {'status': 'passed_isolated_historical_reader_probe', 'source_sha256': sha(source),
            'synthetic_output': cov.tolist(), 'interpretation':
            'At scale1, the inspected June2022 reader assigns each ordered CID+IDSURVEY pair directly; repeated ordered rows overwrite. It does not itself symmetrize, double, or infer a Cepheid host covariance. This source is not certified as the actual release production checkout.'}


def public_search(refresh):
    """Optional live index query; preserve each response, including limits/errors."""
    records = []
    for name, query in PUBLIC_QUERIES.items():
        path = WORK / (name + '.json')
        command = ['gh', 'search', 'code', query, '--limit', '100', '--json', 'repository,path,url']
        if refresh:
            response = subprocess.run(command, capture_output=True, text=True)
            if response.returncode:
                records.append({'query': query, 'status': 'query_failed', 'stderr': response.stderr})
                continue
            # Save the prior snapshot if a later live index gives different bytes.
            body = response.stdout.encode()
            if path.exists() and path.read_bytes() != body:
                path.with_name(name + '-' + sha(path)[:16] + '.json').write_bytes(path.read_bytes())
            path.write_bytes(body)
        if not path.exists():
            records.append({'query': query, 'status': 'not_executed_use_--search-public'})
            continue
        hits = json.loads(path.read_text())
        records.append({'query': query, 'status': 'index_response_not_exhaustive_repository_search',
                        'path': str(path.relative_to(ROOT)), 'sha256': sha(path),
                        'returned_count': len(hits), 'limit': 100, 'may_be_capped': len(hits) >= 100,
                        'matches': hits})
    return records


def local_inventory(archive):
    command = ['rg', '--files', '--hidden', '--no-ignore', str(ROOT), str(archive)]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    files = sorted(set(Path(p) for p in result.stdout.splitlines()))
    # Preserve the filename-only inventory, with portable root labels.
    inventory = '\n'.join(label(p, archive) for p in files) + '\n'
    inventory_path = WORK / 'local-filenames.txt'
    inventory_path.write_text(inventory)
    exact = [label(p, archive) for p in files if p.name == NEEDLE]
    candidates = [label(p, archive) for p in files if NAME_PATTERN.search(p.name)]
    excluded = ('/.git/', '/.venv/', '/.modern-venv/', '/site-packages/', '/node_modules/')
    scanned, hits, skipped_text, tar_records = [], [], [], []
    for path in files:
        text_path = str(path)
        if any(s in text_path for s in excluded) or path.is_relative_to(WORK):
            continue
        # Only original author snapshots and acquisition/provenance metadata;
        # avoid generated point records and unrelated cache contents.
        author = '/sources/repos/PantheonPlusSH0ES__DataRelease' in text_path
        related = re.search(r'pantheon|shoes|covariance', text_path, re.I) is not None
        metadata = METADATA_PATTERN.search(path.name) and path.suffix.lower() in {'.json', '.yaml', '.yml', '.md'}
        selected = (author and path.suffix.lower() in {'.py', '.md', '.yml', '.yaml'}) or (related and metadata)
        if not selected:
            continue
        nbytes = path.stat().st_size
        if nbytes > MAX_TEXT_BYTES:
            skipped_text.append({'path': label(path, archive), 'bytes': nbytes})
            continue
        body = path.read_bytes()
        try:
            text = body.decode('utf-8')
        except UnicodeDecodeError:
            continue
        digest = hashlib.sha256(body).hexdigest()
        scanned.append({'path': label(path, archive), 'bytes': nbytes, 'sha256': digest})
        matches = [{'line': j, 'text': line[:700]} for j, line in enumerate(text.splitlines(), 1) if CONTENT_PATTERN.search(line)]
        if matches:
            hits.append({'path': label(path, archive), 'sha256': digest, 'matches': matches})
    # Inspect only relevant compact bundles. Full SNDATA_ROOT tarballs are
    # explicitly deferred, while their extracted filesystem trees were listed.
    for path in files:
        name = str(path)
        if not re.search(r'\.(tar|tar\.gz|tgz)$', name, re.I):
            continue
        selected = (path.name.startswith('SNDATA_ROOT_') or path.name == 'NewLCs.tar' or
                    (re.search(r'pantheon|shoes|ceph|pplus', path.name, re.I) is not None))
        if not selected:
            continue
        nbytes = path.stat().st_size
        entry = {'path': label(path, archive), 'bytes': nbytes}
        if nbytes > MAX_TAR_BYTES:
            entry.update(status='not_opened_above_declared_compressed_byte_budget', member_names_checked=False)
        else:
            with tarfile.open(path, 'r:*') as tf:
                names = tf.getnames()
            entry.update(status='member_headers_checked_no_extraction', member_names_checked=True,
                         sha256=sha(path), members=len(names),
                         exact_filename_matches=[n for n in names if Path(n).name == NEEDLE],
                         candidate_names=[n for n in names if NAME_PATTERN.search(Path(n).name)])
        tar_records.append(entry)
    detail = WORK / 'local-content-search.json'
    write(detail, {'scanned': scanned, 'hits': hits, 'skipped_text': skipped_text, 'tar_headers': tar_records})
    return {'filename_count': len(files), 'inventory_path': str(inventory_path.relative_to(ROOT)),
            'inventory_sha256': sha(inventory_path), 'exact_filename_matches': exact,
            'candidate_filenames': candidates, 'small_text_files_scanned': len(scanned),
            'small_text_bytes_scanned': sum(x['bytes'] for x in scanned),
            'text_hits': hits, 'skipped_text_count': len(skipped_text), 'tar_headers': tar_records,
            'detail_path': str(detail.relative_to(ROOT)), 'detail_sha256': sha(detail)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--archive', type=Path, default=Path.home() / '.local/share/cosmology-revalidation/archive')
    parser.add_argument('--acquire', action='store_true')
    parser.add_argument('--search-public', action='store_true', help='Refresh bounded public GitHub code searches using gh')
    parser.add_argument('--output', type=Path, default=ROOT / 'studies/unified_cosmology/results/distance_ladder/calibration-input-search.json')
    args = parser.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)
    records, trees, blobs = source_inventory(args.acquire)
    search = public_search(args.search_public)
    local = local_inventory(args.archive.resolve())
    probe = historical_assignment_probe()
    result = {'status': 'bounded_search_complete_construction_input_unrecovered',
              'exact_input': NEEDLE, 'local_search': local, 'public_trees': trees,
              'upstream_sources': records, 'verified_git_blobs': blobs,
              'public_code_search': search,
              'historical_reader_probe': probe,
              'source_interpretation': [
                  'PPLUS.yml lists DUP_SIGINT scale1 as an external file and includes it in NOSYS. Its CALIBRATORS list is not a numerical host-distance covariance.',
                  'The inspected May2022 Pippin wrapper forwards EXTRA_COVS and CALIBRATORS. June2022 SNANA suppresses VPEC/ZSHIFT rows for calibrators and assigns external MU_COV cells; no Cepheid host-to-measurement embedding operation is present in these two generic files.',
                  'The current SNANA builder has later external-covariance changes and cannot certify the production version used in2022.',
                  'The retrieved AS759 notebook explicitly teaches a simplified two-rung ladder without SNe; it is not the Pantheon+ production embedding recipe.'],
              'scope_limits': ['Filename absence is restricted to these enumerated roots and pinned complete public trees.',
                               'The content search covers scoped small author/metadata files, not all files or all public repositories.',
                               'Large compressed archive interiors and historical Git objects were not exhaustively searched; extracted SNDATA_ROOT filesystem filenames are included.',
                               'The exact external duplicate-scatter table, Cepheid host covariance version, and embedding program remain unknown. No inference that DUP_SIGINT is itself the Cepheid covariance is made.',
                               'The near2x sibling pattern remains unresolved; no covariance correction or substitution is justified by this search.'],
              'covariance_changed': False, 'cosmology_fits': 0,
              'source_sha256': {str(p.relative_to(ROOT)): sha(p) for p in [Path(__file__), SOURCE_FILE]}}
    write(args.output, result)
    print(json.dumps({'status': result['status'], 'local_files': local['filename_count'],
                      'text_files': local['small_text_files_scanned'], 'exact_matches': local['exact_filename_matches'],
                      'public_tree_matches': [x['exact_filename_matches'] for x in trees]}))


if __name__ == '__main__':
    main()
