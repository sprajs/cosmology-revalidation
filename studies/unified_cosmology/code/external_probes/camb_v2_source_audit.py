"""Replay a source-only CAMB release audit; never import or execute CAMB."""
import argparse
import hashlib
import json
from pathlib import Path
import urllib.request

ROOT = Path(__file__).resolve().parents[4]
REPORT = ROOT / 'studies/unified_cosmology/results/inference/camb-v2-source-audit.json'


def verify(work, acquire=False):
    report = json.loads(REPORT.read_text())
    for row in report['sources']:
        path = work / row['path']
        if not path.exists() and acquire:
            request = urllib.request.Request(row['url'], headers={'User-Agent': 'cosmology-revalidation-source-audit'})
            data = urllib.request.urlopen(request, timeout=60).read()
            assert hashlib.sha256(data).hexdigest() == row['sha256'], row['url']
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open('xb') as stream:
                stream.write(data)
        data = path.read_bytes()
        assert len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256'], str(path)
    old = (work / 'v166/fortran/cmbmain.f90').read_text().lower()
    sub = (work / 'v166/fortran/subroutines.f90').read_text().lower()
    compact = ''.join(old.split())
    assert 'allocate(state%transfer_times(0:state%num_transfer_redshifts+1))' in compact
    assert 'callstate%timeofzarr(state%transfer_times(1:),' in compact
    assert 'callspline_def(state%transfer_times,scaling,state%num_transfer_redshifts,ddscaling)' in compact
    assert 'real(sp_acc),intent(in)::x(n),y(n)' in ''.join(sub.split())
    fix = json.loads((work / 'head-commit.json').read_text())
    assert fix['sha'] == report['source_commits']['head_spline_fix']
    assert len(fix['files']) == 1 and fix['files'][0]['filename'] == 'fortran/cmbmain.f90'
    assert '+            call cubic_spline_second_derivs(State%transfer_times(1:State%num_transfer_redshifts)' in fix['files'][0]['patch']
    return {'status': 'passed_source_hashes_and_indexing_contract', 'sources': len(report['sources']), 'physical_calls': 0}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work', type=Path, default=ROOT / '.work/unified-cosmology/camb-v2-source-audit')
    parser.add_argument('--acquire', action='store_true', help='Fetch missing exact source bytes; changed upstream bytes refuse.')
    args = parser.parse_args()
    print(json.dumps(verify(args.work, args.acquire)))
