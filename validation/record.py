#!/usr/bin/env python3
"""Verify inputs, recorded run outputs and audit code, with explicit history."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read(path):
    return json.loads(path.read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--name', help='Campaign below results; omit for the published reference results')
    args = parser.parse_args()
    if args.name and (Path(args.name).name != args.name or args.name in {'.', '..'}):
        parser.error('name must be one directory name')
    results = ROOT/'results'
    if args.name:
        results /= args.name
    layout = read(ROOT/'provenance/layout.json')
    historical_code = layout['historical_workflow_code_sha256']
    inputs = read(ROOT/'provenance/inputs.json')['files']
    assert all((ROOT/x['path']).stat().st_size == x['bytes'] and sha(ROOT/x['path']) == x['sha256'] for x in inputs)
    runs, comparisons = {}, {}
    for role in ['baseline', 'alternatives', 'robustness']:
        for path in sorted((results/role).glob('*/run.json')):
            run = read(path)
            assert run['status'] == 'completed', path
            current_lock = sha(ROOT/'uv.lock')
            assert run['lock_sha256'] in {current_lock, layout['historical_lock_sha256']}, path
            previous_code = []
            for p, h in run['code_sha256'].items():
                current = sha(ROOT/p)
                if current != h:
                    assert historical_code.get(p) == h, (p, 'unknown historical code')
                    migration = layout['workflow_path_changes'].get(p)
                    assert migration and migration['current_sha256'] == current, (p, 'unrecorded current change')
                    previous_code.append(p)
            for p, h in run['inputs_sha256'].items():
                assert sha(ROOT/p) == h, p
            for p, h in run['outputs_sha256'].items():
                assert sha(path.parent/p) == h, p
            summary = read(path.parent/'summary.json')
            assert summary.get('valid_for_posterior_summary', True), path
            relative = str(path.parent.relative_to(ROOT))
            runs[relative] = {'workflow': run['workflow'], 'status': run['status'], 'run_sha256': sha(path),
                'outputs_verified': len(run['outputs_sha256']), 'historical_path_code': previous_code,
                'lock_identity': 'current' if run['lock_sha256'] == current_lock else 'historical',
                'posterior_gate': summary.get('valid_for_posterior_summary'),
                'local_array_outputs': [p for p in run['outputs_sha256'] if Path(p).suffix in {'.npz', '.npy'}]}
            if role == 'baseline':
                expected = layout['reference_summary_sha256'][path.parent.name]
                comparisons[path.parent.name] = {'summary_matches_pre_reorganization': sha(path.parent/'summary.json') == expected}
    assert len(comparisons) == 19 and len(runs) == 29, (len(comparisons), len(runs))
    reports = ROOT/'validation/reports'
    for name, key in [('core','audit_code_sha256'),('des','source_sha256'),('raisin','audit_code_sha256'),('hst','script_sha256')]:
        report = read(reports/(name+'.json'))
        assert report[key] == sha(ROOT/'validation'/(name+'_checks.py')), name
        if name in {'des','raisin'}:
            assert report['status'] == 'passed'
        if name == 'hst':
            assert report['checks']['count'] == report['checks']['passed']
    record = {'schema': 'cosmology-validation-v2', 'created_utc': datetime.now(timezone.utc).isoformat(),
        'result_directory': str(results.relative_to(ROOT)), 'input_files_verified': len(inputs),
        'runs': runs, 'baseline_comparisons': comparisons,
        'scope': 'Verified frozen inputs and all outputs of 29 recorded runs, plus targeted independent audits. Historical run records retain their original code hashes; path-only changes are explicit in provenance/layout.json.',
        'audit_sha256': {str(p.relative_to(ROOT)): sha(p) for p in sorted(reports.glob('*.json'))},
        'audit_code_sha256': {str(p.relative_to(ROOT)): sha(p) for p in sorted((ROOT/'validation').glob('*.py'))},
        'publication': 'JSON/CSV scientific records are versioned; bulk data and sampler arrays require the input bundle or regeneration.'}
    out = ROOT/'validation/manifest.json' if not args.name else results/'validation.json'
    out.write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps({'runs_verified': len(runs), 'baseline_summaries_unchanged': sum(x['summary_matches_pre_reorganization'] for x in comparisons.values()), 'inputs_verified': len(inputs), 'manifest': str(out.relative_to(ROOT))}))


if __name__ == '__main__':
    main()
