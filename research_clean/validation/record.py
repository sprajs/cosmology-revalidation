#!/usr/bin/env python3
"""Verify the frozen edition and new run hashes; write an additive audit manifest."""
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
    parser.add_argument('--prefix', default='revalidation-20260926')
    args = parser.parse_args()
    edition = read(ROOT/'provenance/edition.json')
    failures = [p for p, h in edition['files_sha256'].items() if not (ROOT/p).is_file() or sha(ROOT/p) != h]
    assert not failures, f'Frozen edition changed: {failures}'
    inputs = read(ROOT/'provenance/inputs.json')['files']
    assert all((ROOT/x['path']).stat().st_size == x['bytes'] and sha(ROOT/x['path']) == x['sha256'] for x in inputs)
    runs = {}
    comparisons = {}
    for path in sorted((ROOT/'results').glob(args.prefix+'-*/run.json')):
        run = read(path)
        assert run['status'] == 'completed', path
        for p, h in run['code_sha256'].items():
            assert sha(ROOT/p) == h, p
        for p, h in run['inputs_sha256'].items():
            assert sha(ROOT/p) == h, p
        for p, h in run['outputs_sha256'].items():
            assert sha(path.parent/p) == h, p
        summary = read(path.parent/'summary.json')
        assert summary.get('valid_for_posterior_summary', True), path
        runs[path.parent.name] = {'workflow': run['workflow'], 'status': run['status'],
            'record_sha256': sha(path), 'outputs_verified': len(run['outputs_sha256']),
            'posterior_gate': summary.get('valid_for_posterior_summary'),
            'tracked_record_formats': ['.json', '.csv'],
            'local_array_outputs': [p for p in run['outputs_sha256'] if Path(p).suffix in {'.npz', '.npy'}]}
        if path.parent.name.startswith(args.prefix+'-v1-'):
            name = path.parent.name.removeprefix(args.prefix+'-v1-')
            identical = summary == read(ROOT/'results'/name/'summary.json')
            comparisons[name] = {'summary_exactly_equal_to_frozen_edition': identical}
    assert len(comparisons) == 19 and len(runs) == 29, (len(comparisons), len(runs))
    # A later campaign may legitimately change random draws or platform rounding.
    # Record comparison truthfully rather than hiding a mismatch.
    seeds = {}
    for first, second in [('v1-cosmology','lcdm-second-seed'),('joint','joint-second-seed'),('qbins','qbins-second-seed'),('template','template-second-seed')]:
        a = read(ROOT/'results'/f'{args.prefix}-{first}'/'summary.json')
        b = read(ROOT/'results'/f'{args.prefix}-{second}'/'summary.json')
        gaps = {p: (b['posterior'][p]['mean']-v['mean'])/v['sd'] for p,v in a['posterior'].items()}
        seeds[first] = {'parameter_mean_shifts_in_first_run_posterior_sd': gaps,
            'q_mean_shift': b['q0']['mean']-a['q0']['mean'],
            'interpretation': 'Sampler sensitivity check; posterior SD is not Monte Carlo SE or a proof of convergence.'}
    reports = ROOT/'validation/reports'
    for name, key in [('core','audit_code_sha256'),('des','source_sha256'),('raisin','audit_code_sha256'),('hst','script_sha256')]:
        report = read(reports/(name+'.json'))
        assert report[key] == sha(ROOT/'validation'/(name+'_checks.py')), name
        if name in {'des','raisin'}:
            assert report['status'] == 'passed'
        if name == 'hst':
            assert report['checks']['count'] == report['checks']['passed']
    comparison = {'default_runs': comparisons, 'independent_seed_comparisons': seeds}
    (reports/'replay-comparison.json').write_text(json.dumps(comparison, indent=2)+'\n')
    files = [p for p in (ROOT/'validation').rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name != 'manifest.json' and p.suffix != '.log']
    record = {'schema': 'clean-research-revalidation-v1', 'created_utc': datetime.now(timezone.utc).isoformat(),
        'scope': 'Complete replay of the frozen clean edition, all supplied configurations and five independent seeds; targeted independent numerical checks. Not a complete survey reduction or unified cosmology likelihood.',
        'frozen_edition_sha256': sha(ROOT/'provenance/edition.json'),
        'frozen_edition_files_verified_unchanged': len(edition['files_sha256']),
        'input_files_verified': len(inputs), 'runs': runs,
        'validation_files_sha256': {str(p.relative_to(ROOT)): sha(p) for p in sorted(files)},
        'original_workflow_code_changed': False,
        'publication': 'JSON/CSV results and validation code are versioned. Bulk inputs and new array outputs remain local; their SHA-256 values are recorded in each run.json. No claim that a Git clone alone contains all inputs.'}
    (ROOT/'validation/manifest.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps({'runs_verified': len(runs), 'default_summaries_equal': sum(x['summary_exactly_equal_to_frozen_edition'] for x in comparisons.values()), 'frozen_files_unchanged': len(edition['files_sha256']), 'inputs_verified': len(inputs)}))


if __name__ == '__main__':
    main()
