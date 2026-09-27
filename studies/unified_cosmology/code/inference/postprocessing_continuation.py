"""Explicit, finite continuation after the preserved quantile-order failure.

No original record, source, selection, weight or failed receipt is replaced.
Only missing post-processing stages may be requested, once per receipt folder.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
ALLOWED = {'quantile_precision.py', 'expansion_history.py', 'luminosity_history.py',
           'sn_predictive_check.py', 'native_posterior_precision.py',
           'probe_omission.py', 'luminosity_bridge.py'}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_new(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def stopped_plan(folder):
    pp = folder/'finite-pipeline-plan.json'
    sp = folder/'finite-pipeline-status.json'
    plan, status = [json.loads(p.read_text()) for p in [pp, sp]]
    assert status['stage'] == 'stopped_no_retry'
    assert status['plan_sha256'] == digest(pp)
    assert 'quantile_precision.py' in status['error']
    completed = status['completed_results']
    assert [r['script'] for r in completed] == [
        'diagnostics.py', 'spectral_correction.py', 'measurement_summary.py']
    assert completed[-1]['status'] == 'qualified_conditional_measurements'
    assert Path(plan['folder']).resolve() == folder
    for r in completed:
        assert digest(ROOT/r['path']) == r['sha256']
    # A still-running original launcher must retain ownership of its cohort.
    proc = Path('/proc')/str(status['launcher_pid'])
    if proc.exists():
        command = (proc/'cmdline').read_bytes().split(b'\0')
        assert not any(Path(x.decode()).name == 'finite_pipeline.py' for x in command if x)
    return plan, {str(p.relative_to(ROOT)): digest(p) for p in [pp, sp]}, completed


def choose(plan, stages):
    assert stages and len(stages) == len(set(stages)) and set(stages) <= ALLOWED
    rows = []
    for name in stages:
        matches = [dict(r) for r in plan['stages'] if r['script'] == name]
        assert matches, name
        for row in matches:
            assert not row['condition'], 'Conditional native-precision stages need separate review.'
            assert digest(HERE/name) == row['source_sha256']
            if name == 'quantile_precision.py':
                replacement = HERE/'quantile_precision_ordered.py'
                row['preserved_original_script'] = name
                row['preserved_original_source_sha256'] = row['source_sha256']
                row['script'] = replacement.name
                row['source_sha256'] = digest(replacement)
                row['argv'] = [str(replacement) if x == str(HERE/name) else x for x in row['argv']]
            assert not Path(row['output']).exists(), 'Never overwrite a completed output.'
            rows.append(row)
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--folder', type=Path, required=True)
    p.add_argument('--stages', nargs='+', required=True, choices=sorted(ALLOWED))
    p.add_argument('--receipt', type=Path, required=True)
    p.add_argument('--execute', action='store_true')
    a = p.parse_args()
    folder = a.folder.resolve()
    plan, inputs, completed = stopped_plan(folder)
    rows = choose(plan, a.stages)
    correction = ROOT/completed[1]['path']
    from measurement_summary import summarize_run
    qualified = summarize_run(folder, correction)
    record = {'schema': 'explicit-postprocessing-continuation-v1',
              'original_failed_plan_and_status_sha256': inputs,
              'source_sha256': digest(__file__), 'folder': str(folder),
              'target_identity': qualified['target_identity'], 'stages': rows,
              'interpretation': 'Original failure retained; no sampling, correction, or scientific gate is changed.'}
    if not a.execute:
        print(json.dumps(record, indent=2)); return
    a.receipt.mkdir(parents=True, exist_ok=False)
    write_new(a.receipt/'plan.json', record)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               MKL_NUM_THREADS='1', CLIPY_NOJAX='1')
    results = []
    try:
        for index, row in enumerate(rows):
            stopped_plan(folder)
            assert digest(__file__) == record['source_sha256']
            for path, expected in inputs.items():
                assert digest(ROOT/path) == expected
            assert digest(HERE/row['script']) == row['source_sha256']
            assert not Path(row['output']).exists()
            summarize_run(folder, correction)
            write_new(a.receipt/f'{index:02d}-attempt.json', {
                'script': row['script'], 'argv': row['argv'],
                'started_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())})
            with (a.receipt/f'{index:02d}.log').open('x') as log:
                subprocess.run(row['argv'], cwd=ROOT, env=env, stdout=log,
                               stderr=subprocess.STDOUT, check=True)
            output = Path(row['output'])
            value = json.loads(output.read_text())
            result = {'script': row['script'], 'path': str(output.relative_to(ROOT)),
                      'sha256': digest(output), 'status': value['status']}
            write_new(a.receipt/f'{index:02d}-result.json', result)
            results.append(result)
            print(json.dumps(result), flush=True)
        write_new(a.receipt/'completed.json', {'results': results,
                  'status': 'execution_complete_scientific_statuses_retained',
                  'plan_sha256': digest(a.receipt/'plan.json')})
    except BaseException as error:
        write_new(a.receipt/'failure.json', {'error': repr(error), 'completed': results,
                  'plan_sha256': digest(a.receipt/'plan.json'), 'retry_authorized': False})
        raise


if __name__ == '__main__':
    main()
