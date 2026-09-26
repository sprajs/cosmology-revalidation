#!/usr/bin/env python3
"""Separate, hashed binary/path amendment for the frozen CSP 42 baseline.

Only nominal and nominal_copy are executable here. The changed-filter arm
requires its own release and runner after the changed-pilot gate.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'runs/research_2026_09_26/csp_native_filter_response/cohort-execution'
BASE_SOURCE = ROOT / 'scripts/research_2026_09_26/csp_filter_42_runner.py'
BASE_PROTOCOL = OUT / 'execution-protocol.json'
RELEASE = ROOT / 'runs/research_2026_09_26/csp_native_filter_response/root-state-review/full42-baseline-release.json'
NATIVE = ROOT / 'runs/research_2026_09_26/raisin_profile_solver_review/csp-filter-interpretation/native-design/instrumentation/SNANA-v11_04k-output/bin/snlc_fit.exe'
ALIAS = Path('/tmp/snana-csp-audit-8c5f0d')
AMENDMENT = OUT / 'instrumented-baseline-amendment.json'
BASE_PROTOCOL_SHA = 'd1ed7f0a076c277538866125e5816c5b95a9cecac2f106f764d429faf0bfd604'
RELEASE_SHA = '9f23661e1c84eef18c4aa1e63545f561cae3c1514b3126f781a86307b0f358ac'
NATIVE_SHA = '716d5c235468329e878bd56e678cfa55b2a316f333c36ec411e9a320c1770938'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load_base():
    spec = importlib.util.spec_from_file_location('csp_frozen_runner', BASE_SOURCE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def gate():
    assert sha(BASE_PROTOCOL) == BASE_PROTOCOL_SHA
    assert sha(RELEASE) == RELEASE_SHA
    assert sha(NATIVE) == NATIVE_SHA
    release = json.loads(RELEASE.read_text())
    assert release['pass'] is True
    assert release['approved_for_full42_baseline_only'] is True
    assert release['approved_for_full42_changed_response'] is False
    assert release['native_binary_sha256'] == NATIVE_SHA
    assert ALIAS.is_symlink() and ALIAS.resolve() == NATIVE.parent.parent.resolve()
    assert (ALIAS / 'bin/snlc_fit.exe').read_bytes() == NATIVE.read_bytes()
    return release

def freeze(changed_ledger, changed_seconds):
    release = gate()
    assert not AMENDMENT.exists(), 'Do not overwrite an execution amendment'
    changed_ledger = changed_ledger.resolve()
    assert changed_ledger and changed_ledger.is_file()
    assert changed_seconds >= 0
    prior = float(release['completed_native_seconds']) + changed_seconds
    assert prior < 600
    resource = json.loads(changed_ledger.read_text())
    assert abs(float(resource['all_active_native_seconds']) - prior) < 1e-8
    constituents = {}
    for row in resource['ledger']:
        path = ROOT / row['path']
        constituents[row['path']] = sha(path)
    assert abs(sum(float(r['seconds']) for r in resource['ledger']) - prior) < 1e-8
    amend = {
        'status': 'Frozen before any full42 native execution',
        'scope': ['nominal42', 'nominal_copy42'],
        'not_released': 'known_Jdw_to_j42',
        'base_protocol_sha256': BASE_PROTOCOL_SHA,
        'base_source_sha256': sha(BASE_SOURCE),
        'wrapper_source_sha256': sha(Path(__file__)),
        'root_baseline_release_sha256': RELEASE_SHA,
        'instrumented_binary_sha256': NATIVE_SHA,
        'native_alias': str(ALIAS),
        'native_alias_target': str(ALIAS.resolve()),
        'release_completed_native_seconds': release['completed_native_seconds'],
        'changed_pilot_ledger': str(changed_ledger.relative_to(ROOT)),
        'changed_pilot_ledger_sha256': sha(changed_ledger),
        'resource_constituents_sha256': constituents,
        'changed_pilot_active_seconds': changed_seconds,
        'prior_native_active_seconds': prior,
        'total_native_active_seconds_cap': 600,
        'environment': {'SNANA_DIR': str(ALIAS), 'SNDATA_ROOT': str(ROOT/'phase2/official/inputs/SNDATA_ROOT'),
                        'LD_LIBRARY_PATH': str(ROOT/'phase2/official/build/sysroot/usr/lib'),
                        'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'},
        'reason': 'Native v11_04k source-gated instrumented binary emits every object/iteration state and support row. The prior d1ed... input/source protocol remains immutable; no objective or NML parameter change.'
    }
    AMENDMENT.write_text(json.dumps(amend,indent=2)+'\n')
    print('amendment_sha256', sha(AMENDMENT))

def run(arm):
    assert arm in ('nominal','nominal_copy')
    release = gate()
    amend = json.loads(AMENDMENT.read_text())
    assert amend['status'] == 'Frozen before any full42 native execution'
    assert amend['base_protocol_sha256'] == BASE_PROTOCOL_SHA
    assert amend['base_source_sha256'] == sha(BASE_SOURCE)
    assert amend['wrapper_source_sha256'] == sha(Path(__file__))
    assert amend['root_baseline_release_sha256'] == RELEASE_SHA
    assert amend['instrumented_binary_sha256'] == NATIVE_SHA
    assert sha(ROOT / amend['changed_pilot_ledger']) == amend['changed_pilot_ledger_sha256']
    for path, digest in amend['resource_constituents_sha256'].items():
        assert sha(ROOT / path) == digest, path
    prior = float(amend['prior_native_active_seconds'])
    assert prior >= float(release['completed_native_seconds'])
    mod = load_base()
    mod.BIN = NATIVE
    mod.ALREADY_USED = prior
    def instrumented_env():
        env = os.environ.copy()
        env.pop('CSP_POSTINIT_DSHIFT', None)
        env.update(amend['environment'])
        return env
    mod.env_native = instrumented_env
    mod.run_arm(arm, RELEASE, RELEASE_SHA)
    record = {'arm':arm,'amendment_sha256':sha(AMENDMENT),
              'native_binary_sha256':NATIVE_SHA,'root_release_sha256':RELEASE_SHA,
              'native_execution_sha256':sha(OUT/'full/fits'/arm/'execution.json')}
    (OUT/'full/fits'/arm/'instrumented-execution.json').write_text(json.dumps(record,indent=2)+'\n')

if __name__ == '__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('action',choices=['freeze','run'])
    ap.add_argument('--arm',choices=['nominal','nominal_copy'])
    ap.add_argument('--changed-pilot-ledger',type=Path)
    ap.add_argument('--changed-pilot-seconds',type=float)
    a=ap.parse_args()
    if a.action=='freeze':
        assert a.changed_pilot_seconds is not None
        freeze(a.changed_pilot_ledger,a.changed_pilot_seconds)
    else:
        assert a.arm is not None
        run(a.arm)
