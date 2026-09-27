"""Three uncached full-native evaluations at two already fixed control points."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time

import numpy as np

from exact_correction import verify_record
from late_geometry import sample_path
from modern_fast import configuration, identify

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DESIGN = HERE/'native-order-design.json'
CONTROLS = HERE/'expansion-history-native-controls.json'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def differences(left, right):
    assert set(left['loglikes']) == set(right['loglikes'])
    assert set(left['derived']) == set(right['derived'])
    return {'logpost_absolute': abs(left['logpost']-right['logpost']),
            'component_absolute': {k: abs(left['loglikes'][k]-right['loglikes'][k]) for k in left['loglikes']},
            'derived_scaled': {k: abs(left['derived'][k]-right['derived'][k])/(1+abs(right['derived'][k])) for k in left['derived']}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    assert all(os.environ.get(k) == '1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','CLIPY_NOJAX'])
    design = json.loads(DESIGN.read_text())
    controls = json.loads(CONTROLS.read_text())
    info = configuration(model='lcdm', evolution='none', sample='dovekie', calibration='official_planck')
    assert set(info['theory']) == {'camb'}
    target = identify(info, sample_path('dovekie'))
    bindings = {relative(p): sha(p) for p in [Path(__file__), DESIGN, CONTROLS]}
    stored = {}
    for index in [0, 1]:
        entry = controls['controls'][index]
        path = ROOT/entry['native_record_path']
        assert sha(path) == entry['native_record_sha256']
        bindings[relative(path)] = sha(path)
        value = json.loads(path.read_text()); verify_record(value)
        assert value['status'] == 'finite' and value['point'] == entry['point']
        stored[index] = value
    plan = {'design': design, 'target': target, 'input_source_sha256': bindings}
    plan_path = args.output.with_suffix('.design.json')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with plan_path.open('x') as stream:
        json.dump(plan, stream, indent=2); stream.write('\n')
    from cobaya.model import get_model
    records = []
    with get_model(info) as model:
        for index in [0, 1, 0]:
            point = stored[index]['point']
            start = time.monotonic()
            value = model.logposterior(point, cached=False)
            assert np.isfinite(value.logpost)
            spectra = model.provider.get_Cl(ell_factor=True)
            row = {'control': index, 'point': point, 'logpost': float(value.logpost),
                   'logpriors': list(map(float, value.logpriors)),
                   'loglikes': dict(zip(model.likelihood, map(float, value.loglikes))),
                   'derived': dict(zip(model.parameterization.derived_params(), map(float, value.derived))),
                   'spectrum_hashes': {k: hashlib.sha256(np.asarray(spectra[k]).tobytes()).hexdigest() for k in ['tt','ee','te','bb','pp']},
                   'seconds': time.monotonic()-start}
            records.append(row)
    checks = {'cold_versus_revisited': differences(records[0], records[2])}
    for i in [0, 1]:
        source = stored[i]
        reference = {'logpost': source['exact_logpost'], 'loglikes': source['exact_loglikes'], 'derived': source['derived']}
        checks['control_'+str(i)+'_versus_stored'] = differences(records[i], reference)
    failed = []
    for name, check in checks.items():
        if max(check['logpost_absolute'], *check['component_absolute'].values()) > design['absolute_logdensity_and_component_tolerance']:
            failed.append(name+':density')
        if max(check['derived_scaled'].values()) > design['relative_derived_tolerance']:
            failed.append(name+':derived')
    same_spectra = records[0]['spectrum_hashes'] == records[2]['spectrum_hashes']
    if not same_spectra:
        failed.append('revisited_spectra_not_bitwise_identical')
    for path, checksum in bindings.items():
        assert sha(ROOT/path) == checksum
    assert identify(info, sample_path('dovekie')) == target
    result = {'status': 'passed_bounded_native_order_check' if not failed else 'failed_bounded_native_order_check',
              'failed_gates': failed, 'checks': checks, 'records': records, 'native_evaluations': len(records),
              'revisited_spectra_bitwise_identical': same_spectra,
              'target_identity': target, 'input_source_sha256': bindings,
              'execution_design_path': relative(plan_path), 'execution_design_sha256': sha(plan_path),
              'posterior_qualification': False, 'limits': design['limits']}
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': result['status'], 'checks': checks, 'failed_gates': failed}))
    if failed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
