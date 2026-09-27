"""No new CAMB/GPU calls: source density algebra, fixed spectra and sealed caches."""
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['CLIPY_NOJAX'] = '1'
import copy
import gc
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from scipy.special import logsumexp

import joint_lensing_bridge as bridge
from joint_lensing_bridge import ROOT, HERE, DESIGN, GATES, digest, relative, identity, sealed_read, sealed_write


def raw_spectra(path):
    with np.load(path) as archive:
        ell = archive['ell']; dls = archive['spectra']
    factor = np.zeros(len(ell)); factor[1:] = 2*np.pi/(ell[1:]*(ell[1:]+1))
    spectra = {k: value*factor for k, value in zip(['tt', 'ee', 'bb', 'te'], dls[:4])}
    spectra['ell'] = ell
    spectra['pp'] = np.zeros(len(ell))
    spectra['pp'][1:] = dls[4, 1:]*2*np.pi/(ell[1:]*(ell[1:]+1))**2
    return spectra


def main():
    design = json.loads(DESIGN.read_text())
    rng = np.random.default_rng(272919)
    configurations = 0
    for backend in [{}, {'fast_lensing': True}, {'gpu': True, 'fast_lensing': True}]:
        for model in ['lcdm', 'cpl']:
            for evolution in ['none', 'linear', 'smooth01', 'smooth03']:
                for calibration in ['official_planck', 'paper_literal']:
                    settings = dict(backend, model=model, evolution=evolution, calibration=calibration, sample='dovekie')
                    info = bridge.native_configuration(settings)
                    before = bridge.canonical(info)
                    targets = bridge.target_descriptions(info, design, {'test_asset': 'synthetic'})
                    assert bridge.canonical(info) == before
                    for target in targets.values():
                        after = target['configuration']
                        assert after['params'] == before['params']
                        assert after['theory'] == before['theory']
                        assert after['likelihood']['released_sn'] == before['likelihood']['released_sn']
                        assert not any(k in after['likelihood'] for k in design['old_factors'])
                        assert after['params']['A_fg']['prior'] == {'min': 0., 'max': 2.}
                        assert set(after['likelihood']) == (set(before['likelihood'])-set(design['old_factors'])) | {'released_joint_lensing'}
                        configurations += 1
    density_errors = []; constant_errors = []
    for n in range(2, 34):
        matrix = rng.normal(size=(n, n)); covariance = matrix@matrix.T+np.eye(n)
        hartlap = rng.uniform(.8, .99); precision = np.linalg.inv(covariance)*hartlap
        independent = -.5*(np.linalg.slogdet(covariance)[1]-n*np.log(hartlap)+n*np.log(2*np.pi))
        constant_errors.append(abs(bridge.gaussian_constant(precision)-independent))
        parent = {'exact_loglikes': {'oldA': float(rng.normal()), 'oldS': float(rng.normal())}, 'log_weight': float(rng.normal())}
        new = {'kernel_loglike': float(rng.normal()), 'normalized_gaussian_constant': float(rng.normal())}
        constants = {'oldA': 12., 'oldS': 34.}
        result = bridge.bridge_density(parent, new, constants, ['oldA', 'oldS'])
        direct = parent['log_weight']+new['kernel_loglike']-parent['exact_loglikes']['oldA']-parent['exact_loglikes']['oldS']
        density_errors.append(abs(result['target_logweight']-direct))
        same = dict(new, kernel_loglike=sum(parent['exact_loglikes'].values()))
        assert bridge.bridge_density(parent, same, constants, ['oldA', 'oldS'])['target_logweight'] == parent['log_weight']
        # Source closure must include all likelihood factors and priors.
        source = dict(parent, exact_logpost=sum(parent['exact_loglikes'].values())-2.)
        replay = SimpleNamespace(loglikes=list(source['exact_loglikes'].values()), logpriors=[-2.], logpost=source['exact_logpost'])
        assert bridge.source_comparison(source, replay, source['exact_loglikes'], 1e-7)['passed']
        replay.logpriors = [-1.99]
        assert not bridge.source_comparison(source, replay, source['exact_loglikes'], 1e-7)['passed']
    assert max(density_errors) < 1e-13 and max(constant_errors) < 1e-11
    n = 8000; groups = np.repeat(np.arange(4), n//4); x = rng.normal(size=n)
    parents = [{'point': {'x': float(v), 'A_fg': float(a)}, 'derived': {'q0': float(v/2)}, 'log_weight': float(.1*v)}
               for v, a in zip(x, rng.uniform(0, 2, n))]
    rows = [{'status': 'finite_bridge', 'density': {'target_logweight': r['log_weight'], 'fixed_gaussian_normalization_offset': 12.3}} for r in parents]
    gates = json.loads(GATES.read_text())['overlap_gates']
    identical = bridge.summarize_variant(parents, rows, groups, gates)
    assert identical['status'] == 'qualified_conditional_bridge'
    assert 'A_fg' not in identical['posterior']
    assert 'A_fg' in identical['weight_diagnostics_including_auxiliary']['weighted_chain_stability']
    tail = copy.deepcopy(rows); tail[0]['density']['target_logweight'] += 1000
    failure = bridge.summarize_variant(parents, tail, groups, gates)
    assert failure['status'] == 'insufficient_bridge_overlap_or_stability' and failure['posterior'] is None
    assert 'raw_weight_ESS' in failure['failed_gates']
    missing = copy.deepcopy(rows); missing[0] = {'status': 'failed_bridge', 'error': 'synthetic retained failure'}
    assert bridge.summarize_variant(parents, missing, groups, gates)['status'] == 'failed_bridge_evaluation'
    weights = np.exp(x-logsumexp(x)); shifted = np.exp((x+12.3)-logsumexp(x+12.3))
    assert np.max(abs(weights-shifted)) < 1e-15
    checks = {}
    test_plan = {'source_sha256': {}, 'parent_inputs': {}, 'joint_asset_sha256': {},
                 'frozen_target': {'source_sha256': {}, 'versions': {}, 'identity': 'registered'},
                 'settings': {'model': 'cpl', 'evolution': 'none', 'sample': 'dovekie',
                              'calibration': 'official_planck', 'surrogate': 'mock-only'}}
    with patch.object(bridge, 'factory', return_value=(lambda **kwargs: {}, lambda *args: {'identity': 'changed'})):
        try:
            bridge.verify_bindings(test_plan)
            raise RuntimeError('Changed target assets accepted')
        except AssertionError:
            checks['actual_target_identity_rechecked_before_cache_reuse'] = True
    with tempfile.TemporaryDirectory(dir=ROOT/'.work/unified-cosmology') as directory:
        cache = Path(directory); bad_summary = cache/'unqualified.json'
        bad_summary.write_text(json.dumps({'status': 'unqualified'}))
        with patch.object(bridge, 'factory', side_effect=AssertionError('factory must not execute')):
            try:
                bridge.prepare(cache/'no-chain', bad_summary, cache/'must-not-exist')
                raise RuntimeError('Unqualified parent accepted')
            except AssertionError:
                assert not (cache/'must-not-exist').exists()
                checks['unqualified_parent_rejected_before_factory_cache'] = True
        source = cache/'source.json'; source.write_text('{}')
        spectra = cache/'spectra.npz'; np.savez(spectra, tt=np.ones(5))
        plan = {'identity': 'synthetic', 'points': [{'x': 1.}], 'source_native_records': [relative(source)]}
        row = {'identity': 'synthetic', 'index': 0, 'point': {'x': 1.}, 'parent_native_record_sha256': digest(source),
               'spectra_path': relative(spectra), 'spectra_sha256': digest(spectra)}
        sealed_write(cache/'00000.json', row); bridge.verify_cache(cache, plan)
        np.savez(spectra, tt=np.ones(5)*2)
        try:
            bridge.verify_cache(cache, plan); raise RuntimeError('Changed spectra accepted')
        except AssertionError: checks['spectral_numeric_tamper_rejected'] = True
        payload = json.loads((cache/'00000.json').read_text()); payload['point']['x'] = 2.
        (cache/'00000.json').write_text(json.dumps(payload))
        try:
            sealed_read(cache/'00000.json'); raise RuntimeError('Changed payload accepted')
        except AssertionError: checks['record_numeric_tamper_rejected'] = True
        (cache/'00000.json').unlink(); sealed_write(cache/'00000.attempt.json', {'status': 'started'})
        try:
            bridge.verify_cache(cache, plan); raise AssertionError('Incomplete attempt accepted')
        except RuntimeError: checks['incomplete_native_attempt_blocks_retry'] = True
    # Actual stored native spectra are used as data, not recalculated. Include
    # fixed fiducial plus deliberately non-cosmological stress perturbations.
    inputs = {}; closure = []
    index = [0, 1, 2, 3, 4, 32, 64, 96]
    spectra_cases = []
    for i in index:
        path = ROOT/f'.work/unified-cosmology/external-probes/spectral-training/train/{i:04d}.npz'
        inputs[relative(path)] = digest(path)
        spectra_cases.append((str(i), raw_spectra(path)))
    last = spectra_cases[-1][1]
    for i in range(8):
        stress = {k: v.copy() for k, v in last.items()}
        for key in ['tt', 'ee', 'bb', 'te', 'pp']:
            stress[key] *= 1+rng.normal(0, .03)+rng.normal(0, .03)*np.sin(stress['ell']/rng.uniform(50, 2000))
        spectra_cases.append(('synthetic-'+str(i), stress))
    for variant in design['variants']:
        response = bridge.JointResponse(variant)
        for label, spectra in spectra_cases:
            result = response.evaluate(spectra, native_check=True)
            assert result['native_loglike_absolute_difference'] < 1e-8
            closure.append({'variant': variant, 'spectrum': label,
                            'loglike_absolute_difference': result['native_loglike_absolute_difference'],
                            'bandpower_max_absolute_difference': result['native_bandpower_max_absolute_difference']})
        del response; gc.collect()
    sources = [Path(__file__), HERE/'joint_lensing_bridge.py', DESIGN, GATES, HERE/'measurement_summary.py',
               HERE/'spectral_correction.py',
               HERE/'native_posterior_precision.py', HERE/'luminosity_sensitivity.py', HERE/'exact_correction.py',
               HERE.parent/'external_probes/fast_lensing.py', bridge.AUDIT]
    report = {'status': 'passed_fixed_spectrum_and_synthetic_validation', 'new_native_CAMB_calls': 0, 'GPU_calls': 0,
              'target_configuration_cancellation_checks': configurations,
              'density_ratio_cases': len(density_errors), 'maximum_ratio_error': max(density_errors),
              'maximum_effective_covariance_normalization_error': max(constant_errors),
              'identical_target_weight_fixture_points': n, 'identical_target_qualified_synthetic': True,
              'heavy_tail_and_native_failure_withhold_posterior': True, 'cache_and_qualification_checks': checks,
              'actual_stored_spectra': len(index), 'synthetic_spectrum_perturbations': 8,
              'native_generic_vs_compressed': closure,
              'source_sha256': {relative(p): digest(p) for p in sources}, 'fixed_spectrum_sha256': inputs,
              'scope': 'Arithmetic and interface checks only. No actual posterior bridge, no current measurement replacement, and no evidence claim.'}
    output = ROOT/'studies/unified_cosmology/results/inference/joint-lensing-bridge-validation.json'
    output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: report[k] for k in ['status', 'new_native_CAMB_calls', 'target_configuration_cancellation_checks', 'maximum_ratio_error']}))


if __name__ == '__main__':
    main()
