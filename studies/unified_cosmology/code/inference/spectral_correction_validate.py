"""Validate native spectral sidecars without changing density-correction code."""
import os
for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[key] = '1'
os.environ['CLIPY_NOJAX'] = '1'
import argparse
import copy
import json
import multiprocessing as mp
from pathlib import Path
import tempfile
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

import exact_correction as core
import spectral_correction as capture
import joint_lensing_bridge as bridge

ROOT = capture.ROOT
HERE = capture.HERE


class FakeModel:
    def __init__(self, shift=0.):
        self.shift = shift; self.calls = 0; self.point = None
        self.provider = self
        self.parameterization = SimpleNamespace(derived_params=lambda: ['q0'])
        self.likelihood = {'act_dr6_lenslike.ACTDR6LensLike': SimpleNamespace(data={'cinv': np.eye(2)}),
                           'SPT2023_lensing': SimpleNamespace(candl_like=SimpleNamespace(
                               add_logdet=False, priors=[], covariance_chol_dec=np.eye(3)))}

    def logposterior(self, point):
        self.calls += 1; self.point = point.copy()
        return SimpleNamespace(logpost=-point['x']**2+self.shift,
            loglikes=[-.4*point['x']**2+self.shift, -.6*point['x']**2], derived=[-point['x']/10])

    def get_Cl(self, ell_factor, units):
        assert ell_factor is False and units == 'FIRASmuK2'
        ell = np.arange(3201)
        return dict(ell=ell, **{k: (self.point['x']+self.shift)*np.ones(len(ell))
                               for k in ['tt', 'ee', 'te', 'bb', 'pp']})


def initialize_fake():
    core._exact = FakeModel(); core._proposal = FakeModel(.25)


def synthetic_checks():
    checks = {}
    with tempfile.TemporaryDirectory(dir=ROOT/'.work/unified-cosmology') as directory:
        folder = Path(directory)
        tasks = [(i, {'x': float(i+1)}, str(folder), 'synthetic', -(i+1.)**2+.25) for i in range(4)]
        with mp.get_context('spawn').Pool(2, initializer=initialize_fake) as pool:
            rows = list(pool.map(capture.evaluate_with_spectra, tasks))
        for task, row in zip(tasks, rows):
            assert row['status'] == 'finite' and row['log_weight'] == -.25
            side = capture.read_sidecar(folder/f'{task[0]:05d}.json')
            assert side['status'] == 'captured_from_parent_native_evaluation'
            assert side['additional_native_evaluations'] == 0
            with np.load(ROOT/side['spectra_path']) as spectra:
                assert np.array_equal(spectra['tt'], np.ones(3201)*task[1]['x'])
        checks['spawn_workers_preserve_exact_not_proposal_spectra'] = True
        initialize_fake()
        before = core._exact.calls + core._proposal.calls
        assert capture.evaluate_with_spectra(tasks[0]) == rows[0]
        assert core._exact.calls + core._proposal.calls == before
        checks['replay_uses_zero_model_evaluations'] = True
        native = folder/'00000.json'; sidepath = capture.sidecar_path(native)
        side = capture.read_sidecar(native)
        parent_bytes = native.read_bytes(); spectrum = ROOT/side['spectra_path']
        spectrum_bytes = spectrum.read_bytes()
        spectrum.write_bytes(spectrum_bytes+b'changed')
        try:
            capture.read_sidecar(native); raise RuntimeError('Changed spectrum accepted')
        except AssertionError: checks['spectrum_tamper_rejected'] = True
        spectrum.write_bytes(spectrum_bytes)
        bad = json.loads(parent_bytes); bad['exact_logpost'] += 1
        bad['payload_sha256'] = core.record_digest(bad)
        native.write_text(json.dumps(bad))
        try:
            capture.read_sidecar(native); raise RuntimeError('Resealed changed parent accepted')
        except AssertionError: checks['resealed_parent_change_rejected'] = True
        native.write_bytes(parent_bytes)
        with patch.object(capture, 'dependencies', return_value={'changed': 'source'}):
            try:
                capture.read_sidecar(native); raise RuntimeError('Changed capture source accepted')
            except AssertionError: checks['capture_source_change_rejected'] = True
            calls = core._exact.calls
            try:
                capture.evaluate_with_spectra(tasks[0]); raise RuntimeError('Changed startup source accepted')
            except AssertionError:
                assert core._exact.calls == calls
                checks['startup_source_change_blocks_before_any_evaluation'] = True
        plan = {'identity': 'synthetic-plan', 'points': [rows[0]['point']],
                'source_native_records': [capture.relative(native)],
                'captured_spectral_sidecars': [capture.relative(sidepath)],
                'parent_inputs': {capture.relative(p): capture.digest(p) for p in [native, sidepath, spectrum]}}
        reuse = folder/'bridge'; reuse.mkdir()
        bridged = bridge.reuse_captured_spectrum(reuse, plan, 0)
        assert bridged['native_evaluations'] == 0 and bridged['original_exact_proposal_logweight'] == -.25
        assert bridged['spectra_sha256'] == side['spectra_sha256']
        bridge.verify_cache(reuse, plan)
        checks['bridge_reuses_same_hash_bound_native_spectrum_without_model'] = True
        # An older scalar record does not license reading an unrelated current state.
        oldtask = (4, {'x': 5.}, str(folder), 'synthetic', -24.75)
        oldrow = core.evaluate(oldtask)
        count = core._exact.calls
        assert capture.evaluate_with_spectra(oldtask) == oldrow
        assert core._exact.calls == count
        assert capture.read_sidecar(folder/'00004.json')['status'] == 'unavailable_from_preexisting_native_record'
        checks['uncaptured_old_record_does_not_reuse_stale_state_or_recompute'] = True
        with patch.dict(os.environ, {'CLIPY_NOJAX': '0'}):
            badtask = (5, {'x': 6.}, str(folder), 'synthetic', -35.75)
            row = capture.evaluate_with_spectra(badtask)
            assert row['status'] == 'finite'
            assert capture.read_sidecar(folder/'00005.json')['status'] == 'capture_failed_parent_density_retained'
        checks['optional_capture_failure_retains_original_density'] = True
    return checks


def native_checks(output_dir):
    from cobaya.model import get_model
    from modern_fast import configuration
    folder = output_dir.resolve(); folder.mkdir(parents=True, exist_ok=False)
    source = ROOT/'studies/unified_cosmology/results/external_probes/author-configuration.json'
    proposal = json.loads(source.read_text())['models']['lcdm']['proposal_only']
    settings = {'model': 'lcdm', 'evolution': 'none', 'sample': 'dovekie',
                'calibration': 'official_planck', 'fast_lensing': True,
                'surrogate': str(ROOT/'.work/unified-cosmology/inference/surrogate/cubic-0509.npz')}
    info = configuration(**{k: settings[k] for k in ['model', 'evolution', 'sample', 'calibration']})
    names = [k for k, v in info['params'].items() if isinstance(v, dict) and 'prior' in v]
    centre = {k: float(proposal['transformed_mean'][k]) for k in names}
    # Fixed perturbations: no outcome-based choice and no posterior interpretation.
    points = [centre, dict(centre, H0=centre['H0']+.15),
              dict(centre, ombh2=centre['ombh2']+.00003), dict(centre, ns=centre['ns']-.001)]
    core.initialize(settings)
    fresh = get_model(info)
    rows = []
    for i, point in enumerate(points):
        predicted = core._proposal.logposterior(point)
        row = capture.evaluate_with_spectra((i, point, str(folder), 'native-spectral-validation', float(predicted.logpost)))
        assert row['status'] == 'finite'
        side = capture.read_sidecar(folder/f'{i:05d}.json')
        assert side['status'] == 'captured_from_parent_native_evaluation'
        replay = fresh.logposterior(point)
        comparison = bridge.source_comparison(row, replay, fresh.likelihood, 1e-7)
        assert comparison['passed']
        expected = fresh.provider.get_Cl(ell_factor=False, units='FIRASmuK2')
        with np.load(ROOT/side['spectra_path']) as spectra:
            differences = {k: float(np.max(np.abs(spectra[k]-expected[k]))) for k in spectra}
            assert all(np.array_equal(spectra[k], expected[k]) for k in spectra)
        rows.append({'index': i, 'point': point, 'density_closure': comparison,
                     'spectral_max_absolute_differences': differences,
                     'sidecar_path': capture.relative(capture.sidecar_path(folder/f'{i:05d}.json')),
                     'sidecar_sha256': capture.digest(capture.sidecar_path(folder/f'{i:05d}.json'))})
        print(json.dumps({'native_capture_validation': i, 'passed': True}), flush=True)
    fresh.close(); core._exact.close(); core._proposal.close()
    return {'native_logposterior_invocations': 8, 'extra_native_calls_due_to_capture': 0,
            'cases': rows, 'input_sha256': {capture.relative(source): capture.digest(source),
                capture.relative(settings['surrogate']): capture.digest(settings['surrogate'])}}


def main():
    p = argparse.ArgumentParser(); p.add_argument('--native-cache', type=Path)
    p.add_argument('--output', type=Path, required=True); a = p.parse_args()
    report = {'status': 'passed_synthetic_only', 'synthetic_checks': synthetic_checks(),
              'core_correction_source_unchanged': capture.digest(HERE/'exact_correction.py')}
    if a.native_cache:
        report['native_validation'] = native_checks(a.native_cache)
        report['status'] = 'passed_native_and_synthetic_spectral_capture'
    report['source_sha256'] = dict(capture.dependencies(), **{capture.relative(Path(__file__)): capture.digest(Path(__file__))})
    report['scope'] = 'Storage and lineage validation, not a posterior result. Native density correction code, selection and gates remain unchanged.'
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': report['status'], 'synthetic_checks': len(report['synthetic_checks'])}))


if __name__ == '__main__':
    main()
