"""Qualified-parent joint-lensing sensitivity with immutable native spectra.

Default prepares a plan only. --execute-native explicitly launches work later.
No live target, active sampler, or existing correction files are modified.
"""
from __future__ import annotations
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
import argparse
from concurrent.futures import ThreadPoolExecutor
import copy
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from scipy.special import logsumexp

from native_posterior_precision import digest, relative, identity, sealed_read, sealed_write, factory
from target_identity import canonical
from luminosity_sensitivity import weight_diagnostics

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DESIGN = HERE/'joint-lensing-bridge-design.json'
GATES = HERE/'luminosity-sensitivity-design.json'
AUDIT = ROOT/'studies/unified_cosmology/results/inference/joint-lensing-availability.json'
EXTERNAL = HERE.parent/'external_probes'
RELEASE = ROOT/'.work/unified-cosmology/joint-lensing-availability/release'


def gaussian_constant(precision):
    precision = np.asarray(precision)
    sign, logdet = np.linalg.slogdet(precision)
    assert sign > 0
    return float(.5*logdet-.5*len(precision)*np.log(2*np.pi))


def native_configuration(settings):
    configuration, _ = factory(settings)
    kwargs = {k: settings[k] for k in ['model', 'evolution', 'sample', 'calibration']}
    info = configuration(**kwargs)
    assert set(info['theory']) == {'camb'}, 'Only unchanged native CAMB is allowed.'
    design = json.loads(DESIGN.read_text())
    assert all(k in info['likelihood'] for k in design['old_factors'])
    assert info['params']['A_fg']['prior'] == design['auxiliary_parameter']['normalized_prior']
    old = info['likelihood'][design['old_factors'][0]]
    assert old['variant'] == 'actplanck_baseline' and not old['lens_only']
    assert old['apply_hartlap'] and old['nsims_planck'] == 400
    assert info['likelihood']['SPT2023_lensing']['clear_internal_priors']
    return info


def target_descriptions(native, design, joint_hashes):
    original = canonical(native)
    descriptions = {}
    for variant in design['variants']:
        config = copy.deepcopy(original)
        for key in design['old_factors']:
            del config['likelihood'][key]
        config['likelihood']['released_joint_lensing'] = {
            'generic_native_module': 'act_dr6_spt_lenslike', 'variant': variant,
            'lens_only': False, 'like_corrections': True, 'apply_hartlap': True,
            'nsims_act': 792, 'nsims_planck': 400, 'trim_lmax': 2998,
            'act_calib': False, 'indep': False}
        # The auxiliary prior is intentionally retained; this is a density
        # description, not a Cobaya config with an unused parameter.
        restored = copy.deepcopy(config)
        del restored['likelihood']['released_joint_lensing']
        restored['likelihood'].update({k: original['likelihood'][k] for k in design['old_factors']})
        assert restored == original
        description = {'configuration': config, 'auxiliary_parameter': design['auxiliary_parameter'],
                       'joint_asset_sha256': joint_hashes, 'interpretation': design['variant_interpretation'][variant]}
        description['identity'] = identity(description)
        descriptions[variant] = description
    return descriptions


def verify_bindings(plan):
    for mapping in [plan['source_sha256'], plan['parent_inputs'], plan['joint_asset_sha256'], plan['frozen_target']['source_sha256']]:
        assert all(digest(ROOT/p) == h for p, h in mapping.items()), 'Source, parent or likelihood input changed.'
    for package, version in plan['frozen_target']['versions'].items():
        assert importlib.metadata.version(package) == version, 'Frozen environment changed.'
    from late_geometry import sample_path
    configuration, identify_target = factory(plan['settings'])
    kwargs = {k: plan['settings'][k] for k in ['model', 'evolution', 'sample', 'calibration']}
    proposal = configuration(**kwargs, surrogate=Path(plan['settings']['surrogate']))
    assert identify_target(proposal, sample_path(plan['settings']['sample']), Path(plan['settings']['surrogate'])) == plan['frozen_target'], 'Parent target assets changed.'


def prepare(folder, summary_path, cache):
    from measurement_summary import summarize_run
    from late_geometry import sample_path
    # A status string alone is insufficient: the consumer recomputes all gates.
    parent = summarize_run(folder, summary_path)
    summary = json.loads(summary_path.read_text())
    selection_path = ROOT/summary['selection_path']
    selection = json.loads(selection_path.read_text())
    settings = selection['settings']
    frozen = json.loads((folder/'run-0.json').read_text())['target_identity']
    configuration, identify_target = factory(settings)
    kwargs = {k: settings[k] for k in ['model', 'evolution', 'sample', 'calibration']}
    proposal = configuration(**kwargs, surrogate=Path(settings['surrogate']))
    assert identify_target(proposal, sample_path(settings['sample']), Path(settings['surrogate'])) == frozen
    info = native_configuration(settings)
    design = json.loads(DESIGN.read_text())
    audit = json.loads(AUDIT.read_text())
    assert audit['status'] == 'released_likelihood_evaluation_reproduced_no_inference'
    joint = dict(audit['input_sha256'])
    joint.update(audit['source_sha256'])
    from spectral_correction import sidecar_path, read_sidecar
    captured = []
    for i in range(len(selection['points'])):
        native_path = selection_path.parent/f'{i:05d}.json'
        path = sidecar_path(native_path)
        if path.exists():
            row = read_sidecar(native_path)
            parent['input_sha256'][relative(path)] = digest(path)
            if row['status'] == 'captured_from_parent_native_evaluation':
                parent['input_sha256'][row['spectra_path']] = row['spectra_sha256']
                captured.append(relative(path))
                continue
        captured.append(None)
    sources = [Path(__file__), DESIGN, GATES, AUDIT, HERE/'measurement_summary.py',
               HERE/'exact_correction.py', HERE/'native_posterior_precision.py',
               HERE/'spectral_correction.py',
               HERE/'luminosity_sensitivity.py', HERE/'target_identity.py',
               EXTERNAL/'fast_lensing.py']
    plan = {'settings': settings, 'design': design, 'groups': selection['groups'],
            'locations': selection['locations'], 'points': selection['points'],
            'source_native_records': [relative(selection_path.parent/f'{i:05d}.json') for i in range(len(selection['points']))],
            'captured_spectral_sidecars': captured,
            'frozen_target': frozen, 'source_native_configuration': canonical(info),
            'parent_inputs': parent['input_sha256'],
            'source_sha256': {relative(p): digest(p) for p in sources}, 'joint_asset_sha256': joint,
            'target_descriptions': target_descriptions(info, design, joint),
            'correction_summary_path': relative(summary_path), 'correction_summary_sha256': digest(summary_path)}
    plan['identity'] = identity(plan)
    verify_bindings(plan)
    cache = cache.resolve()
    assert cache.is_relative_to(ROOT/'.work'), 'Generated spectra must remain ignored.'
    cache.mkdir(parents=True, exist_ok=True)
    path = cache/'plan.json'
    if path.exists():
        old = sealed_read(path)
        assert {k: v for k, v in old.items() if k != 'payload_sha256'} == plan
    else:
        sealed_write(path, plan)
    verify_cache(cache, plan)
    return sealed_read(path)


def verify_cache(cache, plan):
    manifest = cache/'native-records.json'
    if manifest.exists():
        old = sealed_read(manifest)
        assert old['identity'] == plan['identity']
        for path, expected in old['file_sha256'].items():
            assert digest(ROOT/path) == expected, 'Completed cache file changed.'
    for i, point in enumerate(plan['points']):
        path = cache/f'{i:05d}.json'
        if path.exists():
            row = sealed_read(path)
            assert row['identity'] == plan['identity'] and row['index'] == i and row['point'] == point
            assert row['parent_native_record_sha256'] == digest(ROOT/plan['source_native_records'][i])
            if 'spectra_path' in row:
                assert digest(ROOT/row['spectra_path']) == row['spectra_sha256'], 'Spectral payload changed.'
        elif (cache/f'{i:05d}.attempt.json').exists():
            raise RuntimeError(f'Incomplete native attempt {i}; its evaluation count is unknown. No silent retry.')


def source_comparison(parent, result, names, tolerance):
    values = dict(zip(names, map(float, result.loglikes)))
    assert set(values) == set(parent['exact_loglikes'])
    differences = {k: values[k]-parent['exact_loglikes'][k] for k in values}
    prior_difference = float(sum(result.logpriors)-(parent['exact_logpost']-sum(parent['exact_loglikes'].values())))
    post_difference = float(result.logpost-parent['exact_logpost'])
    assert np.isfinite([*values.values(), *differences.values(), prior_difference, post_difference]).all()
    return {'component_loglikes': values, 'component_closure': differences,
            'prior_closure': prior_difference, 'logpost_closure': post_difference,
            'passed': max(abs(x) for x in [*differences.values(), prior_difference, post_difference]) <= tolerance}


def source_normalization(model, design):
    act = model.likelihood[design['old_factors'][0]]
    spt = model.likelihood[design['old_factors'][1]].candl_like
    assert not spt.add_logdet and len(spt.priors) == 0
    chol = np.asarray(spt.covariance_chol_dec)
    spt_const = float(-np.log(np.diag(chol)).sum()-.5*len(chol)*np.log(2*np.pi))
    return {design['old_factors'][0]: gaussian_constant(act.data['cinv']), design['old_factors'][1]: spt_const}


def reuse_captured_spectrum(cache, plan, i):
    """Use the parent's actual native spectrum without another model evaluation."""
    from spectral_correction import read_sidecar
    from exact_correction import verify_record
    path = cache/f'{i:05d}.json'
    assert not path.exists() and not (cache/f'{i:05d}.attempt.json').exists()
    source = ROOT/plan['source_native_records'][i]
    parent = json.loads(source.read_text()); verify_record(parent)
    assert parent['point'] == plan['points'][i] and parent['status'] == 'finite'
    assert digest(source) == plan['parent_inputs'][relative(source)]
    capture = read_sidecar(source)
    capture_path = ROOT/plan['captured_spectral_sidecars'][i]
    assert digest(capture_path) == plan['parent_inputs'][relative(capture_path)]
    assert capture['status'] == 'captured_from_parent_native_evaluation'
    assert capture['additional_native_evaluations'] == 0
    assert capture['spectra_sha256'] == plan['parent_inputs'][capture['spectra_path']]
    row = {'identity': plan['identity'], 'index': i, 'point': plan['points'][i],
           'parent_native_record_sha256': digest(source), 'native_evaluations': 0,
           'status': 'finite_source_closure', 'seconds': 0.,
           'source_comparison': {'passed': True,
               'method': 'Same native evaluation as the qualified parent; hash-bound spectral capture, not a second density replay.',
               'capture_sidecar_path': relative(capture_path), 'capture_sidecar_sha256': digest(capture_path)},
           'source_normalized_gaussian_constants': capture['source_normalized_gaussian_constants'],
           'original_exact_proposal_logweight': parent['log_weight'],
           'spectra_path': capture['spectra_path'], 'spectra_sha256': capture['spectra_sha256'],
           'spectral_units': capture['spectral_units']}
    sealed_write(path, row)
    return row


def native_worker(plan_path, indices):
    from cobaya.model import get_model
    from exact_correction import verify_record
    plan_path = Path(plan_path).resolve(); cache = plan_path.parent
    plan = sealed_read(plan_path); verify_bindings(plan)
    info = native_configuration(plan['settings'])
    assert canonical(info) == plan['source_native_configuration']
    with get_model(info) as model:
        assert 'camb' in model.theory and set(model.theory) <= {'camb', 'camb.transfers'}
        constants = source_normalization(model, plan['design'])
        for i in indices:
            path = cache/f'{i:05d}.json'; attempt = cache/f'{i:05d}.attempt.json'
            assert not path.exists() and not attempt.exists(), 'No repeat of attempted native work.'
            parent_path = ROOT/plan['source_native_records'][i]
            parent = json.loads(parent_path.read_text()); verify_record(parent)
            point = plan['points'][i]; assert parent['point'] == point and parent['status'] == 'finite'
            assert digest(parent_path) == plan['parent_inputs'][relative(parent_path)]
            row = {'identity': plan['identity'], 'index': i, 'point': point,
                   'parent_native_record_sha256': digest(parent_path), 'native_evaluations': 0}
            sealed_write(attempt, dict(row, status='about_to_evaluate'))
            started = time.monotonic()
            try:
                row['native_evaluations'] = 1
                result = model.logposterior(point)
                closure = source_comparison(parent, result, model.likelihood,
                                            plan['design']['source_closure_absolute_tolerance'])
                row['source_comparison'] = closure
                assert closure['passed'], 'Stored native source density did not reproduce.'
                raw = model.provider.get_Cl(ell_factor=False, units='FIRASmuK2')
                assert set(['ell', 'tt', 'ee', 'te', 'bb', 'pp']) <= set(raw)
                spectra = {k: np.asarray(raw[k]) for k in ['ell', 'tt', 'ee', 'te', 'bb', 'pp']}
                assert all(np.isfinite(v).all() for v in spectra.values())
                assert spectra['ell'][0] == 0 and len(spectra['pp']) >= 3102
                spectrum_path = cache/f'{i:05d}-spectra.npz'
                assert not spectrum_path.exists()
                np.savez_compressed(spectrum_path, **spectra)
                row.update(status='finite_source_closure', source_comparison=closure,
                           source_normalized_gaussian_constants=constants,
                           original_exact_proposal_logweight=parent['log_weight'],
                           spectra_path=relative(spectrum_path), spectra_sha256=digest(spectrum_path),
                           spectral_units='raw C_l: TT/EE/TE/BB FIRASmuK2, pp dimensionless potential')
            except Exception as error:
                row.update(status='failed_native_evaluation', error=repr(error))
            row['seconds'] = time.monotonic()-started
            sealed_write(path, row)
            print(json.dumps({'native_index': i, 'status': row['status']}), flush=True)


class JointResponse:
    """Existing reviewed exact matrix compression, with direct SPT windows."""
    def __init__(self, variant):
        sys.path.insert(0, str(RELEASE)); sys.path.insert(0, str(EXTERNAL))
        import act_dr6_spt_lenslike as native
        assert Path(native.__file__).resolve().is_relative_to(RELEASE.resolve()), 'Different joint release imported.'
        from fast_lensing import BinnedResponse
        self.native = native
        self.data = native.load_data(variant, lens_only=False, like_corrections=True,
                                     nsims_act=792, nsims_planck=400)
        self.compressed = BinnedResponse(self.data)
        self.normalization = gaussian_constant(self.data['cinv'])
        self.variant = variant

    def evaluate(self, spectra, native_check=False):
        n = self.native; ell = spectra['ell']; kk = n.pp_to_kk(spectra['pp'], ell)
        short = n.standardize(ell, kk, 2998)
        cmb = {k: n.standardize(ell, spectra[k], 2998) for k in ['tt', 'ee', 'bb', 'te']}
        prediction = np.r_[self.compressed.predict(short, cmb),
                           self.data['binmat_spt'] @ n.standardize(ell, kk, 3100)]
        r = self.data['data_binned_clkk']-prediction
        ll = float(-.5*r@self.data['cinv']@r)
        result = {'kernel_loglike': ll, 'normalized_gaussian_constant': self.normalization,
                  'normalized_gaussian_loglike': ll+self.normalization, 'bandpowers': prediction.tolist()}
        if native_check:
            original, bp = n.generic_lnlike(self.data, ell, kk, ell, spectra['tt'], spectra['ee'], spectra['te'], spectra['bb'], return_theory=True)
            result.update(native_loglike_absolute_difference=float(abs(ll-original)),
                          native_bandpower_max_absolute_difference=float(abs(prediction-bp).max()))
        return result


def bridge_density(parent, new, constants, old_factors):
    old = sum(parent['exact_loglikes'][key] for key in old_factors)
    delta = float(new['kernel_loglike']-old)
    offset = float(new['normalized_gaussian_constant']-sum(constants.values()))
    return {'lens_kernel_logratio': delta, 'fixed_gaussian_normalization_offset': offset,
            'normalized_density_logratio': delta+offset,
            'target_logweight': float(parent['log_weight']+delta)}


def summarize_variant(parent, rows, groups, gates):
    if any(r['status'] != 'finite_bridge' for r in rows):
        return {'status': 'failed_bridge_evaluation', 'posterior': None,
                'failures': [{'index': i, 'status': r['status'], 'error': r.get('error')} for i, r in enumerate(rows) if r['status'] != 'finite_bridge']}
    lw = np.array([r['density']['target_logweight'] for r in rows])
    names = sorted(set(parent[0]['point']) | set(parent[0]['derived']))
    values = {k: np.array([dict(r['point'], **r['derived'])[k] for r in parent]) for k in names}
    checks = weight_diagnostics(lw, values, np.asarray(groups), gates)
    failures = list(checks['failed_gates'])
    if len(rows) < gates['minimum_exact_points']:
        failures.append('minimum_exact_points')
    weights = np.exp(lw-logsumexp(lw))
    physical = [k for k in names if k != 'A_fg']
    matrix = np.column_stack([values[k] for k in physical]); centered = matrix-weights@matrix
    offsets = np.array([r['density']['fixed_gaussian_normalization_offset'] for r in rows])
    assert np.ptp(offsets) < 1e-9, 'A supposedly fixed Gaussian normalization changed.'
    return {'status': 'qualified_conditional_bridge' if not failures else 'insufficient_bridge_overlap_or_stability',
            'qualified_under_declared_numerical_gates': not failures, 'failed_gates': failures,
            'points': len(rows), 'weight_diagnostics_including_auxiliary': checks,
            'posterior': {k: checks['weighted_summaries_for_diagnostics'][k] for k in physical} if not failures else None,
            'weighted_covariance': {'parameter_order': physical, 'matrix': ((centered*weights[:, None]).T@centered).tolist(),
                                    'qualified': not failures, 'normalization': 'Untrimmed normalized weights; no Bessel correction.'},
            'fixed_normalized_gaussian_offset': float(offsets[0]),
            'auxiliary_A_fg': 'Normalized uniform[0,2] density retained; not reported as a physical measurement.',
            'evidence_or_global_mode_coverage_claimed': False}


def complete(cache, plan, output):
    verify_bindings(plan); verify_cache(cache, plan)
    parents = []; native = []; files = {}
    for i in range(len(plan['points'])):
        path = cache/f'{i:05d}.json'
        assert path.exists(), 'Incomplete native cohort.'
        row = sealed_read(path); native.append(row)
        parents.append(json.loads((ROOT/plan['source_native_records'][i]).read_text()))
        files[relative(path)] = digest(path)
        if 'spectra_path' in row:
            files[row['spectra_path']] = row['spectra_sha256']
    sealed_write(cache/'native-records.json', {'identity': plan['identity'], 'file_sha256': files})
    results = {}
    for variant in plan['design']['variants']:
        response = JointResponse(variant)
        variant_dir = cache/variant; variant_dir.mkdir(exist_ok=True)
        manifest_path = variant_dir/'records.json'
        if manifest_path.exists():
            previous = sealed_read(manifest_path)
            for path, value in previous['file_sha256'].items():
                assert digest(ROOT/path) == value, 'Completed bridge record changed.'
        rows = []; variant_hashes = {}
        for i, (parent, source) in enumerate(zip(parents, native)):
            path = variant_dir/f'{i:05d}.json'
            if path.exists():
                row = sealed_read(path)
                assert row['index'] == i
                assert row['plan_identity'] == plan['identity'] and row['target_identity'] == plan['target_descriptions'][variant]['identity']
                assert row['source_record_sha256'] == files[relative(cache/f'{i:05d}.json')]
            else:
                try:
                    assert source['status'] == 'finite_source_closure', 'Native source failed.'
                    with np.load(ROOT/source['spectra_path'], allow_pickle=False) as spectra:
                        new = response.evaluate(spectra)
                    density = bridge_density(parent, new, source['source_normalized_gaussian_constants'], plan['design']['old_factors'])
                    assert np.isfinite(list(density.values())).all()
                    row = {'status': 'finite_bridge', 'new_lensing': new, 'density': density}
                except Exception as error:
                    row = {'status': 'failed_bridge', 'error': repr(error)}
                row.update(index=i, plan_identity=plan['identity'], target_identity=plan['target_descriptions'][variant]['identity'],
                           source_record_sha256=files[relative(cache/f'{i:05d}.json')])
                sealed_write(path, row); row = sealed_read(path)
            rows.append(row); variant_hashes[relative(path)] = digest(path)
        sealed_write(manifest_path, {'identity': plan['target_descriptions'][variant]['identity'], 'file_sha256': variant_hashes})
        results[variant] = summarize_variant(parents, rows, plan['groups'], json.loads(GATES.read_text())['overlap_gates'])
        results[variant].update(target=plan['target_descriptions'][variant], record_manifest_path=relative(manifest_path), record_manifest_sha256=digest(manifest_path))
        del response
        import gc; gc.collect()
    verify_bindings(plan); verify_cache(cache, plan)
    result = {'status': 'separate_lensing_sensitivities', 'plan_identity': plan['identity'],
              'plan_path': relative(cache/'plan.json'), 'plan_sha256': digest(cache/'plan.json'),
              'native_manifest_path': relative(cache/'native-records.json'), 'native_manifest_sha256': digest(cache/'native-records.json'),
              'native_evaluations': sum(r['native_evaluations'] for r in native), 'variants': results,
              'limitations': plan['design']['limitations'], 'active_target_changed': False}
    output.parent.mkdir(parents=True, exist_ok=True)
    sealed_write(output, result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--chain-folder', type=Path)
    p.add_argument('--correction-summary', type=Path)
    p.add_argument('--cache', type=Path)
    p.add_argument('--output', type=Path)
    p.add_argument('--execute-native', action='store_true')
    p.add_argument('--summarize-only', action='store_true')
    p.add_argument('--workers', type=int, default=1)
    p.add_argument('--worker-plan', type=Path)
    p.add_argument('--indices')
    a = p.parse_args()
    if a.worker_plan:
        native_worker(a.worker_plan, [int(x) for x in a.indices.split(',')]); return
    assert a.chain_folder and a.correction_summary and a.cache
    assert 1 <= a.workers <= 4
    plan = prepare(a.chain_folder.resolve(), a.correction_summary.resolve(), a.cache.resolve())
    if not a.execute_native and not a.summarize_only:
        print(json.dumps({'status': 'qualified_plan_only', 'points': len(plan['points']), 'variants': plan['design']['variants'], 'native_calls': 0})); return
    cache = a.cache.resolve()
    if a.execute_native:
        missing = [i for i in range(len(plan['points'])) if not (cache/f'{i:05d}.json').exists()]
        for i in missing:
            if plan['captured_spectral_sidecars'][i] is not None:
                reuse_captured_spectrum(cache, plan, i)
        missing = [i for i in missing if not (cache/f'{i:05d}.json').exists()]
        def run(indices):
            if not indices: return
            log = cache/('worker-'+str(indices[0])+'.log')
            assert not log.exists(), 'Preserve prior worker log; do not silently restart its attempt.'
            with log.open('w') as stream:
                subprocess.run([sys.executable, str(Path(__file__).resolve()), '--worker-plan', str(cache/'plan.json'),
                                '--indices', ','.join(map(str, indices))], stdout=stream, stderr=subprocess.STDOUT, check=True)
        with ThreadPoolExecutor(max_workers=a.workers) as pool:
            list(pool.map(run, [missing[j::a.workers] for j in range(a.workers)]))
    assert a.output
    result = complete(cache, plan, a.output.resolve())
    print(json.dumps({'status': result['status'], 'variant_status': {k: v['status'] for k, v in result['variants'].items()}}))


if __name__ == '__main__':
    main()
