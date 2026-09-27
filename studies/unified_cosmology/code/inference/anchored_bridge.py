"""Replace Dovekie by the anchored SN factor only after parent qualification."""
import argparse
import copy
import importlib.metadata
import json
import os
from pathlib import Path
import time
from unittest.mock import patch

import numpy as np
from scipy.special import logsumexp

import anchored_adapter as anchored
from expansion_history import native_thermal_parameters
from luminosity_sensitivity import IntegratedLuminosity, weight_diagnostics
from measurement_summary import PARAMETERS, summarize_run, weighted_fraction
from modern_fast import configuration as parent_configuration
from probe_omission import identity, payload, verify_hashes
from target_identity import canonical, digest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DESIGN = HERE/'anchored-bridge-design.json'
GATES = HERE/'luminosity-sensitivity-design.json'


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def configurations(settings, frozen, design):
    assert settings['model'] in design['allowed_models']
    assert settings['evolution'] == design['parent_evolution']
    assert settings['sample'] == design['parent_sample']
    assert settings['calibration'] == design['parent_calibration']
    assert settings.get('fast_lensing') is True and not settings.get('gpu', False)
    options = {k: settings[k] for k in ['model', 'evolution', 'sample', 'calibration']}
    proposal = parent_configuration(**options, surrogate=Path(settings['surrogate']))
    assert canonical(proposal) == frozen['configuration'], 'Parent proposal configuration changed.'
    old = parent_configuration(**options)
    new = anchored.configuration(model=settings['model'])
    left, right = copy.deepcopy(canonical(old)), copy.deepcopy(canonical(new))
    old_sn = left['likelihood'].pop('released_sn')
    right['likelihood'].pop('released_sn')
    assert left == right, 'A non-SN likelihood, prior or physical theory changed.'
    assert old_sn['external'] == 'likelihood.ReleasedDistances'
    assert old_sn.get('smooth_sigma', 0.) == 0. and left['params']['epsilon'] == 0.
    return old, new


def load_SN_interfaces(data_path):
    with np.load(data_path, allow_pickle=False) as data:
        old = IntegratedLuminosity(*(data[k] for k in ['zHD', 'zHEL', 'MU', 'covariance']))
    return old, anchored.ReleasedCalibration()


def background_densities(point, extra, old, new):
    """Single physical background, no theory/provider likelihood construction."""
    import camb
    cosmology = {k: point[k] for k in ['H0', 'ombh2', 'omch2', 'ns', 'tau']}
    cosmology.update(As=1e-10*np.exp(point['logA']), w=point.get('w', -1.), wa=point.get('wa', 0.))
    z = np.concatenate([old.z, new.z_hd_noncalibrator])
    unique, inverse = np.unique(z, return_inverse=True)
    def forbidden(*args, **kwargs):
        raise RuntimeError('Anchored bridge forbids CMB spectra and transfer functions.')
    with patch.object(camb, 'get_results', forbidden), patch.object(camb, 'get_transfer_functions', forbidden):
        parameters = camb.set_params(**cosmology, **extra)
        adjusted = native_thermal_parameters(parameters)
        background = camb.get_background(adjusted)
        distances = np.asarray(background.angular_diameter_distance(unique))[inverse]
    assert distances.shape == z.shape and np.isfinite(distances).all() and np.all(distances > 0)
    prediction = 5*np.log10(distances[:old.n]*(1+old.z)*(1+old.zhel))+25
    original = old.from_prediction(prediction)['baseline_SN_loglike']
    replacement = new.evaluate(distances[old.n:])
    return {'old_SN_reconstructed_loglike': float(original), 'anchored_SN': replacement,
            'unique_background_redshifts': len(unique),
            'old_ordered_rows': old.n, 'new_ordered_noncalibrator_rows': len(new.z_hd_noncalibrator),
            'distance_unit': 'Mpc', 'H0_unit': 'km/s/Mpc',
            'thermal_adapter': {'input_WantTransfer': bool(parameters.WantTransfer),
                                'background_WantTransfer': bool(adjusted.WantTransfer)},
            'CMB_spectrum_calls': 0, 'background_calls': 1}


def replace_SN(record, background, components, design):
    """The actual weight uses the stored native source density, after closure."""
    exact, proposal = record['exact_loglikes'], record['proposal_loglikes']
    assert record['status'] == 'finite' and set(exact) == set(proposal) == set(components)
    tolerance = design['density_accounting_absolute_tolerance']
    assert abs(record['log_weight']-record['exact_logpost']+record['proposal_logpost']) <= tolerance
    old_prior = float(record['exact_logpost']-sum(exact.values()))
    proposal_prior = float(record['proposal_logpost']-sum(proposal.values()))
    assert abs(old_prior-proposal_prior) <= tolerance, 'Parent native/proposal priors differ.'
    source = float(exact['released_sn'])
    check = float(background['old_SN_reconstructed_loglike'])-source
    target = float(background['anchored_SN']['loglike'])
    ratio = target-source
    logweight = float(record['log_weight']+ratio)
    new_components = dict(exact, released_sn=target)
    new_logpost = old_prior+sum(new_components.values())
    closure = new_logpost-record['proposal_logpost']-logweight
    assert abs(closure) <= tolerance
    values = [source, check, target, ratio, logweight, new_logpost, closure, old_prior]
    assert np.isfinite(values).all(), 'Nonfinite bridge density.'
    return {'status': 'finite_anchored_replacement' if abs(check) <= design['source_closure_absolute_tolerance'] else 'source_SN_density_mismatch',
            'source_SN_stored_native_loglike': source, 'source_SN_loglike_closure': check,
            'target_SN_loglike': target, 'SN_log_likelihood_ratio': ratio,
            'parent_native_proposal_logweight': float(record['log_weight']),
            'target_logweight': logweight, 'unchanged_logprior': old_prior,
            'target_logpost_in_inherited_normalization': new_logpost,
            'target_loglikes': new_components, 'density_accounting_closure': float(closure),
            'target_sn_chi2': float(background['anchored_SN']['chi2'])}


def summarize(records, rows, groups, gates):
    assert len(records) == len(rows) == len(groups)
    failures = [{'index': i, 'status': row['status'], 'error': row.get('error')}
                for i, row in enumerate(rows) if row['status'] != 'finite_anchored_replacement']
    empty = dict(posterior=None, conditional_sign_fractions=None, weighted_covariance=None,
                 paired_mean_changes_from_parent=None, qualified_under_declared_numerical_gates=False)
    if failures:
        return dict(empty, status='failed_anchored_bridge_evaluation', failures=failures,
                    failed_gates=sorted({r['status'] for r in failures}))
    groups = np.asarray(groups)
    assert set(groups) == {0, 1, 2, 3}
    names = sorted(set(records[0]['point']) | set(records[0]['derived']))
    values = {name: np.array([dict(r['point'], **r['derived'])[name] for r in records]) for name in names}
    # The old SN chi-square is not the new target's goodness-of-fit diagnostic.
    values['sn_chi2'] = np.array([r['target_sn_chi2'] for r in rows])
    assert all(np.isfinite(x).all() for x in values.values())
    lw = np.array([row['target_logweight'] for row in rows])
    diagnostic = weight_diagnostics(lw, values, groups, gates)
    failed = list(diagnostic['failed_gates'])
    if len(rows) < gates['minimum_exact_points']:
        failed.append('minimum_exact_points')
    result = dict(empty, status='insufficient_anchored_overlap_or_stability', points=len(rows),
                  failed_gates=failed, weight_diagnostics=diagnostic,
                  source_SN_max_absolute_loglike_closure=max(abs(r['source_SN_loglike_closure']) for r in rows),
                  maximum_density_accounting_closure=max(abs(r['density_accounting_closure']) for r in rows),
                  SN_log_likelihood_ratio_quantiles=np.quantile([r['SN_log_likelihood_ratio'] for r in rows], [0,.025,.5,.975,1]).tolist())
    if not failed:
        weights = np.exp(lw-logsumexp(lw))
        parent_lw = np.array([r['log_weight'] for r in records])
        parent_weights = np.exp(parent_lw-logsumexp(parent_lw))
        published = [name for name in PARAMETERS if name in values]
        matrix = np.column_stack([values[name] for name in published])
        centered = matrix-weights@matrix
        result.update(status='qualified_conditional_anchored_bridge', qualified_under_declared_numerical_gates=True,
            posterior={name: diagnostic['weighted_summaries_for_diagnostics'][name] for name in published},
            conditional_sign_fractions={name+'_below_zero': weighted_fraction(values[name] < 0, weights, groups)
                                        for name in ['q0', 'q05', 'q1', 'j0'] if name in values},
            weighted_covariance={'parameter_order': published, 'matrix': ((centered*weights[:, None]).T@centered).tolist(),
                                 'normalization': 'Raw normalized importance weights; no Bessel correction.', 'qualified': True},
            paired_mean_changes_from_parent={name: {'parent_mean': float(parent_weights@values[name]),
                'target_mean': float(weights@values[name]), 'difference': float((weights-parent_weights)@values[name])}
                for name in published})
    return result


def read_record(path, cache_identity, index, point, native_hash):
    row = json.loads(Path(path).read_text())
    checksum = row.pop('payload_sha256')
    assert identity(row) == checksum, 'Cached bridge payload changed.'
    assert row['identity'] == cache_identity and row['index'] == index and row['point'] == point
    assert row['native_record_sha256'] == native_hash
    row['payload_sha256'] = checksum
    return row


def actual(folder, summary_path, cache):
    # No provisional-parent escape hatch: qualification precedes ALL work.
    parent = summarize_run(folder, summary_path)
    assert all(os.environ.get(k) == '1' for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'CLIPY_NOJAX']), 'Use the declared one-thread NumPy runtime.'
    folder, summary_path, cache = map(lambda p: Path(p).resolve(), [folder, summary_path, cache])
    assert cache.is_relative_to(ROOT/'.work')
    verify_hashes(parent['input_sha256'])
    design = json.loads(DESIGN.read_text()); gates = json.loads(GATES.read_text())['overlap_gates']
    summary = json.loads(summary_path.read_text())
    selection_path = ROOT/summary['selection_path']
    assert selection_path.parent.parent == folder
    selection = json.loads(selection_path.read_text()); settings = selection['settings']
    frozen = json.loads((folder/'run-0.json').read_text())['target_identity']
    old_config, new_config = configurations(settings, frozen, design)
    old_path = Path(old_config['likelihood']['released_sn']['data_file']).resolve()
    assert digest(old_path) == frozen['sample_sha256'], 'Parent SN data changed.'
    target = anchored.identify(new_config)
    assert target['assets'] == frozen['assets'], 'A supposedly cancelling external asset changed.'
    assert all(target['versions'][k] == v for k, v in frozen['versions'].items())
    sources = [Path(__file__), DESIGN, GATES, HERE/'anchored_adapter.py', HERE/'anchored-design.json',
               HERE/'expansion_history.py', HERE/'luminosity_sensitivity.py', HERE/'measurement_summary.py',
               HERE/'probe_omission.py', HERE/'exact_correction.py', HERE/'target_identity.py',
               HERE/'modern_fast.py', HERE/'modern_run.py', HERE.parent/'distance_ladder/calibration_interface.py']
    hashes = {relative(p): digest(p) for p in sources}
    lineage = {'qualified_parent_inputs': parent['input_sha256'], 'source_sha256': hashes,
               'parent_proposal_target_identity': frozen['identity'], 'native_target': target,
               'parent_settings': settings, 'parent_SN_sha256': digest(old_path),
               'all_non_SN_factors_and_priors_identical': True}
    cache_identity = identity(lineage)
    old, new = load_SN_interfaces(old_path)
    cache.mkdir(parents=True, exist_ok=True)
    lineage_path = cache/'lineage.json'; payload(lineage_path, lineage)
    ledger = cache/'record-hashes.json'
    if ledger.exists():
        verify_hashes(json.loads(ledger.read_text()))
    records, rows = [], []
    hashes_out = {relative(lineage_path): digest(lineage_path)}
    extra = old_config['theory']['camb']['extra_args']
    start = time.monotonic(); new_calls = 0
    for index, point in enumerate(selection['points']):
        native_path = selection_path.parent/f'{index:05d}.json'
        record = json.loads(native_path.read_text()); records.append(record)
        path = cache/f'{index:05d}.json'
        if path.exists():
            row = read_record(path, cache_identity, index, point, digest(native_path))
        else:
            try:
                new_calls += 1
                background = background_densities(point, extra, old, new)
                row = replace_SN(record, background, old_config['likelihood'], design)
                row['background'] = background
            except Exception as error:
                row = {'status': 'failed_background_or_density', 'error': repr(error)}
            row.update(identity=cache_identity, index=index, point=point, native_record_sha256=digest(native_path))
            row['payload_sha256'] = identity(row)
            payload(path, row)
        rows.append(row); hashes_out[relative(path)] = digest(path)
        if (index+1) % 100 == 0:
            print(json.dumps({'anchored_bridge_points': index+1, 'seconds': time.monotonic()-start}), flush=True)
    payload(ledger, hashes_out)
    result = summarize(records, rows, selection['groups'], gates)
    result.update(parent_settings=parent['settings'], native_target=target,
        source_sha256=hashes, bridge_identity=cache_identity,
        lineage_path=relative(lineage_path), lineage_sha256=digest(lineage_path),
        cache_manifest_path=relative(ledger), cache_manifest_sha256=digest(ledger),
        correction_summary_path=relative(summary_path), correction_summary_sha256=digest(summary_path),
        parent_qualified_under_declared_numerical_gates=True, CMB_spectrum_calls=0,
        background_calls_in_recorded_cohort=sum(r.get('background', {}).get('background_calls', 0) for r in rows),
        runtime_environment={k:os.environ.get(k) for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','CLIPY_NOJAX']},
        cached_or_computed_points=len(rows), limits=design['limits'],
        normalization=design['calibration'],
        units={'H0': 'km/s/Mpc', 'rdrag': 'Mpc', 'w_wa_q_j': 'dimensionless'},
        comparison='Paired conditional mean shifts, not independent-fit significances. Dovekie is replaced, not multiplied by an H0 factor.',
        status_of_failed_weight_summaries='Diagnostics only, never an admitted measurement.')
    for mapping in [hashes, parent['input_sha256'], frozen['source_sha256'],
                    target['source_sha256'], target['anchored_calibration']['input_and_audit_sha256'],
                    {relative(old_path): frozen['sample_sha256']}]:
        verify_hashes(mapping)
    assert all(importlib.metadata.version(k) == v for k, v in target['versions'].items())
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--chain-folder', type=Path, required=True)
    parser.add_argument('--correction-summary', type=Path, required=True)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = actual(args.chain_folder, args.correction_summary, args.cache)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    payload(args.output, result)
    print(json.dumps({'status': result['status'], 'points': result.get('points'), 'CMB_spectrum_calls': 0}))


if __name__ == '__main__':
    main()
