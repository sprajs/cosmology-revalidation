"""Assemble qualified native-CAMB measurements and their numerical limits.

No new data factors, smoothing, clipping or model probabilities are introduced.
Every sign fraction is conditional on the model and luminosity prior named in
its row. Zero sampled failures never become a claim of certain acceleration.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

ROOT = Path(__file__).resolve().parents[4]
PARAMETERS = ['H0', 'omegam', 'ombh2', 'omch2', 'rdrag', 'logA', 'ns', 'tau',
              'w', 'wa', 'epsilon', 'q0', 'q05', 'q1', 'j0', 'j05', 'j1']


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def weighted_fraction(mask, weights, groups):
    value = np.asarray(mask, dtype=float)
    fraction = float(weights @ value)
    chain_fractions = {}
    for group in np.unique(groups):
        keep = groups == group
        mass = weights[keep].sum()
        assert mass > 0
        chain_fractions[str(group)] = float(weights[keep] @ value[keep] / mass)
    batch_mcse = {}
    for count in [10, 20, 40]:
        scores = []
        masses = []
        for group in np.unique(groups):
            for indices in np.array_split(np.flatnonzero(groups == group), count):
                scores.append(weights[indices] @ (value[indices] - fraction))
                masses.append(weights[indices].sum())
        batch_mcse[str(count)] = float(np.sqrt(np.var(scores, ddof=1) / len(scores))
                                       / np.mean(masses))
    return {'fraction': fraction, 'independent_chain_fractions': chain_fractions,
            'batch_mean_mcse': batch_mcse,
            'both_events_observed': bool(value.min() != value.max()),
            'interpretation': 'Conditional Monte Carlo fraction. A zero tail count '
                              'or zero estimated MCSE is not certainty or an exclusion bound.'}


def summarize_run(folder, summary_path):
    folder = Path(folder).resolve()
    summary_path = Path(summary_path).resolve()
    summary = json.loads(summary_path.read_text())
    assert summary['status'] == 'passed_importance_weight_gates'
    selection_path = (ROOT / summary['selection_path']).resolve()
    assert selection_path.name == 'selection.json' and selection_path.parent.parent == folder
    selection = json.loads(selection_path.read_text())
    assert digest(selection_path) == summary['selection_sha256']
    assert relative(selection_path) == summary['selection_path']
    assert digest(Path(__file__).with_name('exact_correction.py')) == summary['code_sha256']
    manifests = [folder / f'run-{rank}.json' for rank in range(4)]
    manifest = json.loads(manifests[0].read_text())
    for path in manifests:
        assert json.loads(path.read_text())['target_identity'] == manifest['target_identity']
    inputs = {relative(p): digest(p) for p in manifests + [selection_path, summary_path]}
    dependencies = selection['correction_dependency_sha256']
    assert dependencies == summary['correction_dependency_sha256']
    for name, expected in dependencies.items():
        assert digest(ROOT / name) == expected, 'Correction dependency changed.'
        inputs[name] = expected
    target = manifest['target_identity']['identity']
    identity = hashlib.sha256(json.dumps(
        {'target_identity': target, 'correction_dependencies': dependencies},
        sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    assert identity == selection['correction_identity'] == summary['correction_identity']
    assert target == selection['proposal_target_identity'] == summary['proposal_target_identity']
    diagnostic_path = ROOT / selection['diagnostics_path']
    assert relative(diagnostic_path) == summary['diagnostics_path']
    assert digest(diagnostic_path) == selection['diagnostics_sha256'] == summary['diagnostics_sha256']
    check = json.loads(diagnostic_path.read_text())
    assert check['status'] == 'passed' and not check['failed_gates']
    assert check['input_sha256'] == selection['chain_sha256']
    assert check['code_sha256'] == digest(Path(__file__).with_name('diagnostics.py'))
    assert check['discard_fraction_after_sampler_burnin'] == .3
    inputs[relative(diagnostic_path)] = digest(diagnostic_path)
    for name, expected in selection['chain_sha256'].items():
        assert digest(ROOT / name) == expected, 'Selected chain bytes changed.'
        inputs[name] = expected
    records = []
    from exact_correction import verify_record
    expected_records = {relative(selection_path.parent / f'{index:05d}.json')
                        for index in range(len(selection['points']))}
    assert set(summary['native_record_sha256']) == expected_records
    for index, point in enumerate(selection['points']):
        path = selection_path.parent / f'{index:05d}.json'
        assert digest(path) == summary['native_record_sha256'][relative(path)], 'Native file changed after correction summary.'
        record = json.loads(path.read_text())
        verify_record(record)
        assert record['index'] == index and record['point'] == point
        assert record['target_identity'] == selection['correction_identity']
        assert record['status'] == 'finite'
        assert abs(record['log_weight'] - record['exact_logpost']
                   + record['proposal_logpost']) < 1e-9
        assert abs(record['proposal_chain_logpost_difference']) <= .01
        inputs[relative(path)] = digest(path)
        records.append(record)
    groups = np.asarray(selection['groups'])
    assert set(groups) == {0, 1, 2, 3} and len(records) >= 2000
    # Recompute all gates and weighted summaries from the actual native records;
    # a status label alone cannot authorize a measurement row.
    from exact_correction import summarize
    recomputed = summarize(records, groups)
    assert recomputed == {key: summary[key] for key in recomputed}
    assert not recomputed['failed_gates']
    lw = np.asarray([row['log_weight'] for row in records])
    weights = np.exp(lw - logsumexp(lw))
    assert abs(1 / (weights @ weights) - summary['raw_weight_ess']) < 1e-7
    values = {name: np.asarray([dict(row['point'], **row['derived'])[name]
                               for row in records])
              for name in summary['posterior']}
    for name, x in values.items():
        assert np.isfinite(x).all()
        assert abs(weights @ x - summary['posterior'][name]['mean']) < 1e-8
    signs = {f'{key}_below_zero': weighted_fraction(values[key] < 0, weights, groups)
             for key in ['q0', 'q05', 'q1', 'j0']}
    signs['accelerating_with_decreasing_scale_factor_acceleration_today'] = weighted_fraction(
        (values['q0'] < 0) & (values['j0'] < 0), weights, groups)
    correlations = {}
    for left, right in [('w', 'wa'), ('w', 'H0'), ('wa', 'H0'),
                        ('epsilon', 'w'), ('epsilon', 'wa'), ('q0', 'j0')]:
        if left not in values or right not in values:
            continue
        x = values[left] - weights @ values[left]
        y = values[right] - weights @ values[right]
        denominator = np.sqrt((weights @ (x*x)) * (weights @ (y*y)))
        correlations[left + ',' + right] = float(weights @ (x*y) / denominator) if denominator else None
    boundaries = {}
    params = manifest['target_identity']['configuration']['params']
    for name, definition in params.items():
        if name not in values or not isinstance(definition, dict):
            continue
        prior = definition.get('prior', {})
        if not isinstance(prior, dict) or not {'min', 'max'} <= set(prior):
            continue
        lo, hi = prior['min'], prior['max']
        margin = .01 * (hi - lo)
        boundaries[name] = {'prior': [lo, hi], 'fraction_within_one_percent_of_lower_bound':
                            float(weights @ (values[name] < lo + margin)),
                            'fraction_within_one_percent_of_upper_bound':
                            float(weights @ (values[name] > hi - margin))}
    support = None
    if 'w' in values and 'wa' in values:
        early_w = values['w'] + values['wa']
        assert np.max(early_w) <= 1e-10
        support = {'constraint': 'w0+wa<=0', 'diagnostic_boundary_width': .05,
                   'fraction_within_005_of_boundary': float(weights @ (early_w > -.05)),
                   'interpretation': 'Diagnostic for the separately imposed CAMB support; '
                                     'not a test of the excluded region.'}
    covariance_names = [name for name in PARAMETERS if name in values]
    matrix = np.column_stack([values[name] for name in covariance_names])
    centered = matrix - weights @ matrix
    covariance = (centered * weights[:, None]).T @ centered
    return {'settings': {key: selection['settings'][key]
                         for key in ['model', 'evolution', 'sample', 'calibration']},
            'target_identity': manifest['target_identity']['identity'],
            'qualified_under_declared_numerical_gates': True,
            'native_correction_summary': relative(summary_path),
            'posterior': {name: summary['posterior'][name] for name in PARAMETERS
                          if name in summary['posterior']},
            'conditional_sign_fractions': signs, 'weighted_correlations': correlations,
            'weighted_covariance': {'parameter_order': covariance_names,
                                    'matrix': covariance.tolist(),
                                    'normalization': 'Posterior sum w_i (x_i-mean)(x_i-mean)^T; weights sum to one.'},
            'scalar_box_prior_boundary_fractions': boundaries,
            'coupled_CAMB_support_boundary': support,
            'native_correction': {key: summary[key] for key in
                                  ['exact_points', 'raw_weight_ess', 'pareto_k',
                                   'largest_normalized_weight', 'log_weight_quantiles']},
            'input_sha256': inputs}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--run', nargs=2, action='append', required=True,
                   metavar=('CHAIN_FOLDER', 'CORRECTION_SUMMARY'))
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    measurements = [summarize_run(*run) for run in a.run]
    keys = [tuple(row['settings'].items()) for row in measurements]
    assert len(keys) == len(set(keys)), 'Duplicate physical target in report.'
    report = {'status': 'qualified_conditional_measurements', 'measurements': measurements,
              'scope': 'Shared flat-GR CMB, BAO and released-supernova cosmology. '
                       'CAMB analytic CPL support w0+wa<=0, fixed neutrino physics, '
                       'released-distance corrections and stated probe factorization. '
                       'Luminosity priors are sensitivities, not measured age corrections.',
              'reporting': 'No Gaussian sigma conversion of sign fractions, evidence ratio, '
                           'global mode guarantee or empirical luminosity-age coefficient.',
              'source_sha256': digest(__file__)}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'status': report['status'], 'targets': [r['settings'] for r in measurements]}, indent=2))


if __name__ == '__main__':
    main()
