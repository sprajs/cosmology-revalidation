"""Qualified-parent component removal using stored exact likelihoods only."""
import argparse
import copy
import hashlib
import importlib.metadata
import json
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from luminosity_sensitivity import weight_diagnostics
from measurement_summary import PARAMETERS, summarize_run, weighted_fraction
from target_identity import canonical

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DESIGN = HERE/'probe-omission-design.json'
GATES = HERE/'luminosity-sensitivity-design.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def identity(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def configuration_factory(settings):
    # Building a configuration does not initialize a likelihood, theory or GPU.
    if settings.get('gpu'):
        from modern_gpu import configuration
    elif settings.get('fast_lensing'):
        from modern_fast import configuration
    else:
        from modern_run import configuration
    return configuration


def target_configuration(settings, frozen, omission, design):
    assert settings['model'] in design['allowed_models']
    assert settings['evolution'] in design['allowed_evolution']
    assert settings['calibration'] in design['allowed_calibration']
    assert omission in design['alternatives'], 'Only the three declared omissions are allowed.'
    factory = configuration_factory(settings)
    options = {key: settings[key] for key in ['model', 'evolution', 'sample', 'calibration']}
    proposal = canonical(factory(**options, surrogate=Path(settings['surrogate'])))
    assert proposal == frozen['configuration'], 'Parent proposal configuration changed.'
    native = canonical(factory(**options))
    expected = set(design['required_parent_components'])
    assert set(native['likelihood']) == expected, 'Unexpected or missing parent probe component.'
    lensing = native['likelihood']['act_dr6_lenslike.ACTDR6LensLike']
    assert lensing['variant'] == 'actplanck_baseline' and lensing['version'] == 'v1.2'
    assert lensing['lens_only'] is False and lensing['no_like_corrections'] is False
    pan = native['likelihood']['SPT2023_lensing']
    assert pan['data_set_file'] == 'candl_data.SPT3G_2018_Lens_and_CMB'
    assert pan['lensing'] is True and pan['clear_internal_priors'] is True
    assert native['likelihood']['released_sn']['external'] == 'likelihood.ReleasedDistances'
    removed = design['alternatives'][omission]['removed_components']
    target = copy.deepcopy(native)
    for name in removed:
        target['likelihood'].pop(name)
    assert target['params'] == native['params'], 'Prior or auxiliary dimension changed.'
    assert target['theory'] == native['theory'], 'Physical theory changed.'
    assert set(target['likelihood']) == expected-set(removed)
    return native, target, list(removed)


def remove_components(record, components, expected_components, tolerance):
    """Pure arithmetic kernel; production components come only from the design."""
    assert len(components) == len(set(components))
    exact = record['exact_loglikes']; proposal = record['proposal_loglikes']
    assert set(exact) == set(proposal) == set(expected_components)
    assert set(components) <= set(exact)
    assert record['status'] == 'finite'
    assert np.isfinite(list(exact.values())+list(proposal.values())+
                       [record[k] for k in ['exact_logpost', 'proposal_logpost', 'log_weight']]).all()
    original = float(record['log_weight'])
    assert abs(original-record['exact_logpost']+record['proposal_logpost']) <= tolerance
    exact_prior = float(record['exact_logpost']-sum(exact.values()))
    proposal_prior = float(record['proposal_logpost']-sum(proposal.values()))
    assert abs(exact_prior-proposal_prior) <= tolerance, 'Native/proposal prior accounting changed.'
    removed = {name: float(exact[name]) for name in components}
    removed_sum = float(sum(removed.values()))
    target_logpost = float(record['exact_logpost']-removed_sum)
    logweight = float(original-removed_sum)
    closure = float(target_logpost-record['proposal_logpost']-logweight)
    assert abs(closure) <= tolerance
    assert np.isfinite([target_logpost, logweight, removed_sum]).all()
    return {'status': 'finite_component_removal', 'removed_native_loglikes': removed,
            'removed_native_loglike_sum': removed_sum,
            'original_native_proposal_logweight': original, 'target_logweight': logweight,
            'target_logpost_in_inherited_normalization': target_logpost,
            'unchanged_inferred_logprior': exact_prior,
            'native_proposal_prior_closure': exact_prior-proposal_prior,
            'target_weight_accounting_closure': closure}


def summarize_omission(records, rows, groups, gates, unused_auxiliary_coordinates=()):
    assert len(records) == len(rows) == len(groups)
    bad = [{'index': i, 'status': row['status'], 'error': row.get('error')}
           for i, row in enumerate(rows) if row['status'] != 'finite_component_removal']
    if bad:
        return {'status': 'failed_component_accounting', 'qualified_under_declared_numerical_gates': False,
                'failures': bad, 'posterior': None, 'conditional_sign_fractions': None,
                'weighted_covariance': None}
    groups = np.asarray(groups)
    assert set(groups) == {0, 1, 2, 3}
    names = sorted(set(records[0]['point']) | set(records[0]['derived']))
    values = {name: np.array([dict(row['point'], **row['derived'])[name] for row in records])
              for name in names}
    assert all(np.isfinite(v).all() for v in values.values())
    logweights = np.array([row['target_logweight'] for row in rows])
    assert np.isfinite(logweights).all()
    diagnostic = weight_diagnostics(logweights, values, groups, gates)
    failed = list(diagnostic['failed_gates'])
    if len(rows) < gates['minimum_exact_points']:
        failed.append('minimum_exact_points')
    weights = np.exp(logweights-logsumexp(logweights))
    parent_lw = np.array([r['log_weight'] for r in records])
    parent_weights = np.exp(parent_lw-logsumexp(parent_lw))
    posterior = signs = covariance = changes = None
    if not failed:
        # Prior-only auxiliaries remain in every overlap/stability check above,
        # but their sampled prior is not published as a measured constraint.
        published = [name for name in PARAMETERS
                     if name in values and name not in unused_auxiliary_coordinates]
        posterior = {name: diagnostic['weighted_summaries_for_diagnostics'][name] for name in published}
        signs = {name+'_below_zero': weighted_fraction(values[name] < 0, weights, groups)
                 for name in ['q0', 'q05', 'q1', 'j0'] if name in values}
        if 'q0' in values and 'j0' in values:
            signs['accelerating_with_decreasing_scale_factor_acceleration_today'] = weighted_fraction(
                (values['q0'] < 0) & (values['j0'] < 0), weights, groups)
        matrix = np.column_stack([values[name] for name in published])
        centered = matrix-weights@matrix
        covariance = {'parameter_order': published, 'matrix': ((centered*weights[:, None]).T@centered).tolist(),
                      'normalization': 'Untrimmed normalized weights, no Bessel correction.'}
        changes = {name: {'parent_mean': float(parent_weights@values[name]),
                          'target_mean': float(weights@values[name]),
                          'difference': float((weights-parent_weights)@values[name])}
                   for name in published}
    return {'status': 'qualified_conditional_probe_omission' if not failed else 'insufficient_omission_overlap_or_stability',
            'qualified_under_declared_numerical_gates': not failed,
            'points': len(rows), 'failed_gates': failed, 'weight_diagnostics': diagnostic,
            'posterior': posterior, 'conditional_sign_fractions': signs,
            'weighted_covariance': covariance, 'paired_mean_changes_from_parent': changes,
            'changes_interpretation': 'Paired point estimates sharing the parent cohort, not independent-fit significances.',
            'maximum_normalized_weight_change_from_parent': float(np.max(abs(weights-parent_weights))),
            'maximum_target_weight_accounting_closure': max(abs(r['target_weight_accounting_closure']) for r in rows),
            'removed_native_loglike_sum_quantiles': np.quantile([r['removed_native_loglike_sum'] for r in rows], [0,.025,.5,.975,1]).tolist()}


def verify_hashes(mapping):
    for path, expected in mapping.items():
        assert digest(ROOT/path) == expected, f'Input/source identity changed: {path}'


def payload(path, value):
    if path.exists():
        assert json.loads(path.read_text()) == value, 'Previously cached omission payload changed.'
    else:
        temporary = path.with_suffix('.part')
        temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
        temporary.replace(path)


def actual(folder, summary_path, omission, cache):
    # This is the sole public entry: no provisional or qualification bypass.
    parent = summarize_run(folder, summary_path)
    folder = Path(folder).resolve(); summary_path = Path(summary_path).resolve()
    cache = Path(cache).resolve()
    assert cache.is_relative_to(ROOT/'.work')
    verify_hashes(parent['input_sha256'])
    design = json.loads(DESIGN.read_text())
    gates = json.loads(GATES.read_text())['overlap_gates']
    summary = json.loads(summary_path.read_text())
    selection_path = ROOT/summary['selection_path']
    assert selection_path.parent.parent == folder
    selection = json.loads(selection_path.read_text()); settings = selection['settings']
    frozen = json.loads((folder/'run-0.json').read_text())['target_identity']
    verify_hashes(frozen['source_sha256'])
    for name, version in frozen['versions'].items():
        assert importlib.metadata.version(name) == version, 'Parent numerical environment changed.'
    native, target, removed = target_configuration(settings, frozen, omission, design)
    data_path = Path(native['likelihood']['released_sn']['data_file']).resolve()
    assert digest(data_path) == frozen['sample_sha256'], 'Parent SN data changed.'
    sources = [Path(__file__), DESIGN, GATES, HERE/'measurement_summary.py',
               HERE/'exact_correction.py', HERE/'luminosity_sensitivity.py', HERE/'target_identity.py',
               HERE/'modern_run.py', HERE/'modern_fast.py']
    if settings.get('gpu'):
        sources.append(HERE/'modern_gpu.py')
    hashes = {relative(path): digest(path) for path in sources}
    native_target = {'configuration': target, 'source_sha256': frozen['source_sha256'],
                     'versions': frozen['versions'], 'parent_likelihood_assets': frozen['assets'],
                     'parent_sample_sha256': frozen['sample_sha256'], 'omission': omission,
                     'removed_components': removed, 'all_parent_prior_dimensions_retained': True,
                     'inherited_support': 'Analytic CAMB CPL w0+wa<=0 when CPL; unchanged flat GR and fixed-neutrino physics.'}
    native_target['identity'] = identity(native_target)
    sampled = {name: value['prior'] for name, value in target['params'].items()
               if isinstance(value, dict) and 'prior' in value}
    auxiliary = []
    if omission == 'sn' and 'epsilon' in sampled:
        auxiliary.append('epsilon')
    if omission == 'lensing' and 'A_fg' in sampled:
        auxiliary.append('A_fg')
    lineage = {'qualified_parent_inputs': parent['input_sha256'], 'source_sha256': hashes,
               'parent_proposal_target_identity': frozen['identity'], 'native_target': native_target,
               'omission': omission, 'parent_settings': parent['settings'],
               'retained_sampled_priors': sampled, 'unused_auxiliary_coordinates': auxiliary,
               'parent_SN_data_sha256': frozen['sample_sha256']}
    cache_identity = identity(lineage)
    cache.mkdir(parents=True, exist_ok=True)
    lineage_path = cache/'lineage.json'; payload(lineage_path, lineage)
    ledger_path = cache/'record-hashes.json'
    if ledger_path.exists():
        verify_hashes(json.loads(ledger_path.read_text()))
    records = []; rows = []; hashes_out = {relative(lineage_path): digest(lineage_path)}
    for index, point in enumerate(selection['points']):
        native_path = selection_path.parent/f'{index:05d}.json'
        record = json.loads(native_path.read_text()); records.append(record)
        try:
            row = remove_components(record, removed, native['likelihood'],
                                    design['logweight_accounting_absolute_tolerance'])
        except (AssertionError, ArithmeticError, KeyError, TypeError, ValueError) as error:
            row = {'status': 'failed_component_accounting', 'error': repr(error)}
        row.update(index=index, point=point, identity=cache_identity,
                   native_record_sha256=digest(native_path))
        path = cache/f'{index:05d}.json'; payload(path, row)
        rows.append(row); hashes_out[relative(path)] = digest(path)
    payload(ledger_path, hashes_out)
    result = summarize_omission(records, rows, selection['groups'], gates, auxiliary)
    result.update(parent_settings=parent['settings'], omission=omission,
                  retained_observations=design['alternatives'][omission]['retained_observations'],
                  removed_components=removed, native_target=native_target, omission_identity=cache_identity,
                  retained_sampled_priors=sampled, unused_auxiliary_coordinates=auxiliary,
                  source_sha256=hashes, lineage_path=relative(lineage_path), lineage_sha256=digest(lineage_path),
                  cache_manifest_path=relative(ledger_path), cache_manifest_sha256=digest(ledger_path),
                  correction_summary_path=relative(summary_path), correction_summary_sha256=digest(summary_path),
                  parent_qualified_under_declared_numerical_gates=True,
                  scope=design['scope'], limitations=design['limitations'],
                  normalization='Original component constants and all proper priors retained. '
                                'A removed SN flat-M normalization constant cancels in normalized weights; no evidence ratio is defined.',
                  CMB_spectrum_calls=0, background_calls=0, GPU_kernel_evaluations=0,
                  hardware_identity_queries='The parent qualifier may query CUDA/device versions for a GPU parent; no GPU arithmetic is performed.',
                  status_of_failed_weight_summaries='Diagnostic only; not an admitted posterior measurement.')
    for mapping in [hashes, frozen['source_sha256'], parent['input_sha256'],
                    {relative(data_path): frozen['sample_sha256']}]:
        verify_hashes(mapping)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--chain-folder', type=Path, required=True)
    parser.add_argument('--correction-summary', type=Path, required=True)
    parser.add_argument('--omit', choices=['sn', 'bao', 'lensing'], required=True)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = actual(args.chain_folder, args.correction_summary, args.omit, args.cache)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': result['status'], 'omission': args.omit,
                      'points': result.get('points'), 'CMB_spectrum_calls': 0, 'background_calls': 0}))


if __name__ == '__main__':
    main()
