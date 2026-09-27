"""Synthetic removal/qualification/cache checks without native physics calls."""
import copy
import importlib.metadata
import json
from pathlib import Path
import tempfile
from unittest.mock import patch

import numpy as np
from scipy.special import logsumexp
from scipy.stats import norm

import probe_omission as audit


def must_fail(operation, contains=None):
    try:
        operation()
    except (AssertionError, ValueError) as error:
        if contains is not None:
            assert contains in str(error), str(error)
        return
    raise AssertionError('Deliberately invalid fixture was accepted.')


def run_tests():
    design = json.loads(audit.DESIGN.read_text())
    gates = json.loads(audit.GATES.read_text())['overlap_gates']
    sources = list(design['required_parent_components'])
    config_cases = 0
    fixture_model = audit.ROOT/'.work/unified-cosmology/inference/quartic-spectral/surrogate-quartic-v1.npz'
    for backend in [{}, {'fast_lensing': True}, {'gpu': True, 'fast_lensing': True}]:
        for model in ['lcdm', 'cpl']:
            for evolution in design['allowed_evolution']:
                for calibration in design['allowed_calibration']:
                    settings = dict(backend, model=model, evolution=evolution, calibration=calibration,
                                    sample='dovekie', surrogate=str(fixture_model))
                    factory = audit.configuration_factory(settings)
                    options = {k: settings[k] for k in ['model', 'evolution', 'sample', 'calibration']}
                    frozen = {'configuration': audit.canonical(factory(**options, surrogate=fixture_model))}
                    for omit in design['alternatives']:
                        native, target, removed = audit.target_configuration(settings, frozen, omit, design)
                        assert target['params'] == native['params'] and target['theory'] == native['theory']
                        assert set(native['likelihood'])-set(target['likelihood']) == set(removed)
                        for name in target['likelihood']:
                            assert target['likelihood'][name] == native['likelihood'][name]
                        if evolution == 'linear':
                            assert target['params']['epsilon']['prior'] == {'min': -.5, 'max': .5}
                        assert target['params']['A_fg']['prior'] == {'min': 0., 'max': 2.}
                        config_cases += 1
                    if model == 'cpl':
                        bad = copy.deepcopy(frozen)
                        bad['configuration']['params']['w']['prior']['min'] = -2.9
                        must_fail(lambda: audit.target_configuration(settings, bad, 'sn', design), 'configuration changed')

    rng = np.random.default_rng(design['synthetic_seed'])
    n = design['synthetic_points']; groups = np.repeat(np.arange(4), n//4)
    factors = {'released_sn': (.3, 2.), 'bao.desi_dr2': (-.15, 3.),
               'act_dr6_lenslike.ACTDR6LensLike': (.2, 4.), 'SPT2023_lensing': (-.25, 4.)}
    parent_precision = 1+sum(1/s**2 for _, s in factors.values())
    parent_mean = sum(m/s**2 for m, s in factors.values())/parent_precision
    proposal_mean = parent_mean+.03; proposal_sd = parent_precision**-.5*1.02
    x = rng.normal(proposal_mean, proposal_sd, n); auxiliary = rng.normal(size=n)
    records = []
    for i in range(n):
        prior = float(norm.logpdf(x[i])+norm.logpdf(auxiliary[i]))
        exact = {name: float(norm.logpdf(x[i], *factors[name])) if name in factors else 0.
                 for name in sources}
        exact_post = prior+sum(exact.values())
        proposal_post = float(norm.logpdf(x[i], proposal_mean, proposal_sd)+norm.logpdf(auxiliary[i]))
        approximate = dict(exact)
        approximate['released_sn'] += proposal_post-exact_post
        records.append({'status': 'finite', 'point': {'w': float(x[i]), 'A_fg': float(auxiliary[i])},
                        'derived': {'q0': float(x[i]-.2), 'q05': float(.5*x[i]),
                                    'q1': float(auxiliary[i]), 'j0': 1.},
                        'exact_logpost': exact_post, 'proposal_logpost': proposal_post,
                        'log_weight': exact_post-proposal_post,
                        'exact_loglikes': exact, 'proposal_loglikes': approximate})
    tolerance = design['logweight_accounting_absolute_tolerance']
    recovery = []; all_rows = {}
    for omit, definition in design['alternatives'].items():
        removed = definition['removed_components']
        rows = [audit.remove_components(r, removed, sources, tolerance) for r in records]
        all_rows[omit] = rows
        result = audit.summarize_omission(records, rows, groups, gates)
        assert result['qualified_under_declared_numerical_gates'], result['failed_gates']
        retained = {k:v for k,v in factors.items() if k not in removed}
        precision = 1+sum(1/s**2 for _, s in retained.values())
        target_mean = sum(m/s**2 for m, s in retained.values())/precision
        expected_sd = precision**-.5
        actual = result['posterior']['w']
        assert abs(actual['mean']-target_mean) < design['synthetic_mean_absolute_tolerance']
        assert abs(actual['sd']-expected_sd) < design['synthetic_sd_absolute_tolerance']
        independent = np.array([norm.logpdf(v)+sum(norm.logpdf(v, *spec) for spec in retained.values())
                                -norm.logpdf(v, proposal_mean, proposal_sd) for v in x])
        computed = np.array([r['target_logweight'] for r in rows])
        # The unchanged independent auxiliary prior cancels exactly.
        density_error = float(np.max(abs(independent-computed)))
        assert density_error < 1e-11
        recovery.append({'omission': omit, 'synthetic_points': n, 'expected_mean': target_mean,
                         'weighted_mean': actual['mean'], 'expected_sd': expected_sd,
                         'weighted_sd': actual['sd'], 'maximum_independent_density_error': density_error,
                         'raw_weight_ESS': result['weight_diagnostics']['raw_weight_ESS'],
                         'Pareto_k': result['weight_diagnostics']['Pareto_k'],
                         'failed_gates': result['failed_gates']})

    null_rows = [audit.remove_components(r, [], sources, tolerance) for r in records]
    assert all(row['target_logweight'] == record['log_weight'] for row, record in zip(null_rows, records))
    null_result = audit.summarize_omission(records, null_rows, groups, gates)
    assert null_result['maximum_normalized_weight_change_from_parent'] == 0.
    original = np.array([r['target_logweight'] for r in all_rows['sn']])
    shifted = []
    for record in records[:256]:
        altered = copy.deepcopy(record)
        for key in ['exact_loglikes', 'proposal_loglikes']:
            altered[key]['released_sn'] += 500.
        altered['exact_logpost'] += 500.; altered['proposal_logpost'] += 500.
        shifted.append(audit.remove_components(altered, ['released_sn'], sources, tolerance)['target_logweight'])
    shifted = np.array(shifted)
    constant_error = float(np.max(abs(np.exp(original[:256]-logsumexp(original[:256]))-
                                     np.exp(shifted-logsumexp(shifted)))))
    assert constant_error < 1e-13
    pair_error = 0.
    for record, row in zip(records[:256], all_rows['lensing'][:256]):
        left = audit.remove_components(record, ['act_dr6_lenslike.ACTDR6LensLike'], sources, tolerance)
        difference = row['target_logweight']-(left['target_logweight']-record['exact_loglikes']['SPT2023_lensing'])
        pair_error = max(pair_error, abs(difference))
    assert pair_error < 1e-12
    must_fail(lambda: audit.remove_components(records[0], ['made_up_probe'], sources, tolerance))
    must_fail(lambda: audit.remove_components(records[0], ['released_sn']*2, sources, tolerance))
    changed = copy.deepcopy(records[0]); changed['proposal_logpost'] += .01; changed['log_weight'] -= .01
    must_fail(lambda: audit.remove_components(changed, ['released_sn'], sources, tolerance), 'prior accounting changed')
    collapsed = copy.deepcopy(all_rows['sn']); collapsed[0]['target_logweight'] += 1000.
    failure = audit.summarize_omission(records, collapsed, groups, gates)
    assert not failure['qualified_under_declared_numerical_gates']
    assert 'raw_weight_ESS' in failure['failed_gates']
    assert all(failure[key] is None for key in ['posterior','conditional_sign_fractions','weighted_covariance','paired_mean_changes_from_parent'])
    assert failure['weight_diagnostics']['per_chain']['1']['weight_fraction'] == 0.
    json.dumps(failure, allow_nan=False)
    short = audit.summarize_omission(records[:1000], all_rows['sn'][:1000], np.repeat(np.arange(4),250), gates)
    assert 'minimum_exact_points' in short['failed_gates'] and short['posterior'] is None

    with tempfile.TemporaryDirectory() as directory:
        cache = Path(directory)/'must-not-create'
        with patch.object(audit, 'summarize_run', side_effect=ValueError('unqualified-parent')):
            with patch.object(audit, 'configuration_factory', side_effect=AssertionError('configuration accessed too early')):
                must_fail(lambda: audit.actual(directory, Path(directory)/'absent', 'sn', cache), 'unqualified-parent')
        assert not cache.exists()

    # The real qualifier is independently validated elsewhere. This explicitly
    # synthetic consumer fixture mocks its return and configuration factory;
    # production has no bypass. Its unchanged input-byte contract is not mocked.
    with tempfile.TemporaryDirectory(dir=audit.ROOT/'.work') as directory:
        directory = Path(directory); folder = directory/'synthetic-chain'; native_dir = folder/'exact'
        native_dir.mkdir(parents=True)
        dataset = directory/'synthetic.npz'; dataset.write_bytes(b'synthetic-hash-only-no-observations')
        fixture_settings = dict(model='cpl',evolution='linear',sample='synthetic',calibration='official_planck',surrogate=str(fixture_model))
        factory = audit.configuration_factory(fixture_settings)
        base = factory(model='cpl',evolution='linear',sample='dovekie')
        base['likelihood']['released_sn']['data_file'] = str(dataset)
        proposal_config = factory(model='cpl',evolution='linear',sample='dovekie',surrogate=fixture_model)
        proposal_config['likelihood']['released_sn']['data_file'] = str(dataset)
        def synthetic_factory(**kwargs):
            return copy.deepcopy(proposal_config if kwargs.get('surrogate') else base)
        frozen = {'identity': 'synthetic-not-a-qualified-cosmology', 'configuration': audit.canonical(proposal_config),
                  'source_sha256': {audit.relative(audit.HERE/'likelihood.py'):audit.digest(audit.HERE/'likelihood.py')},
                  'versions': {name:importlib.metadata.version(name) for name in ['numpy','scipy']},
                  'assets': {'synthetic':'no-native-likelihood-inputs'}, 'sample_sha256':audit.digest(dataset)}
        manifest = folder/'run-0.json'; manifest.write_text(json.dumps({'target_identity': frozen}))
        pick = np.concatenate([np.arange(500)+i*(n//4) for i in range(4)])
        selected = [copy.deepcopy(records[i]) for i in pick]
        for record in selected:
            # The unit-width proper uniform prior adds zero to log density.
            # This fixture tests the auxiliary reporting contract, not an SN fit.
            record['point']['epsilon'] = float(rng.uniform(-.5, .5))
        selection = native_dir/'selection.json'
        selection.write_text(json.dumps({'settings':fixture_settings,'groups':np.repeat(np.arange(4),500).tolist(),
                                        'points':[r['point'] for r in selected]}))
        summary = directory/'synthetic-summary.json'; summary.write_text(json.dumps({'selection_path':audit.relative(selection)}))
        input_paths = [manifest,selection,summary]
        for index, record in enumerate(selected):
            path=native_dir/f'{index:05d}.json'; path.write_text(json.dumps(record)); input_paths.append(path)
        parent = {'settings':{k:fixture_settings[k] for k in ['model','evolution','sample','calibration']},
                  'input_sha256':{audit.relative(p):audit.digest(p) for p in input_paths}}
        with patch.object(audit,'summarize_run',return_value=parent):
            with patch.object(audit,'configuration_factory',return_value=synthetic_factory):
                cache=directory/'cache'
                first=audit.actual(folder,summary,'sn',cache)
                assert first['unused_auxiliary_coordinates']==['epsilon']
                assert first['qualified_under_declared_numerical_gates'], first['failed_gates']
                assert 'epsilon' in first['retained_sampled_priors']
                assert 'epsilon' in first['weight_diagnostics']['weighted_chain_stability']
                assert 'epsilon' in first['weight_diagnostics']['weighted_summaries_for_diagnostics']
                assert 'epsilon' not in first['posterior']
                assert 'epsilon' not in first['weighted_covariance']['parameter_order']
                assert 'epsilon' not in first['paired_mean_changes_from_parent']
                second=audit.actual(folder,summary,'sn',cache);assert first==second
                rowpath=cache/'00000.json';old=rowpath.read_bytes();row=json.loads(old);row['target_logweight']+=.1
                rowpath.write_text(json.dumps(row))
                must_fail(lambda:audit.actual(folder,summary,'sn',cache),'identity changed')
                rowpath.write_bytes(old)
                old=manifest.read_bytes();manifest.write_bytes(old+b' ')
                must_fail(lambda:audit.actual(folder,summary,'sn',directory/'wrong-parent'),'identity changed')
                assert not (directory/'wrong-parent').exists();manifest.write_bytes(old)
                original_data=dataset.read_bytes();dataset.write_bytes(b'changed')
                must_fail(lambda:audit.actual(folder,summary,'sn',directory/'wrong-data'),'Parent SN data changed')
                assert not (directory/'wrong-data').exists();dataset.write_bytes(original_data)
                must_fail(lambda:audit.actual(folder,summary,'lensing',cache),'payload changed')
                fresh=audit.actual(folder,summary,'lensing',directory/'lensing-cache')
                assert fresh['unused_auxiliary_coordinates']==['A_fg']
    paths=[Path(__file__),Path(audit.__file__),audit.DESIGN,audit.GATES,
           audit.HERE/'measurement_summary.py',audit.HERE/'exact_correction.py',
           audit.HERE/'luminosity_sensitivity.py',audit.HERE/'target_identity.py',
           audit.HERE/'modern_run.py',audit.HERE/'modern_fast.py',audit.HERE/'modern_gpu.py',
           audit.HERE.parent/'external_probes/modern_adapter.py']
    return {'status':'passed_synthetic_probe_omission_validation','observational_points_used':0,
            'CMB_spectrum_calls':0,'background_calls':0,'GPU_calls':0,
            'configuration_cases':config_cases,'known_Gaussian_recoveries':recovery,
            'null_removal_preserves_raw_weights_exactly':True,
            'component_constant_normalized_weight_error':constant_error,
            'two_lensing_component_accounting_error':pair_error,
            'unknown_or_duplicate_components_rejected':True,'changed_prior_accounting_rejected':True,
            'insufficient_points_and_concentrated_weights_withhold_measurements':True,
            'unqualified_parent_rejected_before_cache_or_configuration':True,
            'consumer_fixture_points':2000,'consumer_qualification_and_configuration_explicitly_mocked':True,
            'consumer_replay_identical':True,'changed_parent_and_cache_and_data_rejected':True,
            'wrong_target_cache_rejected':True,'unused_epsilon_and_A_fg_retained':True,
            'unused_epsilon_checked_but_not_published_as_measurement':True,
            'source_sha256':{audit.relative(p):audit.digest(p) for p in paths},
            'scope':'Synthetic mathematics and consumer-contract validation only. No reduced-data cosmological target has been evaluated.'}


def main():
    with patch('camb.get_results',side_effect=AssertionError('Forbidden native spectra call')) as spectra:
        with patch('camb.get_background',side_effect=AssertionError('Forbidden background call')) as background:
            with patch('cobaya.model.get_model',side_effect=AssertionError('Forbidden model initialization')) as model:
                report=run_tests()
                assert spectra.call_count==background.call_count==model.call_count==0
    output=audit.ROOT/'studies/unified_cosmology/results/inference/probe-omission-validation.json'
    output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':report['status'],'configuration_cases':report['configuration_cases'],
                      'known_Gaussian_recoveries':report['known_Gaussian_recoveries']},indent=2))


if __name__=='__main__':
    main()
