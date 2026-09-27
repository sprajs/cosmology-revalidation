"""Independent Gaussian-replication and release-density checks; no CAMB calls."""
import hashlib
import importlib.metadata
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from scipy.linalg import cholesky, qr
from scipy.special import gammainc

import sn_predictive_check as audit
from likelihood import ReleasedDistances


def main():
    design = json.loads(audit.DESIGN.read_text())
    rng = np.random.default_rng(design['synthetic_seed'])
    count = design['synthetic_replicates_per_configuration']
    sigma_limit = design['synthetic_moment_and_tail_sigma_tolerance']
    rows = []
    for n in design['synthetic_sample_sizes']:
        z = np.sort(rng.uniform(.01, 1.2, n))
        mixing = rng.normal(size=(n, n))
        covariance = .01*(mixing@mixing.T/n+np.eye(n))
        noise_factor = cholesky(covariance, lower=True)
        for sigma in design['synthetic_sigma_mag']:
            gaussian = audit.ProjectedGaussian(covariance, z, sigma, design)
            # Generate measurement noise plus a single shared coefficient vector
            # per replicate, rather than drawing directly from effective Sigma.
            replicated = noise_factor@rng.normal(size=(n, count))
            replicated += sigma*gaussian.basis@rng.normal(size=(4, count))
            values = gaussian.quadratic(replicated)
            df = n-1
            mean_z = abs(values.mean()-df)/np.sqrt(2*df/count)
            variance_z = abs(values.var(ddof=1)-2*df)/np.sqrt(8*df*(df+6)/count)
            assert mean_z < sigma_limit and variance_z < sigma_limit
            probabilities = gammainc(df/2., values/2.)
            tail_checks = {}
            for level in [.025, .05, .5, .95, .975]:
                fraction = float(np.mean(probabilities <= level))
                score = abs(fraction-level)/np.sqrt(level*(1-level)/count)
                assert score < sigma_limit
                tail_checks[str(level)] = {'observed_fraction': fraction, 'standard_error_units': score}
            assert abs(probabilities.mean()-.5) < sigma_limit/np.sqrt(12*count)
            subset = replicated[:, :256]
            offset_error = float(max(abs(gaussian.quadratic(subset+53.)-values[:256])))
            assert offset_error < 1e-9
            # Independent dense GLS inverse and QR contrast basis both recover
            # the same N-1-dimensional statistic and determinant normalization.
            precision = np.linalg.inv(gaussian.covariance)
            u = precision@np.ones(n)
            direct = np.sum(subset*(precision@subset), axis=0)-(u@subset)**2/u.sum()
            direct_error = float(max(abs(direct-values[:256])))
            assert direct_error < 1e-9
            white = np.linalg.solve(gaussian.factor, subset)
            const = np.linalg.solve(gaussian.factor, np.ones(n))
            orthogonal, _ = qr(const[:, None], mode='full')
            contrasts = orthogonal[:, 1:].T@white
            qr_error = float(max(abs(np.sum(contrasts**2, axis=0)-values[:256])))
            assert qr_error < 1e-9
            independent_norm = np.linalg.slogdet(gaussian.covariance)[1]+np.log(u.sum())+df*np.log(2*np.pi)
            norm_error = abs(independent_norm-gaussian.log_normalization)
            assert norm_error < 1e-10
            # Execute the first-party released SN density on a synthetic distance
            # provider. There is no cosmological background or spectrum call.
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory)/'synthetic.npz'
                observed = 39.-subset[:, 0]
                np.savez(path, zHD=z, zHEL=z, MU=observed, covariance=covariance)
                likelihood = object.__new__(ReleasedDistances)
                likelihood.data_file = str(path); likelihood.smooth_sigma = sigma
                likelihood.initialize()
                distance = np.full(n, 10**((39.-25.)/5))/(1+z)**2
                likelihood.provider = SimpleNamespace(get_angular_diameter_distance=lambda _: distance)
                derived = {}; loglike = likelihood.logp(epsilon=0., _derived=derived)
                source_score_error = abs(derived['sn_chi2']-values[0])
                source_density_error = abs(loglike-gaussian.loglike_from_quadratic(values[0]))
                assert source_score_error < 1e-9 and source_density_error < 1e-9
            rows.append({'n': n, 'sigma_mag': sigma, 'replicates': count, 'replicate_df': df,
                         'mean': float(values.mean()), 'variance': float(values.var(ddof=1)),
                         'mean_standard_error_units': float(mean_z), 'variance_standard_error_units': float(variance_z),
                         'nominal_tail_checks': tail_checks, 'offset_quadratic_maximum_error': offset_error,
                         'independent_GLS_maximum_error': direct_error, 'independent_QR_maximum_error': qr_error,
                         'normalization_error': float(norm_error), 'released_density_error': float(source_density_error),
                         'released_quadratic_error': float(source_score_error)})
    # Numerical weighting/accounting fixtures do not masquerade as cosmological
    # posterior data. All deliberately corrupt records must withhold the output.
    total = 4000; groups = np.repeat(np.arange(4), total//4)
    q = rng.chisquare(gaussian.df, size=total)
    lw = rng.normal(0., .15, total)
    records = [{'derived': {'sn_chi2': float(v)},
                'exact_loglikes': {'released_sn': float(gaussian.loglike_from_quadratic(v))},
                'log_weight': float(w)} for v, w in zip(q, lw)]
    result, replay = audit.check_records(records, groups, gaussian, design)
    weights = np.exp(lw-lw.max()); weights /= weights.sum()
    direct_lower = float(weights@gammainc(gaussian.df/2., q/2.))
    actual_lower = result['posterior_tail_diagnostics']['replicate_quadratic_below_observed']['mean']
    assert abs(direct_lower-actual_lower) < 1e-14
    assert abs(result['posterior_tail_diagnostics']['replicate_quadratic_above_observed']['mean']+actual_lower-1.) < 1e-14
    assert result['replicate_degrees_of_freedom'] == gaussian.n-1
    changed = json.loads(json.dumps(records)); changed[3]['derived']['sn_chi2'] += .01
    failure, unavailable = audit.check_records(changed, groups, gaussian, design)
    assert failure['status'] == 'failed_native_SN_density_check'
    assert failure['posterior_tail_diagnostics'] is None and unavailable is None
    assert failure['failures'][0]['index'] == 3
    changed = json.loads(json.dumps(records)); changed[4]['derived']['sn_chi2'] = -1.
    failure, _ = audit.check_records(changed, groups, gaussian, design)
    assert failure['status'] == 'failed_native_SN_density_check'
    with tempfile.TemporaryDirectory() as directory:
        destination = Path(directory)/'must-not-create'
        with patch.object(audit, 'summarize_run', side_effect=ValueError('unqualified-parent')):
            with patch.object(audit.np, 'load', side_effect=AssertionError('Data read before qualification')):
                try:
                    audit.actual(directory, Path(directory)/'missing.json', destination)
                except ValueError as error:
                    assert str(error) == 'unqualified-parent'
                else:
                    raise AssertionError('Unqualified parent was accepted.')
        assert not destination.exists()
    # Exercise the consumer's data/configuration/lineage/cache contract on an
    # explicitly synthetic parent. Only the already-tested qualification
    # dependency is mocked; production has no bypass option.
    with tempfile.TemporaryDirectory(dir=audit.ROOT/'.work') as directory:
        directory = Path(directory); folder = directory/'synthetic-chain'
        native = folder/'synthetic-correction'; native.mkdir(parents=True)
        dataset = directory/'synthetic-data.npz'
        np.savez(dataset, zHD=z, zHEL=z, MU=np.full(len(z), 39.), covariance=covariance)
        settings = {'model': 'cpl', 'evolution': 'smooth03', 'sample': 'synthetic',
                    'calibration': 'synthetic'}
        frozen = {'identity': 'synthetic-not-a-qualified-cosmology',
                  'configuration': {'params': {'epsilon': 0.}, 'likelihood': {'released_sn': {
                      'external': 'likelihood.ReleasedDistances', 'smooth_sigma': .3, 'data_file': str(dataset)}}},
                  'sample_sha256': audit.digest(dataset),
                  'source_sha256': {audit.relative(audit.HERE/'likelihood.py'): audit.digest(audit.HERE/'likelihood.py')},
                  'versions': {name: importlib.metadata.version(name) for name in ['numpy', 'scipy']}}
        manifest = folder/'run-0.json'; manifest.write_text(json.dumps({'target_identity': frozen}))
        locations = np.repeat(np.arange(4), 500).tolist()
        selection = native/'selection.json'
        selection.write_text(json.dumps({'settings': settings, 'groups': locations, 'points': [{}]*2000}))
        summary = directory/'synthetic-summary.json'
        summary.write_text(json.dumps({'selection_path': audit.relative(selection)}))
        input_paths = [manifest, selection, summary]
        for index, record in enumerate(records[:2000]):
            path = native/f'{index:05d}.json'; path.write_text(json.dumps(record)); input_paths.append(path)
        parent_fixture = {'settings': settings, 'input_sha256': {audit.relative(p): audit.digest(p) for p in input_paths}}
        cache = directory/'cache'
        with patch.object(audit, 'summarize_run', return_value=parent_fixture):
            first = audit.actual(folder, summary, cache)
            second = audit.actual(folder, summary, cache)
            assert first == second
            replay_file = cache/'quadratic-points.json'
            before = replay_file.read_bytes()
            altered = json.loads(before); altered['stored_native_quadratic'][0] += .1
            replay_file.write_text(json.dumps(altered))
            try:
                audit.actual(folder, summary, cache)
            except AssertionError as error:
                assert 'Previously recorded diagnostic payload changed' in str(error)
            else:
                raise AssertionError('Changed cached quadratic accepted.')
            replay_file.write_bytes(before)
            original_data = dataset.read_bytes()
            np.savez(dataset, zHD=z, zHEL=z, MU=np.full(len(z), 39.), covariance=2*covariance)
            try:
                audit.actual(folder, summary, cache)
            except AssertionError as error:
                assert 'Released SN data changed' in str(error)
            else:
                raise AssertionError('Different covariance data accepted.')
            dataset.write_bytes(original_data)
            changed = json.loads(manifest.read_text()); changed['target_identity']['configuration']['likelihood']['released_sn']['smooth_sigma'] = .1
            manifest.write_text(json.dumps(changed))
            try:
                audit.actual(folder, summary, directory/'must-not-create')
            except AssertionError:
                pass
            else:
                raise AssertionError('Wrong evolution/covariance setting accepted.')
            assert not (directory/'must-not-create').exists()
    sources = [Path(__file__), Path(audit.__file__), audit.DESIGN, audit.HERE/'likelihood.py',
               audit.HERE/'measurement_summary.py', audit.HERE/'exact_correction.py']
    report = {'status': 'passed_synthetic_conditional_Gaussian_validation',
              'observational_points_used': 0, 'CMB_spectrum_calls': 0, 'background_calls': 0,
              'configurations': rows, 'synthetic_replicates': count*len(rows),
              'weighted_tail_fixture_points': total, 'weighted_lower_tail_error': abs(direct_lower-actual_lower),
              'normalized_density_mismatch_withholds_all_tails': True,
              'negative_quadratic_rejected': True, 'unqualified_parent_rejected_before_data_or_cache': True,
              'synthetic_consumer_fixture_points': 2000,
              'synthetic_qualification_dependency_mocked_for_consumer_contract_only': True,
              'consumer_replay_identical': True, 'changed_cached_quadratic_rejected': True,
              'changed_released_covariance_data_rejected': True, 'mismatched_evolution_configuration_rejected': True,
              'source_sha256': {audit.relative(p): audit.digest(p) for p in sources},
              'scope': 'Synthetic model-calibration and algebra checks. No qualified cosmological posterior has been evaluated.'}
    path = audit.ROOT/'studies/unified_cosmology/results/inference/sn-predictive-validation.json'
    path.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': report['status'], 'replicates': report['synthetic_replicates'],
                      'max_offset_error': max(r['offset_quadratic_maximum_error'] for r in rows),
                      'max_released_density_error': max(r['released_density_error'] for r in rows)}))


if __name__ == '__main__':
    main()
