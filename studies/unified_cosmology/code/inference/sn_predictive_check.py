"""Conditional marginal-Gaussian SN quadratic check, with no new CAMB calls."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.linalg import cholesky, solve_triangular
from scipy.special import gammainc, gammaincc, logsumexp

from measurement_summary import summarize_run
from exact_correction import weighted_summary

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DESIGN = HERE/'sn-predictive-design.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


class ProjectedGaussian:
    """Keep the full covariance and remove just the unpenalized intercept."""
    def __init__(self, covariance, z, sigma, design):
        covariance = np.asarray(covariance, dtype=float)
        z = np.asarray(z, dtype=float)
        self.n = len(z)
        assert self.n > 1 and covariance.shape == (self.n, self.n)
        assert np.isfinite(covariance).all() and np.isfinite(z).all()
        assert np.all(z > 0) and sigma in set(design['allowed_evolution'].values())
        assert np.allclose(covariance, covariance.T,
                           rtol=design['covariance_symmetry_relative_tolerance'],
                           atol=design['covariance_symmetry_absolute_tolerance']), 'Asymmetric released covariance.'
        self.basis = CubicSpline(design['smooth_knots'], np.eye(5)[:, 1:],
                                 bc_type='natural')(z)
        self.covariance = covariance+sigma**2*self.basis@self.basis.T
        # No symmetrization, diagonal approximation, inverse-error reweighting or
        # empirical covariance rescaling is inserted into the released target.
        self.factor = cholesky(self.covariance, lower=True)
        constant = solve_triangular(self.factor, np.ones(self.n), lower=True)
        self.intercept_information = float(constant@constant)
        self.unit_constant = constant/np.sqrt(self.intercept_information)
        self.logdet = float(2*np.log(np.diag(self.factor)).sum())
        self.df = self.n-1
        self.log_normalization = float(self.logdet+np.log(self.intercept_information)
                                       +self.df*np.log(2*np.pi))

    def quadratic(self, residual):
        """Used for synthetic checks; actual diagnostics reuse native sn_chi2."""
        residual = np.asarray(residual, dtype=float)
        assert residual.shape[0] == self.n
        centered = residual-residual.mean(axis=0)
        white = solve_triangular(self.factor, centered, lower=True)
        if white.ndim == 1:
            projected = white-self.unit_constant*(self.unit_constant@white)
        else:
            projected = white-np.outer(self.unit_constant, self.unit_constant@white)
        return np.sum(projected**2, axis=0)

    def loglike_from_quadratic(self, quadratic):
        return -.5*(np.asarray(quadratic)+self.log_normalization)


def weighted_average(values, weights, groups, counts):
    """Chronological ratio-estimator batch uncertainty; no weight smoothing."""
    values = np.asarray(values); groups = np.asarray(groups)
    average = float(weights@values)
    per_chain = {}
    batches = {}
    for group in np.unique(groups):
        indices = np.flatnonzero(groups == group)
        mass = float(weights[indices].sum())
        assert mass > 0
        per_chain[str(group)] = {'mean': float(weights[indices]@values[indices]/mass),
                                 'weight_mass': mass}
    for count in counts:
        scores = []; masses = []; sizes = []
        for group in np.unique(groups):
            for indices in np.array_split(np.flatnonzero(groups == group), count):
                assert len(indices) >= 10, 'Too few selected points in a diagnostic batch.'
                scores.append(float(weights[indices]@(values[indices]-average)))
                masses.append(float(weights[indices].sum())); sizes.append(len(indices))
        batches[str(count)] = {'mean_mcse': float(np.sqrt(np.var(scores, ddof=1)/len(scores))/np.mean(masses)),
                               'minimum_points_per_batch': min(sizes)}
    return {'mean': average, 'independent_chain_means': per_chain,
            'chronological_batch_mcse': batches,
            'interpretation': 'Monte Carlo precision diagnostic for this weighted expectation; not a coverage guarantee.'}


def check_records(records, groups, gaussian, design):
    groups = np.asarray(groups)
    assert len(records) == len(groups) and set(groups) == {0, 1, 2, 3}
    quadratics = []; logweights = []; closures = []; failures = []
    for index, record in enumerate(records):
        try:
            q = float(record['derived']['sn_chi2'])
            stored_loglike = float(record['exact_loglikes']['released_sn'])
            logweight = float(record['log_weight'])
            assert np.isfinite([q, stored_loglike, logweight]).all() and q >= 0
            closure = float(gaussian.loglike_from_quadratic(q)-stored_loglike)
            assert abs(closure) <= design['normalized_SN_loglike_closure_absolute_tolerance'], 'Native normalized SN density mismatch.'
            quadratics.append(q); logweights.append(logweight); closures.append(closure)
        except (KeyError, ValueError, TypeError, AssertionError) as error:
            failures.append({'index': index, 'error': str(error)})
    if failures:
        return {'status': 'failed_native_SN_density_check', 'failures': failures,
                'posterior_tail_diagnostics': None}, None
    q = np.asarray(quadratics); lw = np.asarray(logweights)
    weights = np.exp(lw-logsumexp(lw))
    lower = gammainc(gaussian.df/2., q/2.)
    upper = gammaincc(gaussian.df/2., q/2.)
    assert np.isfinite(lower).all() and np.isfinite(upper).all()
    assert max(abs(lower+upper-1.)) <= design['tail_complement_absolute_tolerance']
    result = {'status': 'computed_conditional_marginal_Gaussian_quadratic_check',
              'observations': gaussian.n, 'replicate_degrees_of_freedom': gaussian.df,
              'subtracted_unpenalized_intercepts': 1,
              'subtracted_cosmological_or_proper_Gaussian_mode_counts': 0,
              'posterior_tail_diagnostics': {
                  'replicate_quadratic_below_observed': weighted_average(lower, weights, groups, design['batch_counts']),
                  'replicate_quadratic_above_observed': weighted_average(upper, weights, groups, design['batch_counts'])},
              'observed_quadratic': weighted_summary(q, weights),
              'observed_quadratic_per_replicate_df': weighted_summary(q/gaussian.df, weights),
              'maximum_native_normalized_SN_loglike_closure': max(abs(np.asarray(closures))),
              'tail_complement_maximum_error': float(max(abs(lower+upper-1.))),
              'floating_point_zero_tail_counts': {'lower': int(np.sum(lower == 0.)), 'upper': int(np.sum(upper == 0.))},
              'original_raw_weight_ESS': float(1/(weights@weights)),
              'largest_original_normalized_weight': float(weights.max()),
              'effective_covariance_logdet': gaussian.logdet,
              'free_magnitude_information': gaussian.intercept_information,
              'normalized_density_constant': gaussian.log_normalization}
    replay = {'index': list(range(len(records))), 'groups': groups.tolist(),
              'stored_native_quadratic': q.tolist(), 'original_native_proposal_logweight': lw.tolist(),
              'untrimmed_normalized_weight': weights.tolist(), 'normalized_SN_density_closure': closures,
              'conditional_lower_tail': lower.tolist(), 'conditional_upper_tail': upper.tolist()}
    return result, replay


def actual(folder, summary_path, cache):
    # First gate: do not load even the released data for an unqualified parent.
    parent = summarize_run(folder, summary_path)
    folder = Path(folder).resolve(); summary_path = Path(summary_path).resolve()
    cache = Path(cache).resolve(); assert cache.is_relative_to(ROOT/'.work')
    summary = json.loads(summary_path.read_text())
    selection_path = ROOT/summary['selection_path']
    selection = json.loads(selection_path.read_text())
    assert selection_path.parent.parent == folder
    settings = selection['settings']
    frozen = json.loads((folder/'run-0.json').read_text())['target_identity']
    design = json.loads(DESIGN.read_text())
    evolution = settings['evolution']; sigma = design['allowed_evolution'][evolution]
    sn_configuration = frozen['configuration']['likelihood']['released_sn']
    assert sn_configuration['external'] == 'likelihood.ReleasedDistances'
    assert sn_configuration['smooth_sigma'] == sigma
    expected_epsilon = ({'min': -.5, 'max': .5} if evolution == 'linear' else 0.)
    epsilon = frozen['configuration']['params']['epsilon']
    assert (epsilon['prior'] if evolution == 'linear' else epsilon) == expected_epsilon
    data_path = Path(sn_configuration['data_file']).resolve()
    assert data_path.is_relative_to(ROOT)
    assert digest(data_path) == frozen['sample_sha256'], 'Released SN data changed.'
    likelihood_source = HERE/'likelihood.py'
    assert digest(likelihood_source) == frozen['source_sha256'][relative(likelihood_source)], 'Released SN model source changed.'
    for name in ['numpy', 'scipy']:
        assert importlib.metadata.version(name) == frozen['versions'][name]
    with np.load(data_path, allow_pickle=False) as data:
        z = np.asarray(data['zHD']); zhel = np.asarray(data['zHEL']); mu = np.asarray(data['MU'])
        assert z.shape == zhel.shape == mu.shape and np.isfinite(mu).all() and np.isfinite(zhel).all()
        gaussian = ProjectedGaussian(data['covariance'], z, sigma, design)
    records = [json.loads((selection_path.parent/f'{index:05d}.json').read_text())
               for index in range(len(selection['points']))]
    # summarize_run already validates each record's seal, point, target and file
    # hash, the original importance weights, and all numerical qualification gates.
    result, replay = check_records(records, selection['groups'], gaussian, design)
    sources = [Path(__file__), DESIGN, HERE/'measurement_summary.py',
               HERE/'exact_correction.py', likelihood_source]
    hashes = {relative(path): digest(path) for path in sources}
    lineage = {'qualified_parent_inputs': parent['input_sha256'], 'source_sha256': hashes,
               'parent_target_identity': frozen['identity'], 'released_SN_configuration': sn_configuration,
               'SN_data_path': relative(data_path), 'SN_data_sha256': digest(data_path),
               'numerical_versions': {name: importlib.metadata.version(name) for name in ['numpy', 'scipy']}}
    cache.mkdir(parents=True, exist_ok=True)
    lineage_path = cache/'lineage.json'
    if lineage_path.exists():
        assert json.loads(lineage_path.read_text()) == lineage, 'Different diagnostic cache identity.'
    else:
        lineage_path.write_text(json.dumps(lineage, indent=2)+'\n')
    manifest_path = cache/'record-hashes.json'
    if manifest_path.exists():
        for name, expected in json.loads(manifest_path.read_text()).items():
            assert digest(ROOT/name) == expected, 'Previously recorded diagnostic payload changed.'
    replay_path = cache/'quadratic-points.json'
    if replay is not None:
        if replay_path.exists():
            assert json.loads(replay_path.read_text()) == replay, 'Diagnostic replay disagrees with exact parent records.'
        else:
            replay_path.write_text(json.dumps(replay, indent=2, allow_nan=False)+'\n')
    output_hashes = {relative(lineage_path): digest(lineage_path)}
    if replay_path.exists():
        output_hashes[relative(replay_path)] = digest(replay_path)
    manifest_path.write_text(json.dumps(output_hashes, indent=2)+'\n')
    result.update(settings=parent['settings'], sigma_mag=sigma,
                  parent_target_identity=frozen['identity'],
                  parent_qualified_under_declared_numerical_gates=True,
                  conditional_replication=('Fresh Gaussian luminosity coefficients drawn from their declared prior for each dataset; not a shared-curve posterior prediction.'
                                           if sigma else 'Projected Gaussian distance noise at fixed cosmological and sampled nuisance parameters.'),
                  scope=design['scope'], limitations=design['limitations'],
                  flat_magnitude_measure=design['flat_magnitude_measure'],
                  covariance_adjusted_or_rescaled=False, CMB_spectrum_calls=0, background_calls=0,
                  source_sha256=hashes, lineage_path=relative(lineage_path), lineage_sha256=digest(lineage_path),
                  record_manifest_path=relative(manifest_path), record_manifest_sha256=digest(manifest_path),
                  correction_summary_path=relative(summary_path), correction_summary_sha256=digest(summary_path))
    for mapping in [hashes, parent['input_sha256'], {relative(data_path): frozen['sample_sha256']}]:
        assert all(digest(ROOT/name) == expected for name, expected in mapping.items()), 'Inputs changed during diagnostic.'
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
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': result['status'], 'CMB_spectrum_calls': 0, 'background_calls': 0}))


if __name__ == '__main__':
    main()
