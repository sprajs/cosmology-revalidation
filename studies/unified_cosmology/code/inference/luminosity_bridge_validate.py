"""Synthetic exact-density/weight tests for the separately identified bridge."""
import hashlib
import json
import copy
from pathlib import Path
import numpy as np

from luminosity_bridge import (ROOT, HERE, DESIGN, GATES, bridge_record,
    summarize_bridge, identity, configuration_for_settings, backend_source_paths)
from luminosity_sensitivity import IntegratedLuminosity
from target_identity import canonical


def main():
    configurations = 0
    proposal_replays = 0
    fixture = ROOT/'.work/unified-cosmology/inference/surrogate/cubic-0509.npz'
    backends = [({}, 'modern_run'), ({'fast_lensing': True}, 'modern_fast'),
                ({'gpu': True, 'fast_lensing': True}, 'modern_gpu')]
    from modern_run import configuration as original_configuration
    from modern_fast import configuration as fast_configuration
    from modern_gpu import configuration as gpu_configuration
    independent_factories = dict(modern_run=original_configuration,
                                 modern_fast=fast_configuration, modern_gpu=gpu_configuration)
    for settings, module in backends:
        configuration = configuration_for_settings(settings)
        assert configuration is independent_factories[module]
        assert (HERE/'modern_gpu.py' in backend_source_paths(settings)) == bool(settings.get('gpu'))
        for model in ['lcdm', 'cpl']:
            for calibration in ['official_planck', 'paper_literal']:
                for source in ['none', 'smooth01', 'smooth03']:
                    options = dict(model=model, evolution=source, calibration=calibration)
                    before = canonical(configuration(**options))
                    # Rebuild exactly the factory/class recorded by a CPU or GPU
                    # parent without constructing a model or touching the GPU.
                    proposal = canonical(configuration(**options, surrogate=fixture))
                    frozen = canonical(independent_factories[module](**options, surrogate=fixture))
                    assert proposal == frozen
                    actual_class = proposal['theory']['spectral_surrogate']['external']
                    assert ('GPUSpectralSurrogate' in str(actual_class)) == bool(settings.get('gpu'))
                    proposal_replays += 1
                    for target, sigma in [('smooth01', .1), ('smooth03', .3)]:
                        after = canonical(configuration(model=model, evolution=target, calibration=calibration))
                        changed = copy.deepcopy(before)
                        changed['likelihood']['released_sn']['smooth_sigma'] = sigma
                        assert changed == after
                        configurations += 1
    rng = np.random.default_rng(272813)
    sigma = {'none': 0., 'smooth01': .1, 'smooth03': .3}
    rows = []
    for case in range(16):
        n = 24+case
        z = np.sort(rng.uniform(.01, 1.15, n))
        factor = rng.normal(size=(n, n))
        covariance = (factor@factor.T/n+np.eye(n))*.1**2
        residual = rng.normal(size=n)*.2
        sn = IntegratedLuminosity(z, z, -residual, covariance)
        background = sn.from_prediction(np.zeros(n))
        direct = {}
        for name, width in sigma.items():
            c = covariance+width**2*sn.basis@sn.basis.T
            p = np.linalg.inv(c); u = p@np.ones(n)
            score = residual@p@residual-(u@residual)**2/u.sum()
            normalization = np.linalg.slogdet(c)[1]+np.log(u.sum())+(n-1)*np.log(2*np.pi)
            direct[name] = -.5*(score+normalization)
        for source in sigma:
            for target in ['smooth01', 'smooth03']:
                parent = {'log_weight': -.127, 'exact_loglikes': {'released_sn': direct[source]}}
                row = bridge_record(parent, background, source, target, 1e-5)
                assert row['status'] == 'finite_bridge'
                error = abs(row['SN_log_likelihood_ratio']-(direct[target]-direct[source]))
                assert error < 1e-9
                assert abs(row['target_logweight']-parent['log_weight']-row['SN_log_likelihood_ratio']) < 1e-14
                if source == target:
                    assert row['SN_log_likelihood_ratio'] == 0.
                    assert row['target_logweight'] == parent['log_weight']
                rows.append({'case': case, 'source': source, 'target': target,
                             'independent_augmented_covariance_log_ratio_error': float(error),
                             'source_density_closure': row['source_SN_loglike_closure']})
    # A full four-chain known-Gaussian toy checks the unchanged-target path and
    # deliberately inadequate overlap, without any observational/native points.
    n = 8000
    x = rng.normal(size=n); y = rng.normal(size=n)
    groups = np.repeat(np.arange(4), n//4)
    parent = [{'point': {'x': float(a), 'y': float(b)},
               'derived': {'q0': float(a), 'q05': float(.5*a+b), 'q1': float(b), 'j0': 1.},
               'log_weight': float(.15*a-.5*.15**2),
               'exact_loglikes': {'released_sn': -12.}}
              for a, b in zip(x, y)]
    background = {'baseline_SN_loglike': -13.,
                  'log_likelihood_ratios': {'smooth01': 1., 'smooth03': 1.}}
    same = [bridge_record(p, background, 'smooth01', 'smooth01', 1e-5) for p in parent]
    gates = json.loads(GATES.read_text())['overlap_gates']
    result = summarize_bridge(parent, same, groups, gates)
    assert result['status'] == 'qualified_conditional_bridge'
    assert result['maximum_normalized_weight_change_from_parent'] == 0.
    assert not result['failed_gates']
    drifted = []
    for index, row in enumerate(same):
        extra = 1000. if index == 0 else 0.
        drifted.append(dict(row, target_logweight=row['target_logweight']+extra,
                            SN_log_likelihood_ratio=extra))
    failed = summarize_bridge(parent, drifted, groups, gates)
    assert failed['status'] == 'insufficient_bridge_overlap_or_stability'
    assert failed['posterior'] is None and failed['conditional_sign_fractions'] is None
    assert 'raw_weight_ESS' in failed['failed_gates']
    assert failed['weight_diagnostics']['per_chain']['1']['weight_fraction'] == 0.
    json.dumps(failed, allow_nan=False)
    mismatch = bridge_record(parent[0], dict(background, baseline_SN_loglike=-12.), 'smooth01', 'smooth03', 1e-5)
    assert mismatch['status'] == 'source_SN_density_mismatch'
    bad = summarize_bridge(parent[:1], [mismatch], groups[:1], gates)
    assert bad['status'] == 'failed_bridge_evaluation' and bad['posterior'] is None
    # Cache checksum covers numerical contents as well as identities.
    old = {'identity': 'synthetic', 'history': [1., 2.]}
    changed = dict(old, history=[1., 3.])
    assert identity(old) != identity(changed)
    sources = [Path(__file__), HERE/'luminosity_bridge.py', DESIGN, GATES,
               HERE/'luminosity_sensitivity.py', HERE/'measurement_summary.py',
               HERE/'exact_correction.py', HERE/'modern_fast.py', HERE/'modern_run.py', HERE/'modern_gpu.py',
               HERE.parent/'external_probes/modern_adapter.py']
    report = {'status': 'passed_synthetic_bridge_validation', 'observational_points_used': 0,
              'CMB_spectrum_calls': 0, 'independent_augmented_covariance_comparisons': len(rows),
              'identical_non_SN_factor_configuration_comparisons': configurations,
              'configuration_backends': [module for _, module in backends],
              'independent_parent_proposal_configuration_replays': proposal_replays,
              'GPU_backend_source_hash_bound_when_used': True,
              'GPU_calls': 0,
              'maximum_log_ratio_error': max(r['independent_augmented_covariance_log_ratio_error'] for r in rows),
              'maximum_source_density_closure': max(abs(r['source_density_closure']) for r in rows),
              'identical_target_preserves_original_raw_weights_exactly': True,
              'qualified_identical_target_synthetic_points': n,
              'concentrated_target_weights_fail_and_preserve_zero_chain_mass': True,
              'source_density_mismatch_withholds_every_posterior': True,
              'checksum_detects_changed_numerical_payload': True,
              'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
              'scope': 'Synthetic density-ratio and weighted-gate checks, not measured overlap for any cosmological luminosity prior.'}
    path = ROOT/'studies/unified_cosmology/results/inference/luminosity-bridge-validation.json'
    path.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
