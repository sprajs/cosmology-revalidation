"""Configuration/provider checks with native theory and model construction forbidden."""
from copy import deepcopy
import argparse
import ast
import hashlib
import importlib.metadata
import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from scipy import linalg
import anchored_adapter as adapter


def validate(surrogate):
    import camb
    import cobaya.model
    from spectral_surrogate import SpectralSurrogate
    from late_geometry import sample_path
    forbidden_calls = []
    def forbidden(*args, **kwargs):
        forbidden_calls.append('forbidden_model_or_theory_call')
        raise AssertionError('This check may not construct or evaluate a native/proposal cosmology.')
    configurations, identities = [], []
    with patch.object(camb, 'get_background', forbidden), patch.object(camb, 'get_results', forbidden), \
         patch.object(camb, 'get_transfer_functions', forbidden), patch.object(camb, 'set_params', forbidden), \
         patch.object(cobaya.model, 'get_model', forbidden), patch.object(cobaya.model.Model, '__init__', forbidden), \
         patch.object(SpectralSurrogate, 'initialize', forbidden):
        lineage = adapter.calibration_lineage()
        for model in ['lcdm', 'cpl']:
            for numerical in [None, Path(surrogate).resolve()]:
                new = adapter.configuration(model=model, surrogate=numerical)
                old = adapter.original_configuration(model=model, evolution='none', sample='dovekie',
                                                     calibration='official_planck', surrogate=numerical)
                original_sn = old['likelihood'].pop('released_sn')
                changed = deepcopy(new); anchored_sn = changed['likelihood'].pop('released_sn')
                assert adapter.canonical(changed) == adapter.canonical(old), 'A non-SN likelihood/prior/theory setting changed.'
                assert new['params'] == old['params'] and new['params']['epsilon'] == 0.
                assert original_sn['external'] is not anchored_sn['external']
                assert anchored_sn['external'] is adapter.AnchoredReleasedSN
                identity = adapter.identify(new, adapter.DATA, numerical)
                assert identity['anchored_calibration'] == lineage
                assert identity['sample_sha256'] == adapter.digest(adapter.DATA)
                assert identity['sample_semantics']['extra_H0_or_Cepheid_factor'] is False
                assert identity['configuration']['likelihood']['released_sn']['external'] == 'anchored_adapter.AnchoredReleasedSN'
                identities.append(identity)
                configurations.append({'model': model, 'numerical_proposal': numerical is not None,
                    'target_identity': identity['identity'], 'all_non_SN_configuration_identical': True,
                    'sampled_prior_dimensions': [k for k, v in new['params'].items() if isinstance(v, dict) and 'prior' in v],
                    'fixed_epsilon': new['params']['epsilon']})
        assert len({x['identity'] for x in identities}) == 4
        # This calls only the algebraic release class, never Cobaya Model. Test
        # provider redshift ordering/duplication independently of the likelihood.
        sn = object.__new__(adapter.AnchoredReleasedSN)
        sn.data_file = str(adapter.DATA)
        sn.initialize()
        req = sn.get_requirements()['angular_diameter_distance']['z']
        assert np.array_equal(req, np.unique(sn.release.z_hd_noncalibrator))
        assert len(sn.inverse_z) == 1580 and len(req) < 1580
        assert int(sn.release.calibrator.sum()) == 77 and len(sn.release.data) == 1657
        assert abs(lineage['likelihood']['logdetC']-sn.release.logdet_covariance) < 1e-9
        assert abs(lineage['likelihood']['flat_M_precision']-sn.release.flat_m_precision) < 1e-8
        calls, checks = [], []
        covariance = sn.release.C
        inv = linalg.cho_solve((sn.release.chol, True), np.eye(len(covariance)))
        one = np.ones(len(covariance)); u = inv @ one; info = float(one @ u)
        for i in range(8):
            def distance(z):
                return (2800+20*i)*z/((1+z)*(1+(.1+.01*i)*z))
            def provider(z):
                assert np.array_equal(z, req)
                calls.append(z.copy())
                return distance(z)
            sn.provider = SimpleNamespace(get_angular_diameter_distance=provider)
            derived = {}
            observed = sn.release.data.m_b_corr.to_numpy()
            mu = sn.release.data.CEPH_DIST.to_numpy().copy()
            z, zh = sn.release.z_hd_noncalibrator, sn.release.z_hel_noncalibrator
            mu[~sn.release.calibrator] = 5*np.log10((1+z)*(1+zh)*distance(z))+25
            residual = observed-mu
            mean = (u @ residual)/info
            residual -= mean
            chi2 = float(residual @ inv @ residual)
            expected = -.5*(chi2+np.linalg.slogdet(covariance)[1]+np.log(info)+(len(mu)-1)*np.log(2*np.pi))
            actual = sn.logp(_derived=derived)
            assert abs(actual-expected) < 1e-7 and abs(derived['sn_chi2']-chi2) < 1e-7
            checks.append({'case': i, 'loglike_absolute_error': float(abs(actual-expected)),
                           'sn_chi2_absolute_error': float(abs(derived['sn_chi2']-chi2))})
        assert len(calls) == 8
        try:
            sn.logp(epsilon=.01)
        except ValueError:
            pass
        else:
            raise AssertionError('Undeclared luminosity extension accepted')
        assert len(calls) == 8, 'Nonzero epsilon must fail before asking for distances.'
        # Stale evidence tests alter only temporary copies of the audit record.
        work = adapter.ROOT/'.work/unified-cosmology/anchored-validation'
        work.mkdir(parents=True, exist_ok=True)
        raw = json.loads(adapter.CALIBRATION_RECORD.read_text())
        rejected = []
        mutations = {
            'wrong_rows': lambda r: r['selection'].update(selected_rows=1656),
            'extra_Cepheid_factor': lambda r: r['likelihood'].update(extra_SH0ES_or_Cepheid_factor=True),
            'fake_replacement_certification': lambda r: r.update(independent_host_factor_replacement_certified=True),
            'changed_release_data_hash': lambda r: next(iter(r['input_sources'].values())).update(sha256='0'*64),
            'changed_likelihood_source_hash': lambda r: r['source_sha256'].update({next(iter(r['source_sha256'])): '0'*64}),
            'missing_covariance_source': lambda r: r['input_sources'].pop(next(p for p in r['input_sources'] if p.endswith('STAT+SYS.cov'))),
        }
        for name, mutate in mutations.items():
            value = deepcopy(raw); mutate(value)
            path = work/(name+'.json'); path.write_text(json.dumps(value, indent=2)+'\n')
            try:
                adapter.calibration_lineage(path)
            except AssertionError:
                rejected.append(name)
            else:
                raise AssertionError('Stale/altered calibration lineage accepted: '+name)
        for name, kwargs in [('other_sample', {'sample': 'dovekie'}), ('luminosity_drift', {'evolution': 'linear'}),
                             ('other_CMB_calibration', {'calibration': 'paper_literal'})]:
            try:
                adapter.configuration(**kwargs)
            except AssertionError:
                rejected.append(name)
            else:
                raise AssertionError('Undeclared target accepted: '+name)
        try:
            sample_path(adapter.SAMPLE)
        except ValueError:
            rejected.append('original_sample_dispatch_is_incompatible')
        else:
            raise AssertionError('Original dispatch unexpectedly accepted anchored sample')
    assert not forbidden_calls
    source_files = [adapter.HERE/name for name in ['anchored_adapter.py', 'anchored_run.py', 'anchored_validate.py', 'anchored-design.json']]
    source_files += [adapter.HERE.parent/'distance_ladder/calibration_interface.py']
    for p in source_files:
        if p.suffix == '.py':
            ast.parse(p.read_text())
    return {'status': 'passed_configuration_identity_and_provider_checks_no_inference',
            'configurations': configurations, 'provider_checks': checks, 'rejected_mutations_or_unsupported_dispatch': rejected,
            'calibration_audit_sha256': adapter.digest(adapter.CALIBRATION_RECORD),
            'anchored_calibration': lineage,
            'SN_selected_row_indices_sha256': hashlib.sha256(sn.release.original_indices.astype('<i8').tobytes()).hexdigest(),
            'SN_calibrator_mask_sha256': hashlib.sha256(sn.release.calibrator.astype('u1').tobytes()).hexdigest(),
            'SN_covariance_float64_sha256': hashlib.sha256(np.asarray(covariance, dtype='<f8').tobytes()).hexdigest(),
            'source_sha256': {adapter.relative(p): adapter.digest(p) for p in source_files},
            'native_or_background_calls': 0, 'Cobaya_model_constructions': 0, 'proposal_theory_initializations': 0,
            'all_actual_external_assets_verified': True,
            'versions': {p: importlib.metadata.version(p) for p in ['numpy','scipy','pandas','camb','cobaya','candl-like','sacc']},
            'input_sha256': {adapter.relative(surrogate): adapter.digest(surrogate)},
            'inference_ready': False,
            'remaining_integration': 'Explicit anchored sampler/correction/current-target qualifier routing and validation are still required. Existing source-frozen original-sample entrypoints are deliberately not modified or bypassed.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--surrogate', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = validate(args.surrogate.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': result['status'], 'configurations': len(result['configurations']), 'native_calls': 0}))


if __name__ == '__main__':
    main()
