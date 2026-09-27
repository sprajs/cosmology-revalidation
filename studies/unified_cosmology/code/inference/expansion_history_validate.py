"""Synthetic background equations and numerical precision, not cosmology data."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
from ctypes import byref, c_double

import numpy as np
from scipy.integrate import quad

from expansion_history import (DESIGN, ROOT, HERE, no_spectra,
    differentiate_hubble, background_history, evaluate_background, pointwise_summary)
from expansion_history import payload_digest, read_cached_background, verify_cache_manifest
from likelihood import expansion_diagnostics


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def analytic_history(z, omega_m, omega_r, w, wa, h0=70.):
    """Differentiate the exact density-polynomial/CPL expression algebraically."""
    x = 1+np.asarray(z)
    dark = (1-omega_m-omega_r)*x**(3*(1+w+wa))*np.exp(-3*wa*(1-1/x))
    logprime = 3*(1+w+wa)/x-3*wa/x**2
    logsecond = -3*(1+w+wa)/x**2+6*wa/x**3
    f = omega_m*x**3+omega_r*x**4+dark
    fp = 3*omega_m*x**2+4*omega_r*x**3+dark*logprime
    fpp = 6*omega_m*x+12*omega_r*x**2+dark*(logprime**2+logsecond)
    return {'H_km_s_Mpc': h0*np.sqrt(f), 'q': .5*x*fp/f-1,
            'j': 1-x*fp/f+.5*x**2*fpp/f}


def main():
    import camb
    import camb.constants as constants
    sys.path.insert(0, str(HERE.parent/'external_probes'))
    from modern_adapter import modern_info
    design = json.loads(DESIGN.read_text())
    z = np.array(design['redshift_grid']); h = design['finite_difference_step']
    power_rows = []
    for power in [0., .5, 1., 1.5, 2., 2.5]:
        numeric = differentiate_hubble(lambda zz: 70.*(1+np.asarray(zz))**power, z, h)
        q = power-1; j = (power-1)*(2*power-1)
        dq = float(np.max(abs(numeric['q']-q)))
        dj = float(np.max(abs(numeric['j']-j)))
        assert dq < 1e-7 and dj < 2e-5
        power_rows.append({'power': power, 'expected_q': q, 'expected_j': j,
                           'maximum_absolute_q_error': dq, 'maximum_absolute_j_error': dj})
    density_rows = []
    for om, ore, w, wa in [(.3, 9e-5, -1., 0.), (.3, 9e-5, -.8, -.6),
                           (.2, .01, -1.5, .4), (.45, .05, -.5, -.2)]:
        exact = analytic_history(z, om, ore, w, wa)
        numeric = differentiate_hubble(lambda zz: analytic_history(zz, om, ore, w, wa)['H_km_s_Mpc'], z, h)
        errors = {key: float(np.max(abs(numeric[key]-exact[key])/(1+abs(exact[key])))) for key in ['q', 'j']}
        assert errors['q'] < 1e-7 and errors['j'] < 2e-5
        density_rows.append({'synthetic_parameters': [om, ore, w, wa], 'scaled_errors': errors})
    # A synthetic massless-neutrino Lambda+matter+radiation CAMB background has
    # closed density equations, while retaining the actual photon temperature.
    point = {'H0': 70., 'ombh2': .0224, 'omch2': .119, 'ns': .965,
             'tau': .055, 'logA': 3.04, 'w': -1., 'wa': 0.}
    extra = modern_info()['theory']['camb']['extra_args']
    toy_extra = dict(extra, mnu=0., num_massive_neutrinos=0)
    args = {k: point[k] for k in ['H0', 'ombh2', 'omch2', 'ns', 'tau', 'w', 'wa']}
    args['As'] = 1e-10*np.exp(point['logA'])
    with no_spectra():
        background = camb.get_background(camb.set_params(**args, **toy_extra))
        om = background.get_Omega('baryon')+background.get_Omega('cdm')
        ore = background.get_Omega('photon')+background.get_Omega('neutrino')
        exact = analytic_history(z, om, ore, -1., 0.)
        numerical = background_history(background, z, h)
        errors = {key: float(np.max(abs(np.array(numerical['history'][key])-exact[key])/(1+abs(exact[key]))))
                  for key in ['H_km_s_Mpc', 'q', 'j']}
        assert errors['H_km_s_Mpc'] < 1e-12 and errors['q'] < 1e-7 and errors['j'] < 2e-5
        # Independent quadrature of the defining age integral in Gyr.
        conversion = constants.Mpc/1000/constants.Gyr
        age_integral = conversion*quad(lambda a: 1/(a*background.hubble_parameter(1/a-1)),
                                       0., 1., epsabs=1e-11, epsrel=1e-10)[0]
        age_error = abs(age_integral-numerical['scalar']['age_Gyr'])
        strict_native_age = float(background.f_DeltaPhysicalTimeGyr(byref(c_double(0.)), byref(c_double(1.)), byref(c_double(1e-10))))
        strict_age_error = abs(age_integral-strict_native_age)
        assert age_error < 1e-7 and strict_age_error < 1e-7
    # Declared modern physical neutrinos/PPF are retained for these synthetic
    # reference points; no CMB likelihood, chain or observational posterior is read.
    modern_rows = []
    for w, wa, h0 in [(-1., 0., 70.), (-.8, -.6, 67.), (-1.4, .2, 73.), (-.5, -1.5, 65.)]:
        current = dict(point, w=w, wa=wa, H0=h0)
        args = {k: current[k] for k in ['H0', 'ombh2', 'omch2', 'ns', 'tau', 'w', 'wa']}
        args['As'] = 1e-10*np.exp(current['logA'])
        with no_spectra():
            background = camb.get_background(camb.set_params(**args, **extra))
            # Existing native diagnostic definition on its separate forward grid.
            native_grid = np.concatenate([np.arange(5)*.001+c for c in [0., .5, 1.]])
            native = expansion_diagnostics(native_grid, background.hubble_parameter(native_grid))
            native.update(omegam=float(background.Params.omegam), rdrag=float(background.get_derived_params()['rdrag']))
            reference_age = conversion*quad(lambda a: 1/(a*background.hubble_parameter(1/a-1)),
                                            0., 1., epsabs=1e-11, epsrel=1e-10)[0]
        result = evaluate_background(current, extra, design, native)
        assert not result['failed_numerical_gates']
        age_quadrature_error = abs(result['scalar']['age_Gyr']-reference_age)
        assert age_quadrature_error < 1e-5
        modern_rows.append({'synthetic_parameters': {'w': w, 'wa': wa, 'H0': h0},
                            'numerical_checks': result['numerical_checks'],
                            'age_quadrature_absolute_Gyr_error': age_quadrature_error})
    # Independently known weighted pointwise quantiles and acceleration fractions.
    values = np.array([[1., -2.], [3., 2.], [5., 4.], [7., 6.]])
    weights = np.array([.1, .2, .3, .4])
    summary = pointwise_summary(values, weights, [.025, .16, .5, .84, .975])
    manual = [np.interp([.025, .16, .5, .84, .975], [.05, .2, .45, .8], col).tolist() for col in values.T]
    assert np.allclose(summary['quantiles_by_redshift_or_scalar'], manual, rtol=0, atol=1e-14)
    assert np.allclose(summary['mean'], [5., 3.8], rtol=0, atol=1e-14)
    from measurement_summary import weighted_fraction
    repeated = np.tile(np.array([-2., -1., 1., 2.]), 2000)
    groups = np.repeat(np.arange(4), 2000)
    sign = weighted_fraction(repeated < 0., np.full(8000, 1/8000), groups)
    assert abs(sign['fraction']-.5) < 1e-12
    assert all(abs(v-.5) < 1e-12 for v in sign['independent_chain_fractions'].values())
    cache_parent = ROOT/'.work/unified-cosmology/inference/expansion-history-validation'
    cache_parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='cache-test-', dir=cache_parent) as temporary:
        directory = Path(temporary); path = directory/'00000.json'
        payload = {'identity': 'synthetic-only', 'index': 0, 'native_record_sha256': 'synthetic-native',
                   'history': {'q': [-.5, .2]}, 'status': 'synthetic_unit_test'}
        original = dict(payload, payload_sha256=payload_digest(payload))
        path.write_text(json.dumps(original))
        assert read_cached_background(path, 'synthetic-only', 0, 'synthetic-native') == original
        manifest = directory/'background-records.json'
        manifest.write_text(json.dumps({str(path.relative_to(ROOT)): digest(path)}))
        verify_cache_manifest(manifest)
        changed = json.loads(path.read_text()); changed['history']['q'][0] = -.7
        path.write_text(json.dumps(changed))
        rejected = []
        for name, operation in [
            ('partial_payload_numeric_tamper', lambda: read_cached_background(path, 'synthetic-only', 0, 'synthetic-native')),
            ('existing_manifest_numeric_tamper', lambda: verify_cache_manifest(manifest))]:
            try:
                operation()
            except AssertionError:
                rejected.append(name)
            else:
                raise AssertionError('Tampered cache was accepted: '+name)
        path.write_text(json.dumps(original))
        try:
            read_cached_background(path, 'different-parent', 0, 'synthetic-native')
        except AssertionError:
            rejected.append('parent_identity_mismatch')
        else:
            raise AssertionError('Wrong parent cache was accepted.')
    sources = [Path(__file__), HERE/'expansion_history.py', DESIGN, HERE/'likelihood.py',
               HERE/'measurement_summary.py', HERE.parent/'external_probes/modern_adapter.py',
               HERE.parent/'external_probes/adapter.py']
    report = {'status': 'passed_synthetic_background_validation', 'observational_points_used': 0,
              'CMB_spectrum_calls': 0, 'power_law_checks': power_rows, 'analytic_CPL_density_checks': density_rows,
              'massless_Lambda_plus_radiation_CAMB_check': {'scaled_errors': errors,
                  'age_quadrature_absolute_Gyr_error': age_error,
                  'strict_native_Romberg_age_absolute_Gyr_error': strict_age_error,
                  'default_native_age_minus_converged_integral_Gyr': numerical['numerical_checks']['age_CAMB_default_minus_integral_Gyr'],
                  'scope': 'Explicit synthetic mnu=0 test of closed density equations, not the adopted physical target.'},
              'modern_physics_synthetic_checks': modern_rows,
              'known_weighted_quantiles_and_sign_fractions_passed': True,
              'cache_tamper_rejections': rejected,
              'cache_review_history': {
                  'pre_fix_source_sha256': '64935e0dad35940f9c71248ffcd4b7032263f7216da082af18cc4e162990221d',
                  'pre_fix_validation_sha256': '226b161b160ddc24e97c476f89b52fac6d23a6d52d1e87b0a35a065cb1f7b946',
                  'finding': 'Independent review found reusable cache checked identity but could reseal altered numerical values.',
                  'resolution': 'Verify old manifest file hashes and per-point canonical checksums, including partial files, before reuse. No observational cache existed.'},
              'source_sha256': {str(p.relative_to(ROOT)): digest(p) for p in sources},
              'versions': {'camb': camb.__version__, 'numpy': np.__version__},
              'initial_failed_check': {'test': 'CAMB default physical_time versus converged age integral within1e-5Gyr',
                  'observed_absolute_difference_Gyr': 0.0002814665199775135,
                  'cause': 'CAMB1.6.6 results.f90 DeltaPhysicalTimeGyr default Romberg tolerance1e-4 at AccuracyBoost1.',
                  'resolution': 'Report converged integral over unchanged H(a), checked with independent Gauss-Legendre and strict native Romberg; retain default discrepancy.'},
              'scope': 'Equation/numerical tests only; no measured expansion history or cosmic age is reported.'}
    path = ROOT/'studies/unified_cosmology/results/inference/expansion-history-validation.json'
    path.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': report['status'], 'synthetic_cases': len(power_rows)+len(density_rows)+len(modern_rows)+1,
                      'maximum_modern_j_step_scaled_error': max(r['numerical_checks']['j_step_scaled_error'] for r in modern_rows)}, indent=2))


if __name__ == '__main__':
    main()
