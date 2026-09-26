"""Outcome-blind, actual-passband rest-SED response geometry.

Only mean/design/covariance keys are read from archives. No observed flux,
chi-square, residual or residual-score sufficient statistic is accessed.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import platform
import sys

import numpy as np
import pandas as pd
import scipy
from scipy.linalg import solve_triangular, eigh, helmert
from scipy.optimize import brentq
from numpy.polynomial.legendre import legvander, leggauss
import sncosmo
from sncosmo.constants import HC_ERG_AA
from sncosmo.utils import integration_grid

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'runs/research_2026_09_26/astra_design'
OUT = ROOT / 'runs/research_2026_09_26/sed_identification'
PROTOCOL = ROOT / 'docs/research-2026-09-26/sed-identification.md'
SPEC = importlib.util.spec_from_file_location('fr', ROOT / 'scripts/salt_dust_audit/flux_response.py')
fr = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fr)
PAIRS = [(n, 0) for n in range(1, 9)] + [(n, 1) for n in range(9)] + [(n, 2) for n in range(9)]
SMALL = [PAIRS.index(p) for p in [(1, 0), (2, 0), (3, 0), (1, 1)]]
GAUGE = helmert(4, full=False).T  # orthonormal zero-sum griz columns
INPUTS, READ_KEYS = {}, {}
GEOMETRY_COEFFICIENTS = {}


def record(path):
    path = Path(path)
    name = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
    if name not in INPUTS:
        INPUTS[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    return path


def read(path, keys):
    path = record(path)
    banned = ['observed', 'data_flux', 'chi2', '_r', 'discovery_u']
    for key in keys:
        assert not any(key == b or key.endswith('__' + b) for b in banned)
        assert 'observed' not in key and 'data_flux' not in key and 'chi2' not in key
    with np.load(path, allow_pickle=False) as data:
        result = {k: data[k] for k in keys}
    READ_KEYS.setdefault(str(path.relative_to(ROOT)), set()).update(keys)
    return result


def redshift_coordinate(z):
    return np.log((1 + z) / 1.5) / np.log(1.5)


def spectral_basis(wave, phase):
    x = 2 * np.log(np.asarray(wave) / 2000.) / np.log(11000. / 2000.) - 1
    w = legvander(x, 8)
    t = legvander(np.tanh(np.asarray(phase) / 20.), 2)
    return np.column_stack([w[:, n] * t[:, q] for n, q in PAIRS])


def amplitude_metric(indices, drift):
    # Fixed quadrature domain; never inferred from outcome or observer vector.
    x, wx = leggauss(48)
    p, wp = leggauss(40)
    phase = 15 + 30 * p
    xx, pp = np.meshgrid(x, phase, indexing='ij')
    wave = 2000 * np.exp((xx.ravel() + 1) * np.log(5.5) / 2)
    b = spectral_basis(wave, pp.ravel())[:, indices]
    weights = (wx[:, None] * wp[None, :] / 4).ravel()
    m = b.T @ (weights[:, None] * b)
    # Uniform v: E[v]=0, E[v^2]=1/3.
    return np.block([[m, m * 0], [m * 0, m / 3]]) if drift else m


def amplitude_summary(coeff, indices, drift):
    metric = amplitude_metric(indices, drift)
    rms = np.sqrt(np.maximum(np.einsum('ik,ij,jk->k', coeff, metric, coeff), 0))
    xx, pp = np.meshgrid(np.linspace(-1, 1, 181), np.linspace(-15, 45, 121), indexing='ij')
    b = spectral_basis(2000 * np.exp((xx.ravel()+1)*np.log(5.5)/2), pp.ravel())[:, indices]
    if drift:
        half = len(indices)
        vals = np.vstack([b @ (coeff[:half] - coeff[half:]), b @ (coeff[:half] + coeff[half:])])
    else:
        vals = b @ coeff
    return {'rms_mag': rms.tolist(), 'max_abs_mag_on_fixed_grid': np.max(abs(vals), axis=0).tolist()}


def broadband(model, bands, offsets, time, band, z, t0, spacing=5.):
    phase = (time - t0) / (1 + z)
    phase_vander = legvander(np.tanh(phase / 20.), 2)
    base = np.empty(len(time))
    raw = np.empty((len(time), len(PAIRS)))
    negative = 0.
    support_metric = np.zeros((len(PAIRS), len(PAIRS)))
    negative_share, negative_epochs = 0., 0
    ni = np.array([n for n, q in PAIRS]); qi = np.array([q for n, q in PAIRS])
    for name in sorted(set(band)):
        idx = np.flatnonzero(band == name)
        bp = bands[name]
        wave, dw = integration_grid(bp.minwave(), bp.maxwave(), spacing)
        rest = wave / (1 + z)
        assert rest.min() >= 2000 and rest.max() <= 11000
        f = model.flux(time[idx], wave)
        weights = wave * bp(wave) * dw / HC_ERG_AA
        norm = 10**(0.4 * 27.5) / sncosmo.get_magsystem('ab').zpbandflux(bp)
        norm *= 10**(-.4 * (.27 + offsets[name]))
        fw = f * (weights * norm)
        base[idx] = fw.sum(axis=1)
        x = 2 * np.log(rest / 2000.) / np.log(5.5) - 1
        spectral = legvander(x, 8)
        ints = fw @ spectral
        for k, (n, q) in enumerate(PAIRS):
            raw[idx, k] = -fr.K * ints[:, n] * phase_vander[idx, q]
        negative = max(negative, float(np.max(np.sum(abs(np.minimum(fw, 0)), axis=1) / np.maximum(abs(base[idx]), 1e-30))))
        neg = np.sum(abs(np.minimum(fw, 0)), axis=1)/np.maximum(np.sum(abs(fw), axis=1), 1e-30)
        negative_share = max(negative_share, float(neg.max()))
        negative_epochs += int(np.sum(neg > 1e-6))
        # Equal accepted-epoch weights, positive photon-throughput weights.
        moment = spectral.T @ ((weights/weights.sum())[:, None]*spectral)
        phase_columns = phase_vander[idx][:, qi]
        support_metric += (phase_columns.T@phase_columns) * moment[np.ix_(ni, ni)]
    diagnostics = dict(negative_fraction_relative_to_signed_mean=negative,
                       negative_absolute_photon_share_max=negative_share,
                       negative_spectral_epochs_gt1e6=negative_epochs,
                       nonpositive_integrated_flux_epochs=int(np.sum(base<=0)))
    return base, raw, diagnostics, support_metric


def decompose(s):
    u, singular, vt = np.linalg.svd(s, full_matrices=False)
    rank = int(np.sum(singular > singular[0] * 1e-10))
    inverse = (vt[:rank].T / singular[:rank]) @ u[:, :rank].T
    return u[:, :rank], singular, inverse


def constrained_response_match(a, target, rms_limit):
    """Minimize mismatch with Euclidean coefficient norm <= rms_limit.

    a acts on RMS-normalized spectral coefficients. This is a deterministic
    design constraint, not an empirical prior or an observed-residual fit.
    """
    u, singular, vt = np.linalg.svd(a, full_matrices=False)
    score = u.T @ target
    good = singular > singular[0]*1e-10
    unconstrained = vt[good].T @ (score[good]/singular[good])
    if np.linalg.norm(unconstrained) <= rms_limit:
        return unconstrained, 0.
    def coefficients(penalty):
        return vt.T @ (singular*score/(singular**2+penalty))
    def excess(penalty):
        return np.linalg.norm(coefficients(penalty))-rms_limit
    upper = singular[0]**2
    while excess(upper) > 0:
        upper *= 10
    penalty = brentq(excess, 0, upper, xtol=1e-10, rtol=1e-12)
    value = coefficients(penalty)
    assert abs(np.linalg.norm(value)-rms_limit) < 1e-8
    return value, penalty


def analyze_family(cohorts, indices, drift, family):
    metric = amplitude_metric(indices, drift)
    # Physical RMS normalization improves conditioning without changing span.
    root = np.linalg.cholesky(metric)
    invroot = solve_triangular(root.T, np.eye(len(metric)), lower=False)
    prepared, result = {}, {}
    for cohort, d in cohorts.items():
        s = d['S'][:, indices]
        if drift:
            s = np.column_stack([s, s * d['v'][:, None]])
        a = s @ invroot
        qu, singular, inverse = decompose(a)
        o = d['O']
        ou, os, ov = np.linalg.svd(o, full_matrices=False)
        _, rho, angles_vt = np.linalg.svd(qu.T @ ou, full_matrices=False)
        mapping = invroot @ inverse @ o
        error = o - s @ mapping
        retained = eigh(error.T @ error, o.T @ o, eigvals_only=True)
        assert np.allclose(np.sort(retained), np.sort(1-rho**2), atol=1e-8)
        # Each direction's observer griz norm is .02 mag, not unit whitened flux.
        directions = (ov.T / os) @ angles_vt.T
        directions *= .02 / np.linalg.norm(directions, axis=0)
        amplitudes = amplitude_summary(mapping @ directions, indices, drift)
        fixed = d['fixed']
        theta = invroot @ inverse @ fixed
        remainder = fixed - s @ theta
        support_indices = indices + [26+j for j in indices] if drift else indices
        support_metric = d['support_metric'][np.ix_(support_indices, support_indices)]
        support_rms = lambda coeff: np.sqrt(np.maximum(np.einsum('ik,ij,jk->k', coeff, support_metric, coeff), 0)).tolist()
        amplitudes['observed_support_rms_mag'] = support_rms(mapping@directions)
        result[cohort] = {
            'columns': len(metric), 'numerical_rank': qu.shape[1],
            'rms_normalized_singular_values': singular.tolist(),
            'canonical_correlations': rho.tolist(),
            'remaining_observer_information_eigenvalues': retained.tolist(),
            'canonical_griz_dimming_mag': (GAUGE @ directions).T.tolist(),
            'canonical_sed_match': amplitudes,
            'canonical_observer_signal_norm': np.linalg.norm(o @ directions, axis=0).tolist(),
            'canonical_unmatched_signal_norm': np.linalg.norm(error @ directions, axis=0).tolist(),
            'frozen_vector_secondary': {
                'information': float(fixed @ fixed),
                'matched_information_fraction': float(1-remainder@remainder/(fixed@fixed)),
                'remaining_signal_norm': float(np.linalg.norm(remainder)),
                'observed_support_rms_mag': support_rms(theta[:, None]),
                **amplitude_summary(theta[:, None], indices, drift)},
            'fixed_rms_penalty_sensitivities': {}, 'hard_rms_amplitude_frontier': {}}
        for sigma in [.01, .02, .05, .10]:
            regmap = invroot @ np.linalg.solve(a.T @ a + np.eye(a.shape[1])/sigma**2, a.T @ o)
            regerr = o - s @ regmap
            ev = eigh(regerr.T @ regerr, o.T @ o, eigvals_only=True)
            regcoef = regmap @ directions
            result[cohort]['fixed_rms_penalty_sensitivities'][str(sigma)] = {
                'remaining_observer_information_eigenvalues': ev.tolist(),
                'canonical_match_amplitude': amplitude_summary(regcoef, indices, drift)}
            targets = np.column_stack([o@directions, fixed])
            beta_penalty = [constrained_response_match(a, targets[:, j], sigma) for j in range(4)]
            beta = np.column_stack([bp[0] for bp in beta_penalty])
            coeff = invroot @ beta
            mismatch = targets-s@coeff
            info = np.sum(targets**2, axis=0)
            result[cohort]['hard_rms_amplitude_frontier'][str(sigma)] = {
                'target_order': ['canonical1', 'canonical2', 'canonical3', 'frozen_vector_secondary'],
                'remaining_information_fraction': (np.sum(mismatch**2, axis=0)/info).tolist(),
                'remaining_signal_norm': np.linalg.norm(mismatch, axis=0).tolist(),
                'penalty_lagrange_multiplier': [bp[1] for bp in beta_penalty],
                'observed_support_rms_mag': support_rms(coeff),
                **amplitude_summary(coeff, indices, drift)}
            GEOMETRY_COEFFICIENTS[f'{family}_{cohort}_rms{sigma}'] = coeff
        prepared[cohort] = dict(S=s, mapping=mapping, inverse=invroot @ inverse, directions=directions)
    d, v = prepared['discovery'], prepared['validation']
    err = cohorts['validation']['O'] - v['S'] @ d['mapping']
    theta_d = d['inverse'] @ cohorts['discovery']['fixed']
    fixederr = cohorts['validation']['fixed'] - v['S'] @ theta_d
    result['discovery_geometry_mapping_transferred_to_validation'] = {
        'remaining_observer_information_eigenvalues': eigh(err.T @ err, cohorts['validation']['O'].T @ cohorts['validation']['O'], eigvals_only=True).tolist(),
        'frozen_vector_remaining_information_fraction': float(fixederr@fixederr/(cohorts['validation']['fixed']@cohorts['validation']['fixed'])),
        'frozen_vector_remaining_signal_norm': float(np.linalg.norm(fixederr)),
        'qualification': 'Coefficients map response columns to response columns; no observed outcome is fitted.'}
    return result


def finite_checks(model, bands, offsets, ids, design, frozen_griz):
    """Finite SED integration; fixed SALT projection, no nuisance refit."""
    labels = [f'{source}_rms{limit}' for source in ['discovery', 'validation']
              for limit in [.01, .02, .05, .10]]
    coefficients = np.column_stack([GEOMETRY_COEFFICIENTS['broad_drift_'+label][:, 3]
                                   for label in labels])
    exact, linear, maxm = [], [], np.zeros(len(labels))
    distance_rows, nonpositive = [], []
    standardize = np.array([1., .16087, -3.1178, 0.])
    direct_linear_error = 0.
    offset = 0
    for i, cid in enumerate(ids):
        nom = read(BASE/'validation1020/objectives'/f'objective_{cid}.npz',
                   ['MJD', 'band', 'model_flux', 'zHEL', 'MWEBV', 'parameters_x0_x1_c_t0'])
        order = pd.DataFrame({'t': nom['MJD'], 'b': nom['band']}).sort_values(['t', 'b'], kind='stable').index.to_numpy()
        time, band, official = (nom[k][order] for k in ['MJD', 'band', 'model_flux'])
        z, ebv = float(nom['zHEL'][0]), float(nom['MWEBV'][0])
        x0, x1, c, t0 = nom['parameters_x0_x1_c_t0']
        model.set(z=z, t0=t0, x0=x0, x1=x1, c=c, mwebv=ebv, mwrv=3.1, hostebv=0., hostrv=3.1)
        local = coefficients[:26]+redshift_coordinate(z)*coefficients[26:]
        phase_vander = legvander(np.tanh((time-t0)/(1+z)/20), 2)
        delta = np.zeros((len(time), len(labels)))
        delta_linear = np.zeros_like(delta)
        native_base = np.empty(len(time))
        for name in sorted(set(band)):
            idx = np.flatnonzero(band==name); bp = bands[name]
            wave, dw = integration_grid(bp.minwave(), bp.maxwave(), 5.)
            wv = legvander(2*np.log(wave/(1+z)/2000)/np.log(5.5)-1, 8)
            dm = np.zeros((len(idx), len(wave), len(labels)))
            for qdegree in range(3):
                take = [j for j,(n,q) in enumerate(PAIRS) if q==qdegree]
                orders = [PAIRS[j][0] for j in take]
                dspectral = wv[:, orders]@local[take]
                dm += phase_vander[idx, qdegree, None, None]*dspectral[None]
            weights = wave*bp(wave)*dw/HC_ERG_AA
            norm = 10**(.4*27.5)/sncosmo.get_magsystem('ab').zpbandflux(bp)
            norm *= 10**(-.4*(.27+offsets[name]))
            fw = model.flux(time[idx], wave)*(weights*norm)
            base = fw.sum(axis=1)
            native_base[idx] = base
            delta[idx] = (official[idx]/base)[:, None]*np.einsum('tw,twk->tk', fw, np.expm1(-fr.K*dm))
            delta_linear[idx] = (official[idx]/base)[:, None]*np.einsum('tw,twk->tk', fw, -fr.K*dm)
            nonzero = bp(wave)>0
            maxm = np.maximum(maxm, np.max(abs(dm[:, nonzero]), axis=(0,1)))
        old = read(BASE/'validation1020/analysis/objects'/f'{cid}.npz', ['jacobian_flux', 'exact_covariance', 'quoted_error'])
        j = old['jacobian_flux'].copy(); j[:,0] = -fr.K*official
        chol = np.linalg.cholesky(old['exact_covariance'])
        u, singular, vt = np.linalg.svd(solve_triangular(chol,j,lower=True),full_matrices=True)
        q = u[:,4:]
        projected = q.T@solve_triangular(chol,delta,lower=True)
        exact.append(projected)
        ss = design['S'][offset:offset+q.shape[1]]
        linear.append(ss@local)
        direct_linear_error = max(direct_linear_error, float(np.max(abs(q.T@solve_triangular(chol,delta_linear,lower=True)-linear[-1]))))
        griz = np.column_stack([-fr.K*official*(band==b) for b in 'griz'])
        all_delta = np.column_stack([delta_linear, delta, griz@frozen_griz])
        inverse = (vt.T/singular)@u[:,:4].T
        dp = -inverse@solve_triangular(chol,all_delta,lower=True)
        for column, label in enumerate(['linear_'+name for name in labels]+['finite_'+name for name in labels]+['frozen_observer']):
            distance_rows.append(dict(CID=cid,zHEL=z,experiment=label,delta_mB=dp[0,column],delta_x1=dp[1,column],delta_c=dp[2,column],delta_t0=dp[3,column],delta_fixed_reference_preBBC=standardize@dp[:,column]))
        for row in np.flatnonzero((native_base<=0)|(official<=0)):
            nonpositive.append(dict(CID=cid,band=band[row],zHEL=z,phase=(time[row]-t0)/(1+z),native_mean=native_base[row],official_mean=official[row],quoted_flux_error=old['quoted_error'][row],native_abs_mean_over_error=abs(native_base[row])/old['quoted_error'][row],official_abs_mean_over_error=abs(official[row])/old['quoted_error'][row]))
        offset += q.shape[1]
        if i%200==0:
            print('finite check',i+1,'/',len(ids),flush=True)
    exact, linear = np.vstack(exact), np.vstack(linear)
    assert offset==len(design['fixed'])
    target = design['fixed'][:,None]
    report = {}
    supported_rms = np.sqrt(np.einsum('ik,ij,jk->k',coefficients,design['support_metric'],coefficients))
    for k,label in enumerate(labels):
        a,b=exact[:,k],linear[:,k]
        report[label] = {'linear_response_norm':float(np.linalg.norm(b)),
                        'finite_response_norm':float(np.linalg.norm(a)),
                        'finite_minus_linear_norm':float(np.linalg.norm(a-b)),
                        'relative_finite_minus_linear_norm':float(np.linalg.norm(a-b)/np.linalg.norm(b)),
                        'linear_target_remaining_information_fraction':float(np.sum((target[:,0]-b)**2)/np.sum(target**2)),
                        'finite_target_remaining_information_fraction':float(np.sum((target[:,0]-a)**2)/np.sum(target**2)),
                        'validation_throughput_weighted_supported_wavelength_rms_mag':float(supported_rms[k]),
                        'max_abs_sed_change_on_observed_support_mag':float(maxm[k])}
    np.savez_compressed(OUT/'finite-sed-responses.npz', labels=np.array(labels), exact=exact, linear=linear)
    df = pd.DataFrame(distance_rows)
    df.to_csv(OUT/'constructed-distance-responses.csv',index=False)
    pd.DataFrame(nonpositive).to_csv(OUT/'nonpositive-model-epochs.csv',index=False)
    member = df[df.experiment=='frozen_observer'].sort_values(['zHEL','CID'],kind='stable').reset_index(drop=True)
    low,high = set(member.CID.iloc[:255]),set(member.CID.iloc[-255:])
    membership = member[['CID','zHEL']].assign(low_quartile=member.CID.isin(low),high_quartile=member.CID.isin(high))
    prior_membership = ROOT/'runs/research_2026_09_26/shared_distance_response_12/contrast-membership.csv'
    matched_existing = None
    if prior_membership.exists():
        old_member = pd.read_csv(record(prior_membership),dtype={'CID':str}).set_index('CID').sort_index()
        assert np.array_equal(old_member[['low_quartile','high_quartile']],membership.set_index('CID').sort_index()[['low_quartile','high_quartile']])
        matched_existing = True
    membership.to_csv(OUT/'distance-contrast-membership.csv',index=False)
    contrasts = {}
    for name, data in df.groupby('experiment',sort=False):
        contrasts[name] = {'high255_minus_low255_mag':float(data.loc[data.CID.isin(high),'delta_fixed_reference_preBBC'].mean()-data.loc[data.CID.isin(low),'delta_fixed_reference_preBBC'].mean()), 'absolute_parameter_response_quantiles_0_50_95_100': {p: np.quantile(abs(data[p]),[0,.5,.95,1]).tolist() for p in ['delta_mB','delta_x1','delta_c','delta_t0']}}
    assert direct_linear_error < 1e-8
    return {'scope':'Finite emitted SED change with fixed official mean normalization and fixed linear SALT projection. No nonlinear nuisance refit.', 'validation_objects':len(ids), 'experiments':report, 'direct_linear_derivative_check_max_error':direct_linear_error,
            'constructed_distance_response': {'scope':'Optional geometry consequence, not an inferred correction: delta_p=-Jw^+ L^-1 delta_flux, fixed alpha=.16087 beta=3.1178, no BBC/nonlinear fit/retraining/selection.', 'matched_existing_high255_low255_membership':matched_existing,'experiments':contrasts},
            'nonpositive_model_epochs':len(nonpositive)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze', action='store_true')
    parser.add_argument('--compute', action='store_true')
    parser.add_argument('--amplitude-amendment', action='store_true')
    args = parser.parse_args()
    if args.freeze:
        OUT.mkdir(parents=True, exist_ok=False)
        (OUT/'protocol.md').write_bytes(PROTOCOL.read_bytes())
        (OUT/'frozen-protocol.json').write_text(json.dumps({'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'protocol_sha256': hashlib.sha256(PROTOCOL.read_bytes()).hexdigest(), 'pairs': PAIRS, 'small_indices': SMALL}, indent=2)+'\n')
    if not args.compute:
        return
    assert (OUT/'frozen-protocol.json').exists()
    assert not (OUT/'result.json').exists() or args.amplitude_amendment
    record(__file__); record(OUT/'protocol.md'); record(OUT/'frozen-protocol.json')
    record(ROOT/'scripts/salt_dust_audit/flux_response.py')
    record(Path(sncosmo.__file__).parent/'models.py')
    model, bands, paths, zp = fr.build_model()
    for path in paths:
        record(path)
    offsets = {str(r['Filter Name'])[-1]: float(r['Primary Mag']) for r in zp}
    coef = read(BASE/'validation1020/frozen-discovery-coefficients.npz', ['basis_mean', 'gauge_griz'])
    frozen_griz = coef['gauge_griz'] @ coef['basis_mean']
    dc = BASE/'exact43/comparison/matched-matrices.npz'
    with np.load(record(dc)) as data:
        discovery_ids = sorted({k.split('__')[0] for k in data.files})
    validation = pd.read_csv(record(BASE/'validation1020/cohort.csv'), usecols=['CID'], dtype=str)
    cohorts, checks, object_rows = {}, [], []
    for cohort, ids in [('discovery', discovery_ids), ('validation', validation.CID.tolist())]:
        stack = {key: [] for key in ['S', 'S_unmapped', 'O', 'fixed', 'v', 'legacy']}
        cohort_support = np.zeros((52, 52)); epoch_count = 0
        for i, cid in enumerate(ids):
            directory = 'exact43' if cohort == 'discovery' else 'validation1020'
            nom = read(BASE/directory/'objectives'/f'objective_{cid}.npz', ['MJD', 'band', 'model_flux', 'zHEL', 'MWEBV', 'parameters_x0_x1_c_t0'])
            order = pd.DataFrame({'t': nom['MJD'], 'b': nom['band']}).sort_values(['t', 'b'], kind='stable').index.to_numpy()
            time, band, official = (nom[k][order] for k in ['MJD', 'band', 'model_flux'])
            z, ebv = float(nom['zHEL'][0]), float(nom['MWEBV'][0])
            x0, x1, c, t0 = nom['parameters_x0_x1_c_t0']
            model.set(z=z, t0=t0, x0=x0, x1=x1, c=c, mwebv=ebv, mwrv=3.1, hostebv=0., hostrv=3.1)
            if cohort == 'discovery':
                names = ['jacobian_flux', 'exact_covariance', 'flux_model', 'official_mean_exactC_T', 'mjd', 'band']
                old = read(dc, [cid+'__'+k for k in names]); old = {k: old[cid+'__'+k] for k in names}
                native, legacy, oldtime = old['flux_model'], old['official_mean_exactC_T'], old['mjd']
            else:
                names = ['jacobian_flux', 'exact_covariance', 'native_flux_model', 'published_mask_T', 'MJD', 'band']
                old = read(BASE/'validation1020/analysis/objects'/f'{cid}.npz', names)
                native, legacy, oldtime = old['native_flux_model'], old['published_mask_T'], old['MJD']
            assert np.array_equal(time, oldtime) and np.array_equal(band, old['band'])
            base, raw, spectral_checks, support = broadband(model, bands, offsets, time, band, z, t0)
            vz = redshift_coordinate(z)
            cohort_support += np.block([[support, vz*support], [vz*support, vz*vz*support]])
            epoch_count += len(time)
            meanerr = float(np.max(abs(base-native)/np.maximum(abs(native),1e-20)))
            assert meanerr < 1e-10, (cid, meanerr)
            unmapped = raw.copy()
            raw *= (official/base)[:, None]
            jac = old['jacobian_flux'].copy(); jac[:, 0] = -fr.K*official
            chol = np.linalg.cholesky(old['exact_covariance'])
            jw = solve_triangular(chol, jac, lower=True)
            u, sv, vt = np.linalg.svd(jw, full_matrices=True)
            assert np.sum(sv > sv[0]*1e-10) == 4
            q = u[:, 4:]
            proj = lambda value: q.T @ solve_triangular(chol, value, lower=True)
            griz = np.column_stack([-fr.K*official*(band==b) for b in 'griz'])
            obs = proj(griz@GAUGE)
            original_gauge = coef['gauge_griz']
            observer_match = float(np.max(abs(proj(griz@original_gauge)-legacy[:, :3])))
            assert observer_match < 1e-8
            gray = float(np.linalg.norm(proj(-fr.K*official)))
            assert gray < 1e-8
            quaderr = None
            if i % 50 == 0:
                basefine, fine, _, _ = broadband(model, bands, offsets, time, band, z, t0, spacing=2.5)
                fine *= (official/basefine)[:, None]
                quaderr = float(np.linalg.norm(proj(fine-raw))/max(np.linalg.norm(proj(raw)),1e-20))
                assert quaderr < .002, (cid, quaderr)
            stack['S'].append(proj(raw)); stack['S_unmapped'].append(proj(unmapped)); stack['O'].append(obs)
            stack['fixed'].append(proj(griz@frozen_griz)); stack['v'].append(np.full(q.shape[1], redshift_coordinate(z)))
            stack['legacy'].append(legacy[:, 3:])
            checks.append(dict(cohort=cohort, CID=cid, native_mean_relative_error=meanerr, archived_observer_max_error=observer_match, gray_norm=gray, finer_grid_projected_relative_error=quaderr, **spectral_checks))
            object_rows.append(dict(cohort=cohort, CID=cid, zHEL=z, x1=x1, c=c, epochs=len(time), projected_rows=q.shape[1], phase_min=float(((time-t0)/(1+z)).min()), phase_max=float(((time-t0)/(1+z)).max())))
            if i % 100 == 0:
                print(cohort, i+1, '/', len(ids), flush=True)
        cohorts[cohort] = {key: np.concatenate(value, axis=0) for key, value in stack.items()}
        cohorts[cohort]['support_metric'] = cohort_support/epoch_count
    np.savez_compressed(OUT/'design-matrices.npz', **{cohort+'_'+key: value for cohort,d in cohorts.items() for key,value in d.items()})
    pd.DataFrame(checks).to_csv(OUT/'numerical-checks.csv',index=False)
    pd.DataFrame(object_rows).to_csv(OUT/'object-ledger.csv',index=False)
    result = {'scope': 'Local no-residual broadband response geometry; not mechanism identification or a p-value.', 'frozen_griz_mag_secondary': frozen_griz.tolist(), 'families': {}, 'legacy_rest_comparator': {}, 'environment': dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, sncosmo=sncosmo.__version__)}
    for label, indices, drift in [('small_shared', SMALL, False), ('small_drift', SMALL, True), ('broad_shared', list(range(26)), False), ('broad_drift', list(range(26)), True)]:
        result['families'][label] = analyze_family(cohorts, indices, drift, label)
        print('analyzed', label, flush=True)
    np.savez_compressed(OUT/'amplitude-constrained-coefficients.npz', **GEOMETRY_COEFFICIENTS)
    result['native_spectral_normalization_sensitivity'] = {}
    for name,d in cohorts.items():
        o = d['O']; oq, _, _ = decompose(o)
        ss = np.column_stack([d['S_unmapped'], d['S_unmapped']*d['v'][:, None]])
        su, _, _ = decompose(ss)
        rho = np.linalg.svd(su.T@oq, compute_uv=False)
        result['native_spectral_normalization_sensitivity'][name] = {'broad_drift_unmapped_canonical_correlations': rho.tolist(), 'relative_change_in_projected_broad_matrix': float(np.linalg.norm(d['S_unmapped']-d['S'])/np.linalg.norm(d['S']))}
    for name,d in cohorts.items():
        q, singular, inv = decompose(d['legacy'])
        oq, _, _ = decompose(d['O'])
        rho = np.linalg.svd(q.T@oq, compute_uv=False)
        fixed = d['fixed']; rest = fixed-q@(q.T@fixed)
        result['legacy_rest_comparator'][name] = {'canonical_correlations': rho.tolist(), 'remaining_observer_information_eigenvalues': sorted((1-rho**2).tolist()), 'frozen_vector_matched_information_fraction': float(1-rest@rest/(fixed@fixed))}
    if args.amplitude_amendment:
        record(OUT/'amplitude-reporting-amendment.md')
        result['finite_sed_check'] = finite_checks(model,bands,offsets,validation.CID.tolist(),cohorts['validation'],frozen_griz)
    (OUT/'result.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    outputs = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'}
    (OUT/'manifest.json').write_text(json.dumps({'inputs_sha256': INPUTS, 'output_sha256': outputs, 'archive_keys_read': {k: sorted(v) for k,v in READ_KEYS.items()}, 'argv': sys.argv}, indent=2)+'\n')


if __name__ == '__main__':
    main()
