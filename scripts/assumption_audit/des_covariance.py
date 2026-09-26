#!/usr/bin/env python3
"""Inspect released DES-Dovekie uncertainty modes without refitting light curves.

Only writes under runs/assumption_audit/des. Covariances are reconstructed from
the released packed precision matrices before any differences or selections.
"""
from pathlib import Path
import hashlib
import json
import platform
import numpy as np
import pandas as pd
import scipy
from scipy.linalg import cho_factor, cho_solve, eigh
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.special import ndtr

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'sources/repos/des-science__DES-SN5YR/4_DISTANCES_COVMAT'
OUT = ROOT / 'runs/assumption_audit/des'
INPUTS = []


def precision(path):
    INPUTS.append(path)
    with np.load(path) as f:
        assert {'nsn', 'cov'} <= set(f.files), f.files
        n = int(f['nsn'][0])
        p = np.zeros((n, n))
        p[np.triu_indices(n)] = f['cov']
    return p + np.triu(p, 1).T


def inverse(p):
    return cho_solve(cho_factor(p, lower=True), np.eye(len(p)))


def projected(p):
    q = p @ np.ones(len(p))
    return p - np.outer(q, q) / q.sum()


gx, gw = np.polynomial.legendre.leggauss(64)


def mu(om, z, zh):
    zs = z[:, None] * (gx + 1) / 2
    dc = z / 2 * np.sum(gw / np.sqrt(om * (1 + zs)**3 + 1-om), axis=1)
    return 5 * np.log10(299792.458 / 70 * (1 + zh) * dc) + 25


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    hdpath, metapath = BASE/'DES-Dovekie_HD.csv', BASE/'DES-Dovekie_Metadata.csv'
    INPUTS.extend([hdpath, metapath])
    d = pd.read_csv(hdpath, sep=r'\s+', comment='#', dtype={'CID': str})
    raw = pd.read_csv(metapath, sep=r'\s+', dtype={'CID': str})
    keys = pd.MultiIndex.from_frame(d[['CID', 'IDSURVEY']])
    m = raw.set_index(['CID', 'IDSURVEY']).loc[keys].reset_index()
    assert len(keys.unique()) == len(d)
    assert np.max(abs(m.MU - d.MU)) < 1e-10
    ps, pt = precision(BASE/'STATONLY.npz'), precision(BASE/'STAT+SYS.npz')
    cs, ct = inverse(ps), inverse(pt)
    p = projected(pt)
    z, zh = d.zHD.to_numpy(), d.zHEL.to_numpy()
    y = d.MU.to_numpy()
    fit = minimize_scalar(lambda om: (y-mu(om,z,zh))@p@(y-mu(om,z,zh)), bounds=(.1,.6), method='bounded', options={'xatol':1e-12})
    derivative = (mu(fit.x+1e-5,z,zh)-mu(fit.x-1e-5,z,zh))/2e-5
    response = (p @ derivative) / (derivative @ p @ derivative)
    whiten = 1 / np.sqrt(np.diag(cs))
    records = []
    summed = np.zeros_like(cs)
    vectors = {}
    for path in sorted((BASE/'SingleSYS_CovMatrix').glob('*.npz')):
        c = inverse(precision(path)) - cs
        c = (c + c.T) / 2
        summed += c
        # Statistical whitening avoids BEAMS-downweighted rows dominating norms.
        wc = whiten[:,None] * c * whiten[None,:]
        vals, vecs = eigh(wc, driver='evr')
        # Packed precision values are float32. Inversion/subtraction leaves a
        # roughly 1e-7 whitened eigenvalue noise floor, of both signs.
        tolerance = max(2e-6, vals[-1]*1e-6)
        retain = vals > tolerance
        modes = vecs[:,retain] * np.sqrt(vals[retain])
        unwhitened = modes / whiten[:,None]
        vectors[path.stem] = unwhitened
        reconstructed = unwhitened @ unwhitened.T
        sigma_response = float(np.sqrt(max(0, response @ c @ response)))
        records.append({'name':path.stem, 'positive_rank':int(retain.sum()),
            'negative_rank_beyond_float32_noise':int(np.sum(vals < -tolerance)),
            'rank_threshold_stat_whitened':float(tolerance),
            'min_eigenvalue_stat_whitened':float(vals[0]),
            'positive_eigenvalues_stat_whitened':vals[retain].tolist(),
            'rank_reconstruction_max_abs_mag2':float(np.max(abs(c-reconstructed))),
            'median_diagonal_sigma_mag':float(np.median(np.sqrt(np.maximum(0,np.diag(c))))),
            'fixed_total_GLS_Omega_m_linear_sigma':sigma_response})
        print(path.stem, int(retain.sum()), sigma_response, flush=True)
    target = ct-cs
    difference = target-summed
    wdifference = whiten[:,None]*difference*whiten[None,:]
    wt = whiten[:,None]*target*whiten[None,:]
    missing_eig = eigh(wdifference, eigvals_only=True)
    # Mode alignment after common-intercept marginalization. This is geometric
    # confounding, not evidence of physically correlated nuisance priors.
    names = ['MWEBV','COLORLAW','W22','BS21','P24a','P24b','P24c','CAL_SALT3','MASSLOC','GAMMAEVOL','INTRSC_COLOR']
    alignment = []
    for i, a in enumerate(names):
        ua = vectors[a]
        ga = ua.T@p@ua
        for b in names[i+1:]:
            ub = vectors[b]
            gb = ub.T@p@ub
            va, ea = eigh(ga)
            vb, eb = eigh(gb)
            ia, ib = va>max(1e-12,va[-1]*1e-9), vb>max(1e-12,vb[-1]*1e-9)
            qa, qb = ua@ea[:,ia]/np.sqrt(va[ia]), ub@eb[:,ib]/np.sqrt(vb[ib])
            cosines = np.linalg.svd(qa.T@p@qb, compute_uv=False)
            alignment.append({'a':a,'b':b,'max_abs_canonical_cosine':float(cosines[0]) if len(cosines) else None})
    # A Gaussian mass-measurement approximation is a diagnostic only: not a
    # posterior from the SED fitter, and not a proposed extra distance correction.
    mass, masserr = m.HOST_LOGMASS.to_numpy(),m.HOST_LOGMASS_ERR.to_numpy()
    valid = np.isfinite(mass)&np.isfinite(masserr)&(masserr>0)&(masserr<10)&(mass>0)
    phigh = np.full(len(m),np.nan)
    phigh[valid] = ndtr((mass[valid]-10)/masserr[valid])
    phard = (mass>10).astype(float)
    delta = .033*(phigh-phard)
    mass_rows = pd.DataFrame({'CID':d.CID,'IDSURVEY':d.IDSURVEY,'zHD':z,'mass':mass,'masserr':masserr,'gaussian_p_high':phigh,'hard_step_expectation_difference_mag':delta})
    mass_rows.to_csv(OUT/'host_mass_boundary_diagnostic.csv',index=False)
    boundary = []
    for name,mask in [('all',valid),('DES',valid&(d.IDSURVEY==10)),('low_z',valid&(d.IDSURVEY!=10))]:
        boundary.append({'sample':name,'n_valid':int(mask.sum()),'n_with_0p1_lt_p_high_lt_0p9':int(np.sum(mask&(phigh>.1)&(phigh<.9))),
            'mean_abs_step_expectation_difference_mag':float(np.mean(abs(delta[mask]))),
            'mean_signed_step_expectation_difference_mag':float(np.mean(delta[mask]))})
    result = {'interpretation':'Released covariance audit. Low ranks are finite uncertainty directions, not zero systematic uncertainty. CAL_SALT3 includes CALSPEC; see release_checks.json for subtraction and total closure. VPEC is not PSD as an increment. Linear sigmas are fixed-total GLS diagnostics, not posterior errors or published sigma_w.',
        'native_flat_lcdm_Omega_m':float(fit.x),
        'stat_offdiagonal_max_abs_mag2':float(np.max(abs(cs-np.diag(np.diag(cs))))),
        'stat_sigma_vs_MUERR_max_abs':float(np.max(abs(np.sqrt(np.diag(cs))-d.MUERR))),
        'single_systematics':records,
        'single_covariance_closure':{'count':len(records),'max_abs_difference_mag2':float(np.max(abs(difference))),
            'relative_frobenius_difference':float(np.linalg.norm(difference)/np.linalg.norm(target)),
            'relative_stat_whitened_frobenius_difference':float(np.linalg.norm(wdifference)/np.linalg.norm(wt)),
            'whitened_difference_eigenvalue_min':float(missing_eig[0]),'whitened_difference_eigenvalue_max':float(missing_eig[-1]),
            'Omega_m_systematic_variance_from_full':float(response@target@response),
            'Omega_m_systematic_variance_from_sum':float(response@summed@response)},
        'mode_alignment':sorted(alignment,key=lambda row:row['max_abs_canonical_cosine'] or 0,reverse=True),
        'host_mass_boundary':boundary,
        'calibration_covariance_MC_fractional_SD_under_9_iid_zero_mean_Gaussian_draws':float(np.sqrt(2/9)),
        'dust_covariance_MC_fractional_SD_under_iid_zero_mean_Gaussian_draws':float(np.sqrt(2/3))}
    (OUT/'covariance_results.json').write_text(json.dumps(result,indent=2)+'\n')
    pd.DataFrame(records).drop(columns='positive_eigenvalues_stat_whitened').to_csv(OUT/'single_systematics.csv',index=False)
    sha = lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    codeout=OUT/'code'
    codeout.mkdir(exist_ok=True)
    (codeout/Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    manifest = {'script_sha256':sha(Path(__file__)), 'inputs':{str(f.relative_to(ROOT)):sha(f) for f in INPUTS},
        'code_snapshot':str((codeout/Path(__file__).name).relative_to(ROOT)),
        'outputs':{str((OUT/f).relative_to(ROOT)):sha(OUT/f) for f in ['covariance_results.json','single_systematics.csv','host_mass_boundary_diagnostic.csv']},
        'environment':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__},
        'scope':'No changes to release data, standardization scripts, or phase2.'}
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['single_systematics','mode_alignment']},indent=2))


if __name__ == '__main__':
    main()
