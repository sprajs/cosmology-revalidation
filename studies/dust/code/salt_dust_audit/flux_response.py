#!/usr/bin/env python3
"""Conditional broadband response; not an official refit or BBC rerun.

Run with phase2/env-official/bin/python. Uses published accepted epochs and
released DES SALT3/passbands; all perturbations are hypothetical injections.
"""
from pathlib import Path
import argparse
import gzip
import hashlib
import io
import json
import platform
import sys

import numpy as np
import pandas as pd
import scipy
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import least_squares
from astropy.io import fits
import sncosmo
import extinction

ROOT = Path(__file__).resolve().parents[2]
RELEASE = ROOT / 'sources/repos/des-science__DES-SN5YR@1.3'
OUT = ROOT / 'runs/salt_dust_audit/flux_response'
K = 0.4 * np.log(10.)
PARAMS = ['dm_amplitude_mag', 'dx1', 'dc', 'dt0_observer_day']
MODES = ['host_E_RV1.5', 'host_E_RV2.0', 'host_E_RV3.1', 'host_E_RV4.0',
         'MW_E_RV3.1', 'zp_g_mag', 'zp_r_mag', 'zp_i_mag', 'zp_z_mag']


class VariableF99(sncosmo.PropagationEffect):
    _param_names = ['ebv', 'rv']
    param_names_latex = ['E(B-V)', 'R_V']
    _minwave = 909.1
    _maxwave = 60000.

    def __init__(self):
        self._parameters = np.array([0., 3.1])

    def propagate(self, wave, flux, phase=None):
        ebv, rv = self._parameters
        return flux * 10. ** (-0.4 * extinction.fitzpatrick99(wave, ebv * rv, rv))


def sha(path):
    return hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()


def build_model():
    directory = RELEASE / '2_LCFIT_MODEL/SALT3.DES5YR'
    mapping = {'m0file': 'template_0', 'm1file': 'template_1',
               'clfile': 'color_correction', 'cdfile': 'color_dispersion',
               'lcrv00file': 'lc_variance_0', 'lcrv11file': 'lc_variance_1',
               'lcrv01file': 'lc_covariance_01'}
    paths = {k: directory / f'salt3_{v}.dat.gz' for k, v in mapping.items()}
    streams = {k: io.StringIO(gzip.decompress(p.read_bytes()).decode()) for k, p in paths.items()}
    source = sncosmo.SALT3Source(**streams)
    model = sncosmo.Model(source=source, effects=[VariableF99(), VariableF99()],
                         effect_names=['mw', 'host'], effect_frames=['obs', 'rest'])
    kcor = ROOT / 'phase2/official/inputs/SNDATA_ROOT/kcor/DES/DES-SN5YR/calib_DES-SN5YR_DES.fits.gz'
    bands = {}
    with fits.open(kcor) as hd:
        data = hd['FilterTrans'].data
        for b in 'griz':
            bands[b] = sncosmo.Bandpass(data['wavelength (A)'], data['DES-' + b],
                                       name='audit-DES-' + b, trim_level=0.)
        zpoff = hd['ZPoff'].data.copy()
    return model, bands, list(paths.values()) + [kcor], zpoff


def analyze(row, points, header, model, bands, mag_offsets, check_nonlinear):
    # Public HEAD values already agree with the refit's resampled/scaled SFD.
    # Do not apply the Schlafly 0.86 factor for a second time.
    ebv = float(header['MWEBV'])
    model.set(z=row.zHEL, t0=row.PKMJD, x0=row.x0, x1=row.x1, c=row.c,
              mwebv=ebv, mwrv=3.1, hostebv=0., hostrv=3.1)
    points = points.sort_values(['MJD', 'BAND'], kind='stable').copy()
    b = np.array([bands[x] for x in points.BAND], dtype=object)
    t = points.MJD.to_numpy(float)
    err = points.FLUXCALERR.to_numpy(float)
    # SNANA's SALT3.INFO MAG_OFFSET and KCOR AB-primary magnitudes affect
    # the flux normalization, despite dropping out of dmB derivatives.
    fluxscale = np.array([10.**(-.4*(.27+mag_offsets[x])) for x in points.BAND])
    assert np.all(np.isfinite(err) & (err > 0))
    if not np.all(model.bandoverlap(b)):
        raise ValueError('Published accepted band exceeds model support')
    base_pars = model.parameters.copy()

    def flux(theta=None, host_e=0., host_rv=3.1, mw_delta=0.):
        theta = np.zeros(4) if theta is None else theta
        model.parameters[:] = base_pars
        model.set(x0=row.x0 * np.exp(-K * theta[0]), x1=row.x1 + theta[1],
                  c=row.c + theta[2], t0=row.PKMJD + theta[3],
                  hostebv=host_e, hostrv=host_rv, mwebv=ebv + mw_delta)
        return fluxscale * model.bandflux(b, t, zp=27.5, zpsys='ab')

    f = flux()
    _, modelcov = model.bandfluxcov(b, t, zp=27.5, zpsys='ab')
    modelcov *= np.outer(fluxscale, fluxscale)
    modelcov = (modelcov + modelcov.T) / 2
    steps = np.array([1e-4, 1e-3, 1e-4, .01])
    def jacobian(factor):
        cols = []
        for j, h in enumerate(steps * factor):
            delta = np.zeros(4); delta[j] = h
            cols.append((flux(delta) - flux(-delta)) / (2 * h))
        return np.column_stack(cols)
    j = jacobian(1.)
    jhalf = jacobian(.5)
    derivative_relative_error = float(np.linalg.norm(j - jhalf) / np.linalg.norm(jhalf))
    # Negative host E here is only a central derivative at zero, not a dust prior.
    h = 1e-4
    dust = [(flux(host_e=h, host_rv=rv) - flux(host_e=-h, host_rv=rv)) / (2*h)
            for rv in [1.5, 2., 3.1, 4.]]
    dust += [(flux(mw_delta=h) - flux(mw_delta=-h)) / (2*h)]
    g = np.column_stack(dust + [-K * f * (points.BAND.to_numpy() == x) for x in 'griz'])
    records, archive = [], {'jacobian_flux': j, 'nuisance_flux': g, 'model_covariance': modelcov,
                             'mjd': t, 'band': points.BAND.to_numpy(dtype='U1'),
                             'flux_model': f, 'flux_observed': points.FLUXCAL.to_numpy(float),
                             'flux_error': err}
    # SALT amplitude coordinate has constant conversion to mB; intercept drops out.
    standardizer = np.array([1., .16087, -3.11780, 0.])
    checks = []
    for weighting in ['measurement_only', 'measurement_plus_model']:
        cov = np.diag(err**2)
        if weighting.endswith('plus_model'):
            cov += modelcov
        chol = cholesky(cov, lower=True)
        jw = solve_triangular(chol, j, lower=True)
        gw = solve_triangular(chol, g, lower=True)
        # SVD avoids squaring condition number to estimate the response.
        u, s, vt = np.linalg.svd(jw, full_matrices=False)
        if s[-1] <= s[0] * 1e-10:
            raise ValueError('Rank-deficient fitted tangent space')
        inverse_j = (vt.T / s) @ u.T
        response = inverse_j @ gw
        residual = gw - jw @ response
        norm = np.sum(gw**2, axis=0)
        residual_norm = np.sum(residual**2, axis=0)
        gram = gw.T @ gw
        resid_gram = residual.T @ residual
        distance_response = standardizer @ response
        pcov = (vt.T / s**2) @ vt
        archive.update({weighting + '_covariance': cov, weighting + '_response': response,
                        weighting + '_parameter_covariance': pcov,
                        weighting + '_nuisance_gram': gram,
                        weighting + '_residual_gram': resid_gram,
                        weighting + '_distance_response': distance_response,
                        weighting + '_singular_values': s})
        # Analytic amplitude response and exact gray-calibration degeneracy.
        assert np.allclose(j[:, 0], -K*f, rtol=1e-7, atol=1e-7)
        assert np.allclose(response[:, 5:].sum(axis=1), [1,0,0,0], atol=1e-7)
        assert np.linalg.norm(residual[:, 5:].sum(axis=1)) < 1e-7 * max(np.linalg.norm(gw), 1)
        assert np.max(np.abs(jw.T @ residual)) < 1e-7 * max(np.linalg.norm(jw)*np.linalg.norm(gw),1)
        for k, mode in enumerate(MODES):
            records.append(dict(CID=row.CID, zHEL=row.zHEL, weighting=weighting, mode=mode,
                                epochs=len(t), bands=''.join(sorted(set(points.BAND))),
                                absorbed_fraction=1-residual_norm[k]/norm[k] if norm[k]>0 else None,
                                residual_snr_per_unit=float(np.sqrt(residual_norm[k])),
                                dstandardized_mag_per_unit=distance_response[k],
                                **dict(zip(PARAMS, response[:, k]))))
        if check_nonlinear and weighting.endswith('plus_model'):
            for e in [.001, .01, .05]:
                target = flux(host_e=e, host_rv=3.1)
                objective = lambda theta: solve_triangular(chol, flux(theta)-target, lower=True)
                fit = least_squares(objective, response[:,2]*e, xtol=1e-11, ftol=1e-11, gtol=1e-9,
                                    diff_step=1e-3, max_nfev=150)
                checks.append(dict(CID=row.CID, injected_host_E=e, success=bool(fit.success),
                                   linear_distance_shift=float(distance_response[2]*e),
                                   nonlinear_distance_shift=float(standardizer@fit.x),
                                   parameter_solution=fit.x.tolist(),
                                   residual_chi2=float(fit.fun@fit.fun)))
    diagnostic = dict(CID=row.CID, zHEL=row.zHEL, epochs=len(t),
                      conditional_mwebv=ebv, derivative_relative_error=derivative_relative_error,
                      min_modelcov_eigenvalue=float(np.linalg.eigvalsh(modelcov)[0]),
                      median_observed_model_flux_ratio=float(np.median(points.FLUXCAL.to_numpy()[f>0]/f[f>0])))
    return records, archive, checks, diagnostic


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--count', type=int, default=64)
    parser.add_argument('--output', type=Path, default=OUT)
    args = parser.parse_args()
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=True)
    model, bands, paths, zp = build_model()
    mag_offsets = {str(r['Filter Name'])[-1]:float(r['Primary Mag']) for r in zp}
    metadata = RELEASE/'4_DISTANCES_COVMAT/DES-SN5YR_HD+MetaData.csv'
    photpath = RELEASE/'0_DATA/DES5YR_SALT3_LCFIT.LCPLOT.gz'
    headpath = RELEASE/'0_DATA/DES-SN5YR_DES/DES-SN5YR_DES_HEAD.FITS.gz'
    meta = pd.read_csv(metadata, dtype={'CID':str})
    points = pd.read_csv(photpath, sep=r'\s+', comment='#', dtype={'SNID':str})
    points = points[points.DATA_MODEL==1].copy()
    meta = meta[(meta.IDSURVEY==10)&meta.CID.isin(points.SNID)].sort_values(['zHEL','CID'])
    indexes = np.unique(np.rint(np.linspace(0,len(meta)-1,min(args.count,len(meta)))).astype(int))
    selected = meta.iloc[indexes].copy()
    selected.to_csv(out/'selected_objects.csv', index=False)
    heads = fits.getdata(headpath,1)
    heads = {str(h['SNID']).strip():h for h in heads}
    records, nonlinear, diagnostics, errors = [], [], [], []
    arrays = {}
    for i, row in enumerate(selected.itertuples()):
        try:
            r, a, c, d = analyze(row, points[points.SNID==row.CID], heads[row.CID], model, bands, mag_offsets,
                                 check_nonlinear=(i%8==0))
            records.extend(r); nonlinear.extend(c); diagnostics.append(d)
            arrays.update({row.CID+'__'+k:v for k,v in a.items()})
        except Exception as exc:
            errors.append({'CID':row.CID, 'type':type(exc).__name__, 'message':str(exc)})
        if i%8==0: print(f'{i+1}/{len(selected)} processed; errors={len(errors)}',flush=True)
    pd.DataFrame(records).to_csv(out/'response.csv',index=False)
    pd.DataFrame(nonlinear).to_csv(out/'nonlinear_checks.csv',index=False)
    pd.DataFrame(diagnostics).to_csv(out/'diagnostics.csv',index=False)
    np.savez_compressed(out/'matrices.npz', **arrays)
    summary = {'selected_count':len(selected),'analyzed_count':len(diagnostics),'errors':errors,
               'parameter_order':PARAMS,'nuisance_order':MODES,
               'sign':'Positive nuisance adds extinction or dims calibrated flux. Response is fitted minus original; adopted-MW-model error has opposite sign.',
               'units':'Dust unit is 1 mag E(B-V); zeropoint unit is 1 mag. Derivatives are local, not finite 1-mag forecasts.',
               'conditioning':'Published accepted epochs, original trained SALT3 model, published coordinates, header MWEBV already scaled, fixed redshift, AB passbands with released primary magnitudes and SALT MAG_OFFSET=.27; no BBC/selection/retraining or intrinsic-population likelihood.',
               'weights':'Measurement-only or measurement plus sncosmo SALT3 model covariance frozen at reference point. Not asserted identical to SNANA.',
               'calibration_offsets':[{str(k):str(v) for k,v in zip(zp.names,r)} for r in zp],
               'covariance_recipe':'For a specified nuisance prior S, C_parameter=R S R^T; C_distance=D S D^T. Shared modes require common columns across SNe. This script does not invent S.',
               'nonlinear_success':all(x['success'] for x in nonlinear),
               'derivative_relative_error_max':max((x['derivative_relative_error'] for x in diagnostics),default=None)}
    if records:
        frame=pd.DataFrame(records)
        summary['median_by_mode'] = frame.groupby(['weighting','mode'])[['absorbed_fraction','dstandardized_mag_per_unit','residual_snr_per_unit']].median().reset_index().to_dict('records')
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    paths += [metadata,photpath,headpath,Path(__file__), RELEASE/'2_LCFIT_MODEL/SALT3.DES5YR/SALT3.INFO',
              ROOT/'phase2/env-official/lib/python3.12/site-packages/sncosmo/models.py']
    manifest={'inputs':[{'path':str(p.relative_to(ROOT)),'sha256':sha(p)} for p in paths],
              'environment':dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,sncosmo=sncosmo.__version__,extinction=extinction.__version__),
              'command':sys.argv,'outputs':[{'path':str(p.relative_to(ROOT)),'sha256':sha(p)} for p in sorted(out.iterdir()) if p.is_file() and p.name!='manifest.json']}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:summary[k] for k in ['selected_count','analyzed_count','errors','nonlinear_success','derivative_relative_error_max']},indent=2))
    if errors: raise SystemExit('Some objects failed; see summary.json; no silent exclusions.')


if __name__=='__main__':
    main()
