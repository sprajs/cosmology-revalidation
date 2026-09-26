"""Independent algebra and attenuation-derivative review of the saved pilot."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.linalg import cholesky, solve_triangular, lstsq
import extinction
import flux_response as subject

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT/'runs/salt_dust_audit/flux_response'
OUTPUT = ROOT/'docs/salt-dust-audit/flux-response-review.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    arrays = np.load(INPUT/'matrices.npz')
    selected = pd.read_csv(INPUT/'selected_objects.csv', dtype={'CID':str})
    ids = sorted({key.split('__')[0] for key in arrays.files})
    diagnostics = []
    for cid in ids:
        key = lambda suffix: arrays[cid+'__'+suffix]
        j, g = key('jacobian_flux'), key('nuisance_flux')
        for weight in ['measurement_only', 'measurement_plus_model']:
            covariance = key(weight+'_covariance')
            chol = cholesky(covariance, lower=True)
            jw, gw = [solve_triangular(chol, a, lower=True) for a in (j,g)]
            # LAPACK QR with pivoting is independent of the author's manual SVD.
            r, _, rank, _ = lstsq(jw, gw, lapack_driver='gelsy')
            saved = key(weight+'_response')
            res = gw-jw@r
            pcov = key(weight+'_parameter_covariance')
            gram = jw.T@jw
            diagnostics.append(dict(CID=cid, weighting=weight, rank=int(rank),
                max_response_difference=float(np.max(np.abs(r-saved))),
                covariance_inverse_residual=float(np.max(np.abs(gram@pcov-np.eye(4)))),
                gray_amplitude_closure=float(np.max(np.abs(r[:,5:].sum(axis=1)-[1,0,0,0]))),
                gray_residual_norm=float(np.linalg.norm(res[:,5:].sum(axis=1))),
                residual_orthogonality_relative=float(np.linalg.norm(jw.T@res)/max(np.linalg.norm(jw)*np.linalg.norm(gw),1)),
                standardized_response_difference=float(np.max(np.abs(np.array([1,.16087,-3.11780,0])@r-key(weight+'_distance_response'))))))
    model, bands, _, _ = subject.build_model()
    header_path = subject.RELEASE/'0_DATA/DES-SN5YR_DES/DES-SN5YR_DES_HEAD.FITS.gz'
    heads = {str(h['SNID']).strip():h for h in fits.getdata(header_path,1)}
    # Inspect the raw interpolated variance before sncosmo's negative-value floor.
    # A positive final covariance alone would not exclude a triggered repair.
    raw_variances = []
    for row in selected.itertuples():
        cid = str(row.CID)
        times, bnames = arrays[cid+'__mjd'], arrays[cid+'__band']
        for bname in np.unique(bnames):
            phase=(times[bnames==bname]-row.PKMJD)/(1+row.zHEL)
            wave_eff=bands[bname].shifted(1/(1+row.zHEL)).wave_eff
            source=model.source
            a=source._model['LCRV00'](phase,wave_eff)[:,0]
            b=source._model['LCRV11'](phase,wave_eff)[:,0]
            d=source._model['LCRV01'](phase,wave_eff)[:,0]
            raw_variances.extend(a+2*row.x1*d+row.x1**2*b)
    checks = []
    for index in [0,len(selected)//2,len(selected)-1]:
        row = selected.iloc[index]
        cid = str(row.CID)
        model.set(z=row.zHEL,t0=row.PKMJD,x0=row.x0,x1=row.x1,c=row.c,
                  mwebv=float(heads[cid]['MWEBV']),mwrv=3.1,hostebv=0.,hostrv=3.1)
        times, bnames = arrays[cid+'__mjd'], arrays[cid+'__band']
        saved_ratios = arrays[cid+'__nuisance_flux'][:,:5]/arrays[cid+'__flux_model'][:,None]
        for bname in np.unique(bnames):
            band=bands[bname]; mask=bnames==bname
            # Independent dense trapezoid quadrature: derivative / original flux
            # is the photon-weighted mean extinction coefficient times -0.4 ln 10.
            wave=np.linspace(band.minwave(),band.maxwave(),int(np.ceil((band.maxwave()-band.minwave())/.5))+1)
            f=model.flux(times[mask],wave)
            weight=f*wave[None,:]*band(wave)[None,:]
            base=np.trapezoid(weight,wave,axis=1)
            predicted=[]
            for rv in [1.5,2.,3.1,4.]:
                k=extinction.fitzpatrick99(wave/(1+row.zHEL),rv,rv)
                predicted.append(-subject.K*np.trapezoid(weight*k[None,:],wave,axis=1)/base)
            kmw=extinction.fitzpatrick99(wave,3.1,3.1)
            predicted.append(-subject.K*np.trapezoid(weight*kmw[None,:],wave,axis=1)/base)
            predicted=np.array(predicted).T
            original=saved_ratios[mask]
            checks.append(dict(CID=cid,band=str(bname),epochs=int(mask.sum()),
                min_unscaled_integrated_flux=float(base.min()),
                max_absolute_derivative_ratio_difference=float(np.max(np.abs(predicted-original))),
                max_relative_derivative_ratio_difference=float(np.max(np.abs((predicted-original)/predicted)))))
    maxima={name:max(d[name] for d in diagnostics) for name in diagnostics[0] if name not in ['CID','weighting','rank']}
    result={'status':'independent bounded numerical review',
        'reviewed_objects':len(ids),'reviewed_weighting_cases':len(diagnostics),
        'all_parameter_tangents_full_rank':all(d['rank']==4 for d in diagnostics),
        'maxima':maxima,'independent_photon_integral_checks':checks,
        'raw_sncosmo_interpolated_quadratic_variance': {
            'epochs':len(raw_variances),
            'negative_count':int(np.sum(np.asarray(raw_variances)<0)),
            'minimum':float(np.min(raw_variances))},
        'negative_model_photometry_epochs':int(sum(np.sum(arrays[cid+'__flux_model']<=0) for cid in ids)),
        'limits':['Noise is frozen; this is not full parameter-dependent covariance likelihood verification.',
                  'Published coordinates and accepted epochs are conditional inputs, not a new SALT fit.',
                  'Quadrature uses the same physical extinction function but an analytic derivative and independent wavelength integration.',
                  'Three redshift-ranked objects receive independent integral checks; algebra is checked for every saved object.',
                  'A flux ratio near the data is not a pixel or SNANA calibration reproduction.'],
        'sha256':{str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),ROOT/'scripts/salt_dust_audit/flux_response.py',INPUT/'manifest.json',INPUT/'matrices.npz',INPUT/'selected_objects.csv']}}
    if not result['all_parameter_tangents_full_rank'] or maxima['max_response_difference']>1e-7:
        raise AssertionError(result)
    if any(c['max_relative_derivative_ratio_difference']>1e-3 for c in checks):
        raise AssertionError(checks)
    OUTPUT.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['reviewed_objects','maxima','negative_model_photometry_epochs']},indent=2))
    print('Independent photon-integral max relative difference:',max(c['max_relative_derivative_ratio_difference'] for c in checks))


if __name__=='__main__':
    main()
