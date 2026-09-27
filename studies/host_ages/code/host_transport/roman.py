#!/usr/bin/env python3
"""Preserve native Roman2018 flux likelihood and honest calibration metadata."""
from pathlib import Path
import hashlib,json,datetime
import numpy as np
import pandas as pd
from astropy.io import fits,ascii
ROOT=Path(__file__).resolve().parents[4];W=ROOT/'.work/host-transport';OUT=ROOT/'studies/host_ages/results/host_transport'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    p=W/'roman/snprop.fits.gz'
    with fits.open(p) as f:
        d=pd.DataFrame({name:np.array(f[1].data[name]).astype(str) if f[1].data[name].dtype.kind in 'US' else np.array(f[1].data[name]).astype(float) for name in f[1].columns.names})
    assert len(d)==881 and d.name.is_unique
    text=ascii.read(W/'roman/snprop.dat',readme=str(W/'roman/ReadMe'),format='cds',fill_values=[('---','0'),('--','0'),('-','0')]).to_pandas()
    assert d.name.tolist()==text.Name.tolist()
    for b in 'ugriz':
        valid=d['eflux_local_'+b]>0
        np.testing.assert_allclose(d.loc[valid,'flux_local_'+b],text.loc[valid,'F'+b+'l'],rtol=0,atol=5.1e-7,equal_nan=True)
    rows=[];checks={}
    for _,r in d.iterrows():
        snls=r.survey=='SNLS5'
        for aperture in ['local','global']:
            for b in 'ugriz':
                mag=float(r[f'mag_{aperture}_{b}']);me=float(r[f'emag_{aperture}_{b}'])
                valid_mag=np.isfinite(mag) and np.isfinite(me) and 0<mag<50 and 0<me<10
                if aperture=='local':
                    flux=float(r['flux_local_'+b]);error=float(r['eflux_local_'+b]);available=np.isfinite(flux) and np.isfinite(error) and error>0
                    origin='released aperture flux';likelihood='Gaussian native flux, conditional on diagonal statistical errors'
                else:
                    zp=30 if snls else 22.5
                    flux=10**(.4*(zp-mag)) if valid_mag else np.nan
                    error=np.log(10)/2.5*flux*me if valid_mag else np.nan
                    available=valid_mag;origin='derived from released magnitude';likelihood='Gaussian released magnitude; flux error is first-order approximation'
                rows.append(dict(sn_id=r['name'],survey=r.survey,redshift=r.redshift,sn_ra_deg=r.ra,sn_dec_deg=r.dec,host_id=None,
                    aperture=aperture,aperture_radius_kpc=3. if aperture=='local' else None,
                    band=b,filter_family='CFHT_MegaCam_SNLS' if snls else 'SDSS',
                    native_flux=flux if available else None,native_flux_error=error if available else None,
                    native_flux_unit='SNLS Vega calibrated image counts ZP30' if snls else ('AB maggies' if aperture=='local' else 'AB nanomaggies'),
                    native_magnitude=mag if valid_mag else None,native_magnitude_error=me if valid_mag else None,
                    zero_point=30 if snls else (0 if aperture=='local' else 22.5),magnitude_system='Vega' if snls else 'AB',
                    flux_ab_nanomaggies=flux*(1e9 if aperture=='local' else 1) if available and not snls else None,flux_error_ab_nanomaggies=error*(1e9 if aperture=='local' else 1) if available and not snls else None,
                    measurement_available=bool(available),negative_flux=bool(available and flux<0),upper_limit_only=False,
                    covariance_available=False,flux_origin=origin,likelihood=likelihood,
                    mw_correction_state='No catalogue EBV; paper applies Planck MW extinction during SED fitting. Treat these as observed photometry, with MW nuisance unresolved.',
                    calibration_state='SNLS absolute AB conversion unresolved; native Vega ZP30 verified from flux/magnitude' if snls else 'FITS local flux/magnitude identity gives AB maggies (ZP0), contrary to CDS nanomaggies label; convert by 1e9. Global magnitudes are AB.'))
    q=pd.DataFrame(rows);path=W/'roman-photometry.csv';q.to_csv(path,index=False);d.to_csv(W/'roman-properties.csv',index=False)
    for survey,g in q[q.aperture=='local'].groupby('survey'):
        t=g[g.measurement_available & (g.native_flux>0) & g.native_magnitude.notna()]
        zp=t.native_magnitude+2.5*np.log10(t.native_flux)
        expected_zp=30 if survey=='SNLS5' else 0
        assert np.max(np.abs(zp-expected_zp))<1e-5
        checks[survey]={'events':int(g.sn_id.nunique()),'available_measurements':int(g.measurement_available.sum()),'negative_fluxes':int(g.negative_flux.sum()),'implied_zp_median':float(np.median(zp)), 'implied_zp_range':np.quantile(zp,[0,.5,1]).tolist()}
    result=dict(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),objects=881,photometry_rows=len(q),survey_checks=checks,
        output_path=str(path.relative_to(ROOT)),output_sha256=sha(path),properties_path='.work/host-transport/roman-properties.csv',properties_sha256=sha(W/'roman-properties.csv'),
        native_flux_fits_ascii_agree_atol=5.1e-7,ascii_precision_warning='SDSS-family local maggies are rounded to zero at six-decimal ASCII precision; only FITS retains usable flux precision.',sdss_calibration='All positive local SDSS/CfA/CSP flux-magnitude pairs give ZP0 to <1e-9 mag; native FITS fluxes are AB maggies, contrary to the ReadMe nanomaggies label. Multiply by 1e9 for AB nanomaggies.',upper_limit_semantics='No censoring flag released: finite flux including negative values with positive error is a measurement; missing/zero errors are unavailable, not upper limits.',
        covariance_semantics='No interband or local-global covariance released; local/global aperture data overlap and must not be multiplied as independent likelihoods.',
        mw_semantics='§3.5 paper describes applying Planck MW extinction in fitting, not a per-object EBV column. Explicit foreground modeling or external map required.',
        snls_calibration='Paper §3.1.1 states Vega image ZP30, verified numerically; CDS DeltaZP conversion text is not reconciled to physical AB magnitudes. Preserve native flux and block absolute AB conversion until resolved.',
        code_sha256={str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))},input_sha256={str(p.relative_to(ROOT)):sha(p),'.work/host-transport/roman/ReadMe':sha(W/'roman/ReadMe'),'.work/host-transport/roman/snprop.dat':sha(W/'roman/snprop.dat')})
    (OUT/'roman-interface.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
