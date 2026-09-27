"""Identity, support, signed-flux, unit and covariance checks for DESI hosts."""
import datetime, json
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from acquire import ROOT, WORK, RESULT, CODE, sha
from bands import native_vector, measure


def main():
    checked=0
    result_names=['acquisition.json','spectra.json','bands.json','sources.json','independent-review.json']
    for filename in result_names:
        result=json.loads((RESULT/filename).read_text())
        for key,entries in result.items():
            if key.endswith('_sha256') and isinstance(entries,dict):
                for p,h in entries.items():
                    assert sha(ROOT/p)==h,p
                    checked+=1
        for entry in result.get('sources',[]):
            assert sha(ROOT/entry['path'])==entry['sha256']
            checked+=1
    clean=pd.read_csv(WORK/'clean-matches.csv',dtype={'targetid':str})
    hosts=pd.read_csv(WORK/'input-hosts.csv').set_index('oz_OzDES_ID')
    statuses=pd.read_csv(WORK/'host-match-status.csv')
    assert len(statuses)==1088 and statuses.oz_OzDES_ID.nunique()==1088
    assert len(clean)==clean.targetid.nunique()==clean.oz_OzDES_ID.nunique()==55
    assert clean.redshift_delta.abs().max()<=.003
    a=SkyCoord(clean.mean_fiber_ra.to_numpy()*u.deg,clean.mean_fiber_dec.to_numpy()*u.deg)
    h=hosts.loc[clean.oz_OzDES_ID]
    b=SkyCoord(h.oz_RA.to_numpy()*u.deg,h.oz_DEC.to_numpy()*u.deg)
    sep=a.separation(b).arcsec
    assert max(sep)<=1 and np.allclose(sep,clean.separation_arcsec,atol=1e-10,rtol=0)
    records=json.loads((WORK/'band-records.json').read_text())
    min_eigen=1.;gradient_error=0.;scale_error=0.;max_snr=0.
    rng=np.random.default_rng(927611)
    for record in records:
        cov=np.array(record['band_covariance_formal_diagonal_ivar']);assert np.allclose(cov,cov.T,atol=1e-13)
        scale=np.sqrt(np.diag(cov));corr=cov/np.outer(scale,scale)
        min_eigen=min(min_eigen,float(np.linalg.eigvalsh(corr).min()))
        assert np.linalg.eigvalsh(corr).min()>-1e-10
        assert cov[1,2]!=0 and cov[1,3]!=0, 'Shared Dn/Hdelta photons lost from covariance.'
        target=record['targetid'];meta=json.loads((WORK/f'spectra/{target}.json').read_text())
        assert meta['native_coordinate_gate'] and meta['native_objtype']=='TGT' and meta['native_coadd_fiberstatus']==0
        assert set(meta['flux_units'].values())=={'10**-17 erg/(s cm2 Angstrom)'}
        arrays=np.load(ROOT/meta['array_file']);v=native_vector(arrays);z=meta['desi_z']
        actual=measure(*v,z)
        max_snr=max(max_snr,max(abs(np.array(actual['signed_pixel_snr_percentiles']))))
        scaled=measure(v[0],v[1]*3,v[2]/9,v[3],v[4],z)
        for key in ['Dn4000','Hdelta_native_A','Dn4000_formal_sigma','Hdelta_formal_sigma_A']:
            if actual[key] is not None:
                scale_error=max(scale_error,abs(actual[key]-scaled[key]));assert np.isclose(actual[key],scaled[key],atol=1e-11,rtol=1e-10)
        direction=rng.normal(size=len(v[0]));epsilon=1e-5
        plus=measure(v[0],v[1]+epsilon*direction,v[2],v[3],v[4],z)
        minus=measure(v[0],v[1]-epsilon*direction,v[2],v[3],v[4],z)
        for k,key in enumerate(['Dn4000','Hdelta_native_A']):
            if actual[key] is None:continue
            finite=(plus[key]-minus[key])/(2*epsilon)
            exact=float(actual['gradient_arrays'][k]@direction)
            err=abs(finite-exact)/max(1,abs(exact));gradient_error=max(gradient_error,err)
            assert err<1e-7,(target,key,finite,exact)
        inverse=measure(v[0],-v[1],v[2],v[3],v[4],z)
        assert np.allclose(inverse['band_flux'],-np.array(actual['band_flux']),atol=1e-12)
    wave=np.arange(3600.,9824.1,.8);ivar=np.full(len(wave),1e6);mask=np.zeros(len(wave),dtype=int);arm=np.full(len(wave),'T');z=.4
    flat_nu=measure(wave,4000**2/wave**2,ivar,mask,arm,z)
    flat_lambda=measure(wave,np.ones(len(wave)),ivar,mask,arm,z)
    assert abs(flat_nu['Dn4000']-1)<1e-12
    assert abs(flat_lambda['Hdelta_native_A'])<1e-12
    report={'schema':'calibrated-host-validation-v1','status':'passed','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'hashes_checked':checked,'unique_matched_hosts':len(clean),'all_target_and_native_fiber_quality_gates':True,'all_parent_hosts_accounted':len(statuses),'maximum_coordinate_recalculation_error_arcsec':float(max(abs(sep-clean.separation_arcsec))),'maximum_index_directional_gradient_error':gradient_error,'maximum_global_flux_rescaling_invariance_error':scale_error,'minimum_joint_band_correlation_eigenvalue':min_eigen,'maximum_absolute_pixel_snr_in_index_regions':max_snr,'flat_Fnu_Dn4000':flat_nu['Dn4000'],'flat_Flambda_Hdelta_A':flat_lambda['Hdelta_native_A'],'signed_flux_linearity':True,'cross_index_photon_covariance_retained':True,'error_scope':'Validation of bookkeeping and the stated diagonal-IVAR operator, not proof of complete spectral calibration/noise model.','source_sha256':{str(p.relative_to(ROOT)):sha(p)for p in sorted(CODE.glob('*'))if p.is_file()},'result_sha256':{str((RESULT/p).relative_to(ROOT)):sha(RESULT/p)for p in result_names},'note_sha256':sha(ROOT/'studies/unified_cosmology/notes/calibrated-hosts.md')}
    (RESULT/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items()if not k.endswith('sha256')},indent=2))


if __name__=='__main__':main()
