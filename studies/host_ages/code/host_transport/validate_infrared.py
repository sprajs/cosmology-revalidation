#!/usr/bin/env python3
"""Independent TAN projection, raw-pixel and covariance closure for SPIRE."""
from pathlib import Path
import datetime,hashlib,json
import numpy as np,pandas as pd
from astropy.io import fits
ROOT=Path(__file__).resolve().parents[4];W=ROOT/'.work/host-transport';O=ROOT/'studies/host_ages/results/host_transport'
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    phot=W/'des-spire-photometry.csv';cp=W/'des-spire-covariance.json';npz=W/'des-spire-sky-samples.npz';d=pd.read_csv(phot,dtype={'sn_id':str},float_precision='round_trip');sources=[];checks=[]
    for field,prefix in [('CDFS','CDFS-SWIRE-NEST'),('XMM-LSS','XMM-LSS-NEST')]:
        for band in [250,350,500]:
            p=W/'help'/f'{prefix}_image_{band}_SMAP_v6.0.fits';sources.append(p);t=d[(d.field==field)&(d.band_um==band)&d.valid_map_measurement]
            with fits.open(p) as f:
                hdr=f['image'].header;ra,dec=np.deg2rad(t[['host_ra_deg','host_dec_deg']].to_numpy()).T;ra0,dec0=np.deg2rad([hdr['CRVAL1'],hdr['CRVAL2']]);cd=np.array([[hdr['CD1_1'],hdr['CD1_2']],[hdr['CD2_1'],hdr['CD2_2']]])
                denom=np.sin(dec)*np.sin(dec0)+np.cos(dec)*np.cos(dec0)*np.cos(ra-ra0)
                xi=np.cos(dec)*np.sin(ra-ra0)/denom;eta=(np.sin(dec)*np.cos(dec0)-np.cos(dec)*np.sin(dec0)*np.cos(ra-ra0))/denom
                pixel=np.linalg.solve(cd,np.rad2deg(np.stack([xi,eta]))).T+np.array([hdr['CRPIX1']-1,hdr['CRPIX2']-1]);err=float(abs(pixel-t[['pixel_x','pixel_y']].to_numpy()).max());assert err<1e-8
                x,y=np.rint(pixel).astype(int).T;native=f['image'].data[y,x];noise=f['error'].data[y,x]
                assert np.array_equal(native,t.native_peak_Jy_per_beam.to_numpy());assert np.array_equal(noise,t.instrument_error_Jy_per_beam.to_numpy());assert (f['flag'].data[y,x]==0).all()
                assert np.max(abs(t.background_subtracted_peak_Jy-(t.native_peak_Jy_per_beam-t.background_mean_Jy_per_beam)))<1e-15
                checks.append({'field':field,'band_um':band,'valid_pixels':len(t),'manual_TAN_max_error_pixels':err,'native_flux_error_and_mask_exact':True})
    z=np.load(npz);m=json.loads(cp.read_text());worst=0;mineig=1.;objects=0
    assert len(m)==len(z['sn_id'])
    for i,r in enumerate(m):
        assert r['sn_id']==str(z['sn_id'][i]) and r['field']==str(z['field'][i]);x=z['flux_Jy'][i];x=x[np.isfinite(x).all(axis=1)];n=len(x);assert n==r['n_common_offsets']
        if n<64:continue
        centered=x-x.sum(axis=0)/n;cov=np.einsum('ij,ik->jk',centered,centered)/(n-1);recorded=np.array(r['covariance_Jy2']);err=float(abs(cov-recorded).max());worst=max(worst,err);assert err<1e-16
        mineig=min(mineig,float(np.linalg.eigvalsh(recorded).min()));objects+=1
        t=d[d.sn_id==r['sn_id']].sort_values('band_um');assert len(t)==3
        assert np.max(abs(t.background_mean_Jy_per_beam.to_numpy()-x.mean(0)))<1e-15
        assert np.max(abs(t.sky_distribution_sigma_Jy.to_numpy()-np.sqrt(np.diag(cov))))<1e-15
    assert mineig>0
    subset=d[d.eightband_331];assert subset.sn_id.nunique()==265 and subset.valid_map_measurement.all()
    r={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'passed':True,'independent_image_checks':checks,'covariances_recomputed':objects,'max_covariance_difference_Jy2':worst,'minimum_covariance_eigenvalue_Jy2':mineig,'eightband_overlap_hosts':int(subset.sn_id.nunique()),'limitations':'Algebraic and provenance validation does not turn spatially correlated confusion samples into an exact Gaussian host likelihood or calibrate the assumed Gaussian beam.','input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [phot,cp,npz,*sources]},'code_sha256':{str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))}}
    (O/'infrared-validation.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
if __name__=='__main__':main()
