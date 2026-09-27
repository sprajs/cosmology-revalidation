#!/usr/bin/env python3
"""Signed forced SPIRE map samples with empirical confusion diagnostics."""
from pathlib import Path
import datetime,hashlib,json,sys
import numpy as np,pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import pyarrow.parquet as pq
ROOT=Path(__file__).resolve().parents[4];W=ROOT/'.work/host-transport';O=ROOT/'studies/host_ages/results/host_transport';H=Path(__file__).parent
BANDS=[250,350,500];FWHM=np.array([18.2,24.9,36.3])
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def sample(h,ra,dec):
    wcs=WCS(h['image'].header);x,y=wcs.world_to_pixel_values(ra,dec);ny,nx=h['image'].data.shape
    inside=np.isfinite(x)&np.isfinite(y)&(x>=-.5)&(x<nx-.5)&(y>=-.5)&(y<ny-.5)
    ix=np.where(inside,np.rint(x),0).astype(int);iy=np.where(inside,np.rint(y),0).astype(int)
    f=np.array(h['image'].data[iy,ix],float);e=np.array(h['error'].data[iy,ix],float);flag=np.array(h['flag'].data[iy,ix],int)
    good=inside&np.isfinite(f)&np.isfinite(e)&(e>0)&(flag==0)
    f=np.where(good,f,np.nan);e=np.where(good,e,np.nan)
    rp,dp=wcs.pixel_to_world_values(ix,iy);dist=SkyCoord(ra*u.deg,dec*u.deg).separation(SkyCoord(rp*u.deg,dp*u.deg)).arcsec
    return f,e,good,dist,x,y,flag,inside

def main():
    sm=W/'des-smp-hosts.csv';sp=pd.read_csv(sm,dtype={'SNID':str});sp=sp[sp.selected_Dovekie].copy();assert len(sp)==1623
    qp=W/'des-deep-crosswalk.csv';q=pd.read_csv(qp,dtype={'SNID':str});q=q[q.primary_match&q.selected_Dovekie&q.deep_quality&q.host_dlr_lt4];optical={r.SNID:int(r.deep_ID) for _,r in q.iterrows()}
    # Optical contaminant census. An optical neighbour is not proof of FIR contamination.
    neighbors={}
    for field in ['C3','X3']:
        p=W/f'des/Y3_DEEP_FIELDS_PHOTOM-SN-{field}-0000.parquet';a=pq.read_table(p,columns=['ID','RA','DEC','FLAGS','MASK_FLAGS','KNN_CLASS']).to_pandas();a=a[(a.FLAGS==0)&(a.MASK_FLAGS==0)&(a.KNN_CLASS==1)]
        cat=SkyCoord(a.RA.to_numpy()*u.deg,a.DEC.to_numpy()*u.deg)
        for _,r in q[q.field==field].iterrows():
            dist=SkyCoord(r.HOSTGAL_RA*u.deg,r.HOSTGAL_DEC*u.deg).separation(cat).arcsec;other=a.ID.to_numpy()!=int(r.deep_ID)
            neighbors[r.SNID]={str(b):{'within_half_FWHM':int(np.sum(other&(dist<=fw/2))),'within_FWHM':int(np.sum(other&(dist<=fw)))} for b,fw in zip(BANDS,FWHM)}
    rows=[];matrices=[];coverage=[];sources=[];offsets_out=[]
    j=np.arange(128);r=np.sqrt(60**2+(j+.5)/128*(180**2-60**2));phi=j*np.pi*(3-np.sqrt(5))
    for field,prefix in [('CDFS','CDFS-SWIRE-NEST'),('XMM-LSS','XMM-LSS-NEST')]:
        candidates=sp[(sp.HOSTGAL_RA>=0)&(sp.HOSTGAL_RA<360)&(abs(sp.HOSTGAL_DEC)<=90)].copy()
        # Declared spherical annulus, same offsets at all wavelengths.
        centers=SkyCoord(candidates.HOSTGAL_RA.to_numpy()*u.deg,candidates.HOSTGAL_DEC.to_numpy()*u.deg)
        offset=centers[:,None].directional_offset_by(phi[None,:]*u.rad,r[None,:]*u.arcsec)
        vals=[];errs=[];good=[];dist=[];blank=[];pixels=[];masks=[];inside=[]
        for b in BANDS:
            p=W/'help'/f'{prefix}_image_{b}_SMAP_v6.0.fits';sources.append({'path':str(p.relative_to(ROOT)),'sha256':sha(p)})
            with fits.open(p) as f:
                assert f['image'].header['BUNIT'].replace(' ','')=='Jy/beam' and f['error'].header['BUNIT'].replace(' ','')=='Jy/beam'
                assert not f[0].header['MATCHFLT']
                z=sample(f,centers.ra.deg,centers.dec.deg);v,e,g,ds,x,y,flag,ins=z
                vals.append(v);errs.append(e);good.append(g);dist.append(ds);pixels.append(np.column_stack([x,y]));masks.append(flag);inside.append(ins)
                bs=sample(f,offset.ra.deg,offset.dec.deg)[0];blank.append(bs)
                coverage.append({'field':field,'band_um':b,'inside':int(ins.sum()),'valid':int(g.sum()),'masked_inside':int(np.sum(ins&(flag!=0))),'invalid_flux_error_inside_unmasked':int(np.sum(ins&(flag==0)&~g)),'shape':list(f['image'].data.shape),'unit':f['image'].header['BUNIT']})
        vals=np.array(vals).T;errs=np.array(errs).T;good=np.array(good).T;dist=np.array(dist).T;blank=np.stack(blank,axis=2);allinside=np.any(np.array(inside),axis=0)
        for k,(_,s) in enumerate(candidates.iterrows()):
            if not allinside[k]:continue
            nvalid=np.isfinite(blank[k]).all(axis=1);samples=blank[k,nvalid];n=len(samples);cov=np.cov(samples,rowvar=False,ddof=1) if n>=64 else np.full((3,3),np.nan)
            bg=samples.mean(0) if n else np.full(3,np.nan);med=np.median(samples,axis=0) if n else np.full(3,np.nan)
            response=np.exp(-4*np.log(2)*(dist[k]/FWHM)**2)
            snid=s.SNID;matrices.append({'sn_id':snid,'field':field,'n_common_offsets':n,'covariance_Jy2':cov.tolist(),'meaning':'Empirical3band distribution of common nearby sky positions, including confusion/instrument noise; not divided by N and not a known exact host likelihood. Neighbouring sky samples spatially correlated.'})
            offsets_out.append({'sn_id':snid,'field':field,'flux':blank[k]})
            for ib,b in enumerate(BANDS):
                ns=neighbors.get(snid,{}).get(str(b),{});x,y=pixels[ib][k]
                rows.append({'sn_id':snid,'redshift':s.redshift,'host_id':int(s.HOSTGAL_OBJID),'host_ra_deg':s.HOSTGAL_RA,'host_dec_deg':s.HOSTGAL_DEC,'host_dlr':s.HOSTGAL_DDLR,'host_nmatch':s.HOSTGAL_NMATCH,'field':field,'band_um':b,'eightband_331':snid in optical,'valid_map_measurement':bool(good[k,ib]),'native_peak_Jy_per_beam':vals[k,ib],'instrument_error_Jy_per_beam':errs[k,ib],'background_mean_Jy_per_beam':bg[ib],'background_median_Jy_per_beam':med[ib],'background_subtracted_peak_Jy':vals[k,ib]-bg[ib],'background_subtracted_peak_median_Jy':vals[k,ib]-med[ib],'sky_distribution_sigma_Jy':np.sqrt(cov[ib,ib]),'sky_q95_above_mean_Jy':np.quantile(samples[:,ib]-bg[ib],.95) if n else np.nan,'common_sky_offsets':n,'gaussian_centering_response':response[ib],'gaussian_centering_corrected_Jy':(vals[k,ib]-bg[ib])/response[ib],'FWHM_arcsec':FWHM[ib],'pixel_center_separation_arcsec':dist[k,ib],'pixel_x':x,'pixel_y':y,'map_flag':int(masks[ib][k]),'optical_neighbors_within_half_FWHM':ns.get('within_half_FWHM'),'optical_neighbors_within_FWHM':ns.get('within_FWHM'),'deblended':False,'censored_upper_limit':False})
    df=pd.DataFrame(rows);assert not df.duplicated(['sn_id','band_um']).any();df.to_csv(W/'des-spire-photometry.csv',index=False)
    # Null -> missing covariance, never a fabricated zero variance.
    def clean(x):
        if isinstance(x,list):return [clean(v) for v in x]
        if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
        if isinstance(x,float) and not np.isfinite(x):return None
        return x
    (W/'des-spire-covariance.json').write_text(json.dumps(clean(matrices),indent=2,allow_nan=False)+'\n')
    np.savez_compressed(W/'des-spire-sky-samples.npz',sn_id=np.array([r['sn_id'] for r in offsets_out]),field=np.array([r['field'] for r in offsets_out]),flux_Jy=np.array([r['flux'] for r in offsets_out]),radius_arcsec=r,position_angle_rad=phi)
    stats=[]
    for scope in ['all_selected','eightband331']:
        for b in BANDS:
            t=df[(df.band_um==b)&df.valid_map_measurement&(df.common_sky_offsets>=64)&((df.eightband_331) if scope=='eightband331' else True)]
            v=t.background_subtracted_peak_Jy;sig=t.sky_distribution_sigma_Jy;nb=t.optical_neighbors_within_half_FWHM.dropna()
            stats.append({'scope':scope,'band_um':b,'n':len(t),'negative_background_subtracted':int((v<0).sum()),'median_flux_mJy':float(v.median()*1000),'median_empirical_sky_sigma_mJy':float(sig.median()*1000),'median_instrument_sigma_mJy':float(t.instrument_error_Jy_per_beam.median()*1000),'above3_empirical_sigma':int((v>3*sig).sum()),'below_minus3_empirical_sigma':int((v<-3*sig).sum()),'median_gaussian_response':float(t.gaussian_centering_response.median()),'optical_neighbor_counts_available':len(nb),'at_least1_neighbor_inside_half_FWHM':int((nb>=1).sum())})
    corr=[]
    for m in matrices:
        c=np.array(m['covariance_Jy2']);v=np.sqrt(np.diag(c))
        if np.isfinite(c).all():corr.append(c/np.outer(v,v))
    record={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'design':json.loads((H/'infrared-design.json').read_text()),'coverage':coverage,'statistics':stats,'empirical_local_correlation_median':np.median(corr,axis=0).tolist(),'objects_with_any_map_position':int(df.sn_id.nunique()),'interpretation':'Signed samples at known optical host coordinates are not deblended host FIR flux. Correlated confusion, missing shorter wavelength constraints, dust temperature/emissivity, geometry and AGN remain. No exact Ramaiya501 reproduction and no empirical stellar age or age-brightness correction claimed. Gaussian beam sensitivity is an approximation, not a calibration.', 'sources':sources,'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [sm,qp,H/'infrared-design.json']},'code_sha256':{str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))},'output_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [W/'des-spire-photometry.csv',W/'des-spire-covariance.json',W/'des-spire-sky-samples.npz']}}
    (O/'infrared-interface.json').write_text(json.dumps(clean(record),indent=2,allow_nan=False)+'\n');print(json.dumps(stats,indent=2))
if __name__=='__main__':main()
