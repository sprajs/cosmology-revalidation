#!/usr/bin/env python3
"""Observed DES selected-host likelihoods from public eight-band deep fields."""
from pathlib import Path
import sys,datetime,hashlib,json
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT))
from lib.records import fitres
W=ROOT/'.work/host-transport';O=ROOT/'studies/host_ages/results/host_transport';BANDS=['U','G','R','I','Z','J','H','KS']
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    hp=W/'des/DES-SN5YR_DIFFIMG_HEAD.FITS.gz';mp=W/'des/DES-Dovekie_Metadata.csv'
    with fits.open(hp) as f:
        a=f[1].data
        keep=['SNID','IAUC','HOSTGAL_OBJID','HOSTGAL_RA','HOSTGAL_DEC','HOSTGAL_DDLR','HOSTGAL_NMATCH','HOSTGAL_NMATCH2','HOSTGAL_CONFUSION','HOSTGAL_SPECZ','REDSHIFT_HELIO','MWEBV','PEAKMJD']+[f'HOSTGAL_{k}_{b}' for k in ['MAG','MAGERR'] for b in 'griz']
        d=pd.DataFrame({k:np.char.strip(np.array(a[k]).astype(str)) if a[k].dtype.kind in 'US' else np.array(a[k]).astype(float) for k in keep})
    meta=fitres(mp);selected=set(meta.index.astype(str));d['selected_Dovekie']=d.SNID.isin(selected)
    assert set(d.loc[d.selected_Dovekie,'SNID'])==set(meta[meta.IDSURVEY==10].index), 'DES metadata must match original candidate HEAD exactly'
    assert d.SNID.is_unique
    d=d[(d.HOSTGAL_OBJID>0)&(abs(d.HOSTGAL_DEC)<=90)&(d.HOSTGAL_RA>=0)&(d.HOSTGAL_RA<360)].copy()
    matches=[];source=[]
    extra=['ID','RA','DEC','TILENAME','EBV_SFD98','FLAGS','FLAGS_NIR','MASK_FLAGS','MASK_FLAGS_NIR','KNN_CLASS','FOF_SIZE']
    phot=[f'BDF_{v}_{b}' for v in ['FLUX_CALIB','FLUX_ERR_CALIB','FLUX_DERED_CALIB','FLUX_ERR_DERED_CALIB','MAG_DERED_CALIB'] for b in BANDS]
    for field in ['C3','E2','X3']:
        p=W/f'des/Y3_DEEP_FIELDS_PHOTOM-SN-{field}-0000.parquet'
        frame=pq.read_table(p,columns=extra+phot).to_pandas()
        assert frame.ID.is_unique and np.isfinite(frame[['RA','DEC']]).all().all()
        cat=SkyCoord(frame.RA.to_numpy()*u.deg,frame.DEC.to_numpy()*u.deg)
        # All hosts queried; membership relies on individual counterpart, never field proximity.
        coords=SkyCoord(d.HOSTGAL_RA.to_numpy()*u.deg,d.HOSTGAL_DEC.to_numpy()*u.deg)
        ix,sep,_=coords.match_to_catalog_sky(cat);_,second,_=coords.match_to_catalog_sky(cat,nthneighbor=2)
        take=sep.arcsec<=2
        q=pd.concat([d[take].reset_index(drop=True),frame.iloc[ix[take]].reset_index(drop=True).add_prefix('deep_')],axis=1)
        q['match_sep_arcsec']=sep.arcsec[take];q['second_sep_arcsec']=second.arcsec[take];q['field']=field
        matches.append(q);source.append({'field':field,'galaxies':len(frame),'path':str(p.relative_to(ROOT)),'sha256':sha(p)})
        print(field,len(frame),'candidate matches <=2arcsec',len(q),flush=True)
    q=pd.concat(matches,ignore_index=True);assert q.SNID.is_unique
    q['primary_match']=(q.match_sep_arcsec<=1)&(q.second_sep_arcsec>1)
    q['deep_quality']=(q[['deep_FLAGS','deep_FLAGS_NIR','deep_MASK_FLAGS','deep_MASK_FLAGS_NIR']]==0).all(axis=1)&(q.deep_KNN_CLASS==1)
    q['host_dlr_lt4']=(q.HOSTGAL_DDLR>=0)&(q.HOSTGAL_DDLR<4)
    # Released distance z and classification belong to selected objects only.
    q=q.merge(meta.reset_index()[['CID','zHEL','FIELD','PROB_SNNV19','HOST_DDLR','HOST_NMATCH']],left_on='SNID',right_on='CID',how='left',validate='one_to_one')
    q['redshift']=np.where(q.selected_Dovekie,q.zHEL,q.REDSHIFT_HELIO)
    q.to_csv(W/'des-deep-crosswalk.csv',index=False)
    rows=[]
    for _,r in q[q.primary_match].iterrows():
        for b in BANDS:
            f=float(r[f'deep_BDF_FLUX_CALIB_{b}']);e=float(r[f'deep_BDF_FLUX_ERR_CALIB_{b}']);fd=float(r[f'deep_BDF_FLUX_DERED_CALIB_{b}']);ed=float(r[f'deep_BDF_FLUX_ERR_DERED_CALIB_{b}'])
            valid=np.isfinite([f,e,fd,ed]).all() and e>0 and ed>0 and all(abs(v)!=9.999e9 for v in [f,e,fd,ed])
            rows.append(dict(sn_id=r.SNID,host_id=int(r.deep_ID),redshift=r.redshift,field=r.field,selected_Dovekie=bool(r.selected_Dovekie),host_dlr_lt4=bool(r.host_dlr_lt4),host_dlr=r.HOSTGAL_DDLR,host_nmatch=r.HOSTGAL_NMATCH,deep_quality=bool(r.deep_quality),host_ra_deg=r.HOSTGAL_RA,host_dec_deg=r.HOSTGAL_DEC,match_sep_arcsec=r.match_sep_arcsec,aperture='global bulge+disk model',band=b.lower(),filter_family='DECam' if b in ['U','G','R','I','Z'] else 'VISTA_VIRCAM',flux_ab_nanomaggies=f/1000 if valid else None,flux_error_ab_nanomaggies=e/1000 if valid else None,flux_dered_ab_nanomaggies=fd/1000 if valid else None,flux_error_dered_ab_nanomaggies=ed/1000 if valid else None,ebv_sfd98=r.deep_EBV_SFD98,measurement_available=bool(valid),negative_flux=bool(valid and f<0),upper_limit_only=False,covariance_available=False,native_zero_point_ab=30,foreground='calibrated observed and released dereddened variants both retained; EBV_SFD98 provided'))
    photframe=pd.DataFrame(rows);photframe.to_csv(W/'des-deep-photometry.csv',index=False)
    selected_phot=photframe[photframe.selected_Dovekie&photframe.deep_quality&photframe.host_dlr_lt4]
    selected_phot.to_csv(W/'des-deep-selected-photometry.csv',index=False)
    counts=[]
    for radius in [.5,1,2]:
        m=(q.match_sep_arcsec<=radius)&(q.second_sep_arcsec>radius)
        for subset in ['all_candidates','Dovekie']:
            use=m&(q.selected_Dovekie if subset=='Dovekie' else True)
            t=q[use];quality=t.deep_quality&t.host_dlr_lt4
            counts.append(dict(radius_arcsec=radius,subset=subset,unique_matches=len(t),ambiguous=int(((q.match_sep_arcsec<=radius)&(q.second_sep_arcsec<=radius)&(q.selected_Dovekie if subset=='Dovekie' else True)).sum()),quality_and_dlr_lt4=int(quality.sum()),distinct_hosts=int(t.deep_ID.nunique()),distinct_quality_hosts=int(t[quality].deep_ID.nunique())))
    coverage=[]
    for field in ['all','C3','E2','X3']:
        for lo,hi in [(.06,.3),(.3,.6),(.6,1.2)]:
            t=q[q.primary_match&q.selected_Dovekie&q.deep_quality&q.host_dlr_lt4&(q.redshift>=lo)&(q.redshift<hi)&((q.field==field) if field!='all' else True)]
            bs={}
            for b in BANDS:
                f=t[f'deep_BDF_FLUX_CALIB_{b}'].to_numpy();e=t[f'deep_BDF_FLUX_ERR_CALIB_{b}'].to_numpy();ok=np.isfinite(f)&np.isfinite(e)&(e>0)&(abs(f)!=9.999e9)&(abs(e)!=9.999e9)
                bs[b]={'measurements':int(ok.sum()),'negative':int(np.sum(ok&(f<0))),'SNR_lt3':int(np.sum(ok&(f/e<3))),'median_snr':float(np.median(f[ok]/e[ok])) if ok.any() else None}
            coverage.append({'field':field,'z_bin':[lo,hi],'events':len(t),'hosts':int(t.deep_ID.nunique()),'bands':bs})
    discrepancy=[];identity=[]
    for b in BANDS:
        f=q[f'deep_BDF_FLUX_DERED_CALIB_{b}'];m=q[f'deep_BDF_MAG_DERED_CALIB_{b}'];ok=(f>0)&np.isfinite(m)&(m<90)
        diff=m[ok]-(30-2.5*np.log10(f[ok]));identity.append({'band':b,'n':int(ok.sum()),'max_abs_mag_error':float(abs(diff).max())})
    for b in 'griz':
        f=q[f'deep_BDF_FLUX_CALIB_{b.upper()}'];m=q[f'HOSTGAL_MAG_{b}'];ok=q.primary_match&q.selected_Dovekie&q.deep_quality&(f>0)&(m>0)&(m<50)
        dif=30-2.5*np.log10(f[ok])-m[ok];discrepancy.append({'band':b,'n':int(ok.sum()),'deep_minus_SN_host_mag_quantiles':np.quantile(dif,[.025,.5,.975]).tolist() if len(dif) else None,'abs_difference_gt_1mag':int(np.sum(abs(dif)>1))})
    ra1=np.deg2rad(q.HOSTGAL_RA);de1=np.deg2rad(q.HOSTGAL_DEC);ra2=np.deg2rad(q.deep_RA);de2=np.deg2rad(q.deep_DEC)
    sep=2*np.arcsin(np.sqrt(np.sin((de2-de1)/2)**2+np.cos(de1)*np.cos(de2)*np.sin((ra2-ra1)/2)**2))*180/np.pi*3600
    err=float(np.max(abs(sep-q.match_sep_arcsec)));assert err<1e-7
    result=dict(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),method=json.loads((Path(__file__).parent/'deep-design.json').read_text()),source_catalogues=source,matching_counts=counts,selected_coverage=coverage,photometric_identity=identity,independent_sky_separation_max_error_arcsec=err,photometry_release_comparison=discrepancy,missingness='Native +/-9.999e9 flux/error placeholders are missing, not signed measurements; all other finite flux with positive error retained.',selection_limit='Unselected parent rows are heterogeneous candidates, not a known rejected-SNIa sample. Catalogue/footprint completeness is not the SNIa selection function.',aperture_limit='Y3 deep coadds differ from transient-free SN host coadds, and shared aperture shape/deblending induces unavailable covariance. Griz comparison is diagnostic, not a selection cut or contamination correction.',output_sha256={str(p.relative_to(ROOT)):sha(p) for p in [W/'des-deep-crosswalk.csv',W/'des-deep-photometry.csv',W/'des-deep-selected-photometry.csv']},input_sha256={str(p.relative_to(ROOT)):sha(p) for p in [hp,mp]},code_sha256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),Path(__file__).parent/'deep-design.json',ROOT/'lib/records.py']})
    (O/'des-deep-interface.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'counts':counts,'coverage':coverage[:3]},indent=2))
if __name__=='__main__':main()
