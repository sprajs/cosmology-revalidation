#!/usr/bin/env python3
"""Measure conditional deblending information using actual positions and beams.

No image flux, dust luminosity or age correction is inferred by this diagnostic.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime,timezone
import hashlib,json,time
from pathlib import Path
import numpy as np,pandas as pd
import pyarrow.parquet as pq
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.ndimage import map_coordinates
from scipy.optimize import least_squares
from scipy.spatial import cKDTree
from scipy.linalg import svd

ROOT=Path(__file__).resolve().parents[4];WORK=ROOT/'.work/infrared-resolution'
INPUT=ROOT/'.work/host-transport';OUT=ROOT/'studies/host_ages/results/infrared_resolution'
HERE=Path(__file__).parent
BANDS=[250,350,500];FWHM={250:18.2,350:24.9,500:36.3}
PSF={};MAPS={}

def sha(p):
    with p.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

def xyz(ra,dec):
    r,d=np.deg2rad(ra),np.deg2rad(dec)
    return np.column_stack([np.cos(d)*np.cos(r),np.cos(d)*np.sin(r),np.sin(d)])

def beams():
    result={};records=[]
    for band,name in zip(BANDS,['PSW','PMW','PLW']):
        path=WORK/f'0x5000241aL_{name}_bgmod10_1arcsec.fits.gz'
        with fits.open(path) as hdus:
            data=hdus['image'].data;head=hdus['image'].header
            assert np.isclose(abs(head['CDELT1'])*3600,1.) and np.isclose(head['CDELT2']*3600,1.)
            # A translation does not change pair overlap. Centre the radial
            # approximation on an independently fitted high-S/N beam core.
            py,px=np.unravel_index(np.nanargmax(data),data.shape);fw=FWHM[band]
            yy,xx=np.mgrid[py-int(fw):py+int(fw)+1,px-int(fw):px+int(fw)+1]
            xx=xx-px;yy=yy-py;values=data[yy+py,xx+px];use=xx**2+yy**2<(.8*fw)**2
            def residual(p):
                amp,x0,y0,sx,sy,rho,bg=p
                dx=(xx[use]-x0)/sx;dy=(yy[use]-y0)/sy
                return amp*np.exp(-.5*(dx*dx+dy*dy-2*rho*dx*dy)/(1-rho*rho))+bg-values[use]
            initial=[float(data[py,px]),0,0,fw/2.355,fw/2.355,0,0]
            fit=least_squares(residual,initial,bounds=([.5,-3,-3,fw/5,fw/5,-.7,-.1],[2,3,3,fw,fw,.7,.1]),xtol=1e-12,ftol=1e-12,gtol=1e-12)
            assert fit.success
            cx,cy=px+fit.x[1],py+fit.x[2]
            radius=np.arange(0,4*fw+.25,.25);theta=np.arange(720)*2*np.pi/720
            samples=map_coordinates(data,[cy+radius[:,None]*np.sin(theta),cx+radius[:,None]*np.cos(theta)],order=1,mode='constant',cval=0)
            profile=samples.mean(axis=1);profile/=profile[0]
            half=np.where(profile<.5)[0][0];emp_fwhm=2*np.interp(.5,profile[half-1:half+1][::-1],radius[half-1:half+1][::-1])
            r=int(np.ceil(4*fw+8));yc,xc=np.mgrid[-r:r+1,-r:r+1]
            template=map_coordinates(data,[cy+yc,cx+xc],order=1,mode='constant',cval=0)
            template/=float(map_coordinates(data,[[cy],[cx]],order=1)[0])
            template[np.hypot(xc,yc)>4*fw]=0
            result[band]={'radius':radius,'profile':profile,'template':template,'center':r}
            records.append(dict(band_um=band,empirical_radial_FWHM_arcsec=float(emp_fwhm),nominal_FWHM_arcsec=fw,
                header_reference_xy=[float(head['CRPIX1']-1),float(head['CRPIX2']-1)],fitted_center_xy=[float(cx),float(cy)],
                position_angle_deg=float(hdus[0].header['POSANGLE']),core_fit_RMS=float(np.sqrt(np.mean(residual(fit.x)**2))),source_sha256=sha(path)))
    return result,records

def project_out(matrix,sky):
    return matrix-sky@(sky.T@matrix)

def inflation(host,nuisance,tolerance=1e-10):
    if nuisance.shape[1]==0:return 1.,0,1.
    U,s,_=svd(nuisance,full_matrices=False,check_finite=False)
    rank=int(np.sum(s>tolerance*s[0])) if len(s) else 0
    resid=host-U[:,:rank]@(U[:,:rank].T@host)
    relative=float(np.linalg.norm(resid)/np.linalg.norm(host))
    return (1/relative if relative>1e-10 else np.nan),rank,relative

def initialize():
    global PSF,MAPS
    PSF,_=beams()
    for field,prefix in [('C3','CDFS-SWIRE-NEST'),('X3','XMM-LSS-NEST')]:
        for band in BANDS:
            path=INPUT/'help'/f'{prefix}_image_{band}_SMAP_v6.0.fits'
            hdus=fits.open(path,memmap=True);h=hdus['image'].header;wcs=WCS(h)
            scale=np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix)))*3600
            assert np.allclose(wcs.pixel_scale_matrix,np.diag([-scale,scale])/3600,rtol=1e-8,atol=1e-12)
            MAPS[field,band]=(hdus,wcs,scale)

def worker(job):
    host,neighbours=job;rows=[];checks=[]
    coord=SkyCoord(host['ra']*u.deg,host['dec']*u.deg)
    other=SkyCoord(neighbours['ra']*u.deg,neighbours['dec']*u.deg)
    distance=coord.separation(other).arcsec;order=np.argsort(distance)
    distance=distance[order];other=other[order]
    for band in BANDS:
        hdus,wcs,scale=MAPS[host['field'],band];fw=FWHM[band];psf=PSF[band]
        take=distance<=2*fw;d=distance[take];near=other[take]
        hx,hy=wcs.world_to_pixel(coord);sx,sy=wcs.world_to_pixel(near)
        r=int(np.ceil(4*fw/scale));iy,ix=np.mgrid[int(np.rint(hy))-r:int(np.rint(hy))+r+1,int(np.rint(hx))-r:int(np.rint(hx))+r+1]
        shape=hdus['image'].data.shape;inside=(ix>=0)&(ix<shape[1])&(iy>=0)&(iy<shape[0]);ix=ix[inside];iy=iy[inside]
        use=(hdus['flag'].data[iy,ix]==0)&np.isfinite(hdus['image'].data[iy,ix])&np.isfinite(hdus['error'].data[iy,ix])&(hdus['error'].data[iy,ix]>0)
        ix,iy=ix[use],iy[use];x=(ix-hx)*scale;y=(iy-hy)*scale
        circle=x*x+y*y<=(4*fw)**2;x,y=x[circle],y[circle]
        sky=np.linalg.qr(np.column_stack([np.ones(len(x)),x/fw,y/fw]))[0]
        xsrc=np.r_[0,(sx-hx)*scale];ysrc=np.r_[0,(sy-hy)*scale]
        delta=np.hypot(x[:,None]-xsrc,y[:,None]-ysrc)
        a=np.interp(delta,psf['radius'],psf['profile'],right=0)
        p=project_out(a,sky);h=p[:,0];ngal=len(d)
        all_vif,rank,rel=inflation(h,p[:,1:])
        pair_vif,pairrank,_=inflation(h,p[:,1:2])
        correlations=[];orient_vif=[]
        for angle in np.arange(4)*np.pi/4:
            dx=x[:,None]-xsrc[:2];dy=y[:,None]-ysrc[:2]
            rx=dx*np.cos(angle)-dy*np.sin(angle);ry=dx*np.sin(angle)+dy*np.cos(angle)
            native=map_coordinates(psf['template'],[psf['center']+ry,psf['center']+rx],order=1,mode='constant',cval=0)
            native=project_out(native,sky);v,_,_=inflation(native[:,0],native[:,1:]);orient_vif.append(v)
        if ngal:
            rho=float(h@p[:,1]/np.linalg.norm(h)/np.linalg.norm(p[:,1]));direct=1/np.sqrt(max(1-rho*rho,1e-300))
            checks.append(abs(direct/pair_vif-1))
            gaussian=1/np.sqrt(1-np.exp(-4*np.log(2)*(d[0]/fw)**2))
            resolution=float(d[0]*np.sqrt(4*np.log(2)/(-np.log(.75))))
        else:rho=0.;gaussian=1.;resolution=None
        # Independent full-rank normal-matrix check where conditioning permits.
        gram=p.T@p;condition=np.linalg.cond(gram)
        normal_difference=None
        if np.isfinite(condition) and condition<1e7:
            direct=np.sqrt(np.linalg.inv(gram)[0,0]*(h@h));normal_difference=float(abs(direct/all_vif-1));assert normal_difference<2e-6
        rows.append(dict(sn_id=host['sn_id'],host_id=host['host_id'],field=host['field'],band_um=band,pixels=len(x),neighbours_within_2FWHM=ngal,
            nearest_separation_arcsec=float(d[0]) if ngal else None,nearest_pair_radial_noise_inflation=pair_vif,all_candidates_noise_inflation=all_vif,
            independent_neighbour_rank=rank,host_residual_template_norm_fraction=rel,nearest_pair_native_orientation_min=float(np.min(orient_vif)),nearest_pair_native_orientation_max=float(np.max(orient_vif)),
            nearest_pair_template_correlation=rho,Gaussian_continuous_pair_inflation=float(gaussian),Gaussian_FWHM_for_pair_inflation2_arcsec=resolution,
            normal_matrix_condition=float(condition),independent_normal_matrix_relative_difference=normal_difference))
    return rows,checks

def main():
    start=time.perf_counter();_,beam_records=beams()
    qpath=INPUT/'des-deep-crosswalk.csv';q=pd.read_csv(qpath,dtype={'SNID':str,'deep_ID':'Int64'});q=q[q.primary_match&q.selected_Dovekie&q.deep_quality&q.host_dlr_lt4]
    ipath=INPUT/'des-spire-photometry.csv';ir=pd.read_csv(ipath,dtype={'sn_id':str});ir=ir[ir.eightband_331&ir.valid_map_measurement&(ir.common_sky_offsets>=64)]
    ids=ir.groupby('sn_id').band_um.nunique();ids=set(ids[ids==3].index);q=q[q.SNID.isin(ids)];assert len(q)==265 and q.deep_ID.nunique()==265
    jobs=[];dependencies=[qpath,ipath,HERE/'design.json']
    for field in ['C3','X3']:
        path=INPUT/'des'/f'Y3_DEEP_FIELDS_PHOTOM-SN-{field}-0000.parquet';dependencies.append(path)
        cat=pq.read_table(path,columns=['ID','RA','DEC','FLAGS','MASK_FLAGS','KNN_CLASS']).to_pandas();cat=cat[(cat.FLAGS==0)&(cat.MASK_FLAGS==0)&(cat.KNN_CLASS==1)]
        assert cat.ID.is_unique
        tree=cKDTree(xyz(cat.RA,cat.DEC));angle=2*max(FWHM.values())/3600*np.pi/180;radius=2*np.sin(angle/2)
        for _,r in q[q.field==field].sort_values('SNID').iterrows():
            ix=tree.query_ball_point(xyz([r.deep_RA],[r.deep_DEC])[0],radius);other=cat.iloc[ix];other=other[other.ID!=int(r.deep_ID)]
            jobs.append((dict(sn_id=r.SNID,host_id=int(r.deep_ID),field=field,ra=r.deep_RA,dec=r.deep_DEC),dict(ra=other.RA.to_numpy(),dec=other.DEC.to_numpy())))
    allrows=[];checks=[]
    with ProcessPoolExecutor(2,initializer=initialize) as pool:
        for i,(rows,local) in enumerate(pool.map(worker,jobs)):
            allrows+=rows;checks+=local
            if (i+1)%25==0:print(f'{i+1}/{len(jobs)} hosts',flush=True)
    frame=pd.DataFrame(allrows);output=WORK/'host-beam-information.csv';frame.to_csv(output,index=False)
    summary=[]
    for band,g in frame.groupby('band_um'):
        row={'band_um':int(band),'hosts':len(g),'numerically_unresolved':int(g.all_candidates_noise_inflation.isna().sum())}
        for key in ['nearest_separation_arcsec','neighbours_within_2FWHM','nearest_pair_radial_noise_inflation','all_candidates_noise_inflation','nearest_pair_native_orientation_min','nearest_pair_native_orientation_max','Gaussian_FWHM_for_pair_inflation2_arcsec']:
            row[key]={k:float(v) for k,v in zip(['p05','median','p95'],g[key].quantile([.05,.5,.95]))}
        row['all_candidates_inflation_gt10']=int((g.all_candidates_noise_inflation>10).sum());row['nearest_pair_inflation_gt2']=int((g.nearest_pair_radial_noise_inflation>2).sum());summary.append(row)
    assert max(checks)<1e-8
    dependencies += sorted(WORK.glob('*.fits.gz'))+sorted((INPUT/'help').glob('*NEST_image*fits'))
    result=dict(completed_utc=datetime.now(timezone.utc).isoformat(),seconds=time.perf_counter()-start,objects=265,rows=len(frame),beam_records=beam_records,summaries=summary,
        pair_formula_max_relative_difference=max(checks),independent_normal_matrix_comparisons=int(frame.independent_normal_matrix_relative_difference.notna().sum()),
        independent_normal_matrix_max_relative_difference=float(frame.independent_normal_matrix_relative_difference.max()),
        scope='Conditional linear resolution/information diagnostic with independent equal-variance pixels, empirical radial beam and unregularized signed nuisance amplitudes. No image flux or intrinsic host dust luminosity is inferred. Optical candidates need not emit; missing IR-only sources, source extension, coadd PSF, SED-dependent beam and correlated confusion remain.',
        dependencies_sha256={str(p.relative_to(ROOT)):sha(p) for p in dependencies},code_sha256=sha(Path(__file__)),output=str(output.relative_to(ROOT)),output_sha256=sha(output))
    (OUT/'information.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
