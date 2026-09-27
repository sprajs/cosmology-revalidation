"""Signed host spectral-band measurements with mask and covariance ledgers."""
from __future__ import annotations
from datetime import datetime,timezone
import json,re
from pathlib import Path
import numpy as np,pandas as pd
from astropy.io import fits
from scipy.linalg import toeplitz
from acquire import ROOT,WORK,OUT,sha


def band_operator(wave,z,valid,windows):
    step=float(np.median(np.diff(wave)));rest=wave/(1+z);half=step/(1+z)/2
    weights=[];coverage=[]
    for lo,hi in windows:
        overlap=np.maximum(0,np.minimum(rest+half,hi)-np.maximum(rest-half,lo))*valid
        coverage.append(float(overlap.sum()/(hi-lo)))
        weights.append(overlap/overlap.sum()*(rest/4000)**2 if overlap.sum()>0 else np.zeros(len(wave)))
    return np.asarray(weights),np.asarray(coverage)


def propagated_covariance(operator,variance,rho=0.):
    good=np.any(operator!=0,axis=0);indices=np.where(good)[0]
    root=operator[:,good]*np.sqrt(variance[good])[None,:]
    if rho==0:return root@root.T
    correlation=rho**np.abs(indices[:,None]-indices[None,:])
    return root@correlation@root.T


def ratio_profile(flux,covariance,grid):
    """Profile positive band amplitude at each candidate red/blue ratio.

    Signed observations and their complete2x2 working covariance enter directly;
    no Gaussian approximation to flux division is used.
    """
    inverse=np.linalg.inv(covariance);direction=np.column_stack([np.ones(len(grid)),grid])
    denom=np.einsum('ni,ij,nj->n',direction,inverse,direction)
    amplitude=np.maximum(0,(direction@inverse@flux)/denom)
    residual=flux-direction*amplitude[:,None]
    chi2=np.einsum('ni,ij,nj->n',residual,inverse,residual)
    delta=chi2-chi2.min();keep=delta<=3.841458820694124
    return dict(mle=float(grid[np.argmin(chi2)]),lower=float(grid[keep][0]),upper=float(grid[keep][-1]),
                lower_grid_limited=bool(keep[0]),upper_grid_limited=bool(keep[-1]),
                amplitude_at_mle=float(amplitude[np.argmin(chi2)]),minimum_chi2=float(chi2.min())),delta


def main():
    design_path=Path(__file__).with_name('ozdes-design.json');design=json.loads(design_path.read_text())
    acquisition_path=OUT/'ozdes-acquisition.json';acquisition=json.loads(acquisition_path.read_text())
    cross_path=ROOT/'.work/unified-cosmology/survey-selection/followup/dovekie-ozdes-clean.csv'
    cross=pd.read_csv(cross_path,dtype={'SNID':str,'oz_OzDES_ID':str}).drop_duplicates('oz_OzDES_ID').set_index('oz_OzDES_ID')
    names=list(design['rest_windows_A']);windows=list(design['rest_windows_A'].values())
    ib,ir=names.index('Dn_blue'),names.index('Dn_red');grid=np.linspace(0,4,1601)
    rows=[];vectors=[];covs=[];coverages=[];profiles=[];kernel_ids={};validation=[];failures=[]
    folder=WORK/'ozdes-band-kernels';folder.mkdir(exist_ok=True)
    for number,source in enumerate(acquisition['records']):
        if source['status']!='downloaded':failures.append(source);continue
        target=source['target'];meta=cross.loc[target]
        path=ROOT/source['path'];assert sha(path)==source['sha256']
        with fits.open(path) as hdus:
            h=hdus[0].header;flux=np.asarray(hdus[0].data,float);variance=np.asarray(hdus['VARIANCE'].data,float);mask=np.asarray(hdus['BADPIX'].data)
            z=float(h['Z']);wave=h['CRVAL1']+(np.arange(len(flux))+1-h['CRPIX1'])*h['CDELT1']
            assert h.get('CUNIT1','').strip().lower()=='angstrom' and h['SOURCE']==target and abs(z-meta.oz_z)<1e-5
            dichroics=set()
            used_exposures=0;flags=[]
            for epoch in hdus[3::3]:
                if epoch.header.get('INCOADD',True):
                    used_exposures+=1;flags.append(epoch.header.get('QC','unknown'))
                    dichroics.update(map(int,re.findall(r'[xX](5700|6700)',str(epoch.header.get('BLUSENSF','')))))
        valid=(mask==0)&np.isfinite(flux)&np.isfinite(variance)&(variance>0)
        safe_flux=np.where(valid,flux,0.);safe_variance=np.where(valid,variance,0.)
        operator,coverage=band_operator(wave,z,valid,windows)
        value=operator@safe_flux
        covariance=np.asarray([propagated_covariance(operator,safe_variance,rho) for rho in [0.,.25,.5]])
        eligible=bool(np.all(coverage[[ib,ir]]>=.8) and np.all(np.diag(covariance[0])[[ib,ir]]>0))
        if eligible:
            prof,delta=ratio_profile(value[[ib,ir]],covariance[0][np.ix_([ib,ir],[ib,ir])],grid)
            prof25,_=ratio_profile(value[[ib,ir]],covariance[1][np.ix_([ib,ir],[ib,ir])],grid)
            prof50,_=ratio_profile(value[[ib,ir]],covariance[2][np.ix_([ib,ir],[ib,ir])],grid)
        else:prof=prof25=prof50=None;delta=np.full(len(grid),np.nan)
        sensitivity={}
        for shift in [-2.,2.]:
            shifted,_=band_operator(wave+shift,z,valid,windows);sensitivity[str(shift)]=(shifted@safe_flux-value).tolist()
        known_dichroics=sorted(dichroics) or [5700,6700]
        crosses=any(3850*(1+z)<=d+150 and 4100*(1+z)>=d-150 for d in known_dichroics)
        row=dict(target=target,z=z,source=str(path.relative_to(ROOT)),good_pixels=int(valid.sum()),
                 Dn_blue=float(value[ib]),Dn_red=float(value[ir]),var_blue=float(covariance[0,ib,ib]),
                 cov_blue_red=float(covariance[0,ib,ir]),var_red=float(covariance[0,ir,ir]),
                 blue_coverage=float(coverage[ib]),red_coverage=float(coverage[ir]),profile_eligible=eligible,
                 blue_SNR=float(value[ib]/np.sqrt(covariance[0,ib,ib])) if covariance[0,ib,ib]>0 else None,
                 red_SNR=float(value[ir]/np.sqrt(covariance[0,ir,ir])) if covariance[0,ir,ir]>0 else None,
                 dichroic_overlap=crosses,dichroics=known_dichroics,coadd_exposures=used_exposures,
                 native_QC_counts=pd.Series(flags).value_counts().to_dict(),native_DO_HELIO=h.get('DO_HELIO'),
                 ratio_profile=prof,ratio_profile_rho025=prof25,ratio_profile_rho05=prof50,
                 wavelength_shift_band_changes=sensitivity)
        rows.append(row);vectors.append(value);covs.append(covariance);coverages.append(coverage);profiles.append(delta)
        indices=np.where(np.any(operator!=0,axis=0))[0]
        kernel=folder/(target+'.npz')
        np.savez_compressed(kernel,pixel_index=indices,wavelength_A=wave[indices],operator=operator[:,indices],
                            pixel_flux=flux[indices],pixel_variance=variance[indices],z=z)
        kernel_ids[str(kernel.relative_to(ROOT))]=sha(kernel)
        if number%100==0:
            # Independent explicit diagonal matrix calculation and overlap check.
            direct=operator[:,indices]@np.diag(variance[indices])@operator[:,indices].T
            error=float(np.max(np.abs(direct-covariance[0]))/max(1e-30,np.max(np.abs(direct))))
            assert error<1e-12
            validation.append(dict(target=target,relative_covariance_difference=error,
                                   overlapping_window_covariance=float(covariance[0,ir,names.index('Hd_blue')])))
            print('measured',number,flush=True)
    scalar=[]
    for row in rows:
        r={k:v for k,v in row.items() if not isinstance(v,(dict,list))}
        if row['ratio_profile']:r.update({'ratio_'+k:v for k,v in row['ratio_profile'].items()})
        scalar.append(r)
    pd.DataFrame(scalar).to_csv(WORK/'ozdes-band-ledger.csv',index=False)
    (WORK/'ozdes-band-records.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(WORK/'ozdes-band-likelihood.npz',target=np.asarray([r['target'] for r in rows]),
         band_names=np.asarray(names),flux=np.asarray(vectors),covariance=np.asarray(covs),coverage=np.asarray(coverages),
         rho=np.asarray([0.,.25,.5]),ratio_grid=grid,ratio_delta_chi2=np.asarray(profiles))
    redshift=[]
    for lo,hi in [(0,.3),(.3,.5),(.5,.7),(.7,.9),(.9,1.3)]:
        subset=[r for r in rows if lo<=r['z']<hi];eligible=[r for r in subset if r['profile_eligible']]
        identified=[r for r in eligible if not r['ratio_profile']['lower_grid_limited'] and not r['ratio_profile']['upper_grid_limited']]
        redshift.append(dict(zmin=lo,zmax=hi,hosts=len(subset),profile_eligible=len(eligible),
             finite_two_sided_95_ratio_profiles=len(identified),dichroic_overlap=sum(r['dichroic_overlap'] for r in subset),
             blue_nonpositive=sum(r['Dn_blue']<=0 for r in eligible),red_nonpositive=sum(r['Dn_red']<=0 for r in eligible),
             median_blue_SNR=float(np.median([r['blue_SNR'] for r in eligible])) if eligible else None))
    result=dict(completed_utc=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),design_sha256=sha(design_path),
       inputs_sha256={str(p.relative_to(ROOT)):sha(p) for p in [acquisition_path,cross_path]},
       outputs_sha256={str((WORK/n).relative_to(ROOT)):sha(WORK/n) for n in ['ozdes-band-ledger.csv','ozdes-band-records.json','ozdes-band-likelihood.npz']},
       kernel_sha256=kernel_ids,hosts=len(rows),failed_acquisitions=failures,band_order=names,covariance_rhos=[0,.25,.5],redshift_support=redshift,
       independent_covariance_checks=validation,
       scope='Measured signed count-spectrum windows with covariance under stated pixel-noise assumptions. Ratio likelihood profiles are conditional on local relative response; grid-limited support is not a physical bound. No ages, SN attenuation, selection-complete host population or empirical B(z)correction inferred.',
       physical_identification='Unknown spectrophotometric response, coadd resampling covariance, central aperture, selected redshift success and host/progenitor mapping prevent treating these indices as independent absolute ages. Shared target enters host likelihood once.')
    (OUT/'ozdes-band-interface.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'hosts':len(rows),'redshift_support':redshift}))

if __name__=='__main__':main()
