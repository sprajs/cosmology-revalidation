"""Observed exposure repeatability and local continuum-cancelled signed features."""
from __future__ import annotations
import hashlib,json
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from astropy.io import fits
from scipy.optimize import minimize_scalar
from scipy.stats import chi2
from acquire import ROOT,WORK,OUT,sha
from ozdes_bands import band_operator,propagated_covariance

WINDOWS=[[3850,3950],[4000,4100]]
HDWINDOWS=[[4041.6,4079.75],[4083.5,4122.25],[4128.5,4161]]

def common_ratio(flux,cov):
    inverse=np.linalg.inv(cov)
    def evaluate(r,details=False):
        direction=np.array([1.,r]);denom=np.einsum('i,nij,j->n',direction,inverse,direction)
        amp=np.maximum(0,np.einsum('i,nij,nj->n',direction,inverse,flux)/denom)
        res=flux-amp[:,None]*direction
        val=float(np.einsum('ni,nij,nj->',res,inverse,res))
        return (val,amp) if details else val
    grid=np.linspace(0,4,81);values=np.array([evaluate(r) for r in grid]);i=int(values.argmin())
    opt=minimize_scalar(evaluate,bounds=(grid[max(0,i-1)],grid[min(80,i+1)]),method='bounded',options={'xatol':1e-9})
    options=[(float(opt.fun),float(opt.x)),(values[0],0.),(values[-1],4.)];value,r=min(options)
    _,amp=evaluate(r,True)
    return dict(ratio=r,chi2=value,amplitudes=amp,grid_boundary=bool(r<1e-6 or r>4-1e-6))

def arrays(hdus,index,z):
    h=hdus[index].header;f=np.asarray(hdus[index].data,float);v=np.asarray(hdus[index+1].data,float);m=np.asarray(hdus[index+2].data)
    wave=h['CRVAL1']+(np.arange(len(f))+1-h['CRPIX1'])*h['CDELT1']
    valid=(m==0)&np.isfinite(f)&np.isfinite(v)&(v>0)
    return wave,np.where(valid,f,0),np.where(valid,v,0),valid

def main():
    acquisition=OUT/'ozdes-acquisition.json';design=Path(__file__).with_name('ozdes-repeat-design.json')
    records=json.loads(acquisition.read_text())['records'];rows=[];features=[];kernels={};repeatfixtures={}
    folder=WORK/'ozdes-repeat';folder.mkdir(exist_ok=True)
    for number,source in enumerate(records):
        if source['status']!='downloaded':continue
        target=source['target'];path=ROOT/source['path'];assert sha(path)==source['sha256']
        with fits.open(path) as hdus:
            z=float(hdus[0].header['Z']);wave,flux,var,valid=arrays(hdus,0,z)
            wfnu,coverage=band_operator(wave,z,valid,HDWINDOWS)
            rest=wave/(1+z);w=wfnu/(rest/4000)[None,:]**2
            means=w@rest;t=(means[1]-means[0])/(means[2]-means[0])
            transform=np.array([[1-t,0,t],[1-t,-1,t]])
            op=transform@w;observed=op@flux;cov=propagated_covariance(op,var)
            eligible=bool(np.all(coverage>=.8) and np.all(np.diag(cov)>0))
            moment0=float(op[1].sum());moment1=float(op[1]@rest)
            curvature=float(op[1]@((rest-4102.875)/50.)**2)
            feature=dict(target=target,z=z,continuum=float(observed[0]),absorption_contrast=float(observed[1]),
               variance_continuum=float(cov[0,0]),covariance_continuum_contrast=float(cov[0,1]),variance_contrast=float(cov[1,1]),
               coverage=coverage.tolist(),eligible=eligible,linear_moment0=moment0,linear_moment1_A=moment1,
               unit_quadratic_curvature_response=curvature)
            features.append(feature)
            pix=np.where(np.any(op!=0,axis=0))[0];kernel=folder/(target+'-Hdelta.npz')
            np.savez_compressed(kernel,pixel_index=pix,rest_wavelength_A=rest[pix],operator=op[:,pix],flux=flux[pix],variance=var[pix])
            kernels[str(kernel.relative_to(ROOT))]=sha(kernel)
            ff=[];cc=[];metadata=[];invalid=0
            for index in range(3,len(hdus),3):
                h=hdus[index].header
                if not h.get('INCOADD',True):continue
                ew,ef,ev,good=arrays(hdus,index,z);operator,cover=band_operator(ew,z,good,WINDOWS)
                covariance=propagated_covariance(operator,ev)
                if not(np.all(cover>=.8) and np.all(np.linalg.eigvalsh(covariance)>0)):invalid+=1;continue
                ff.append(operator@ef);cc.append(covariance);metadata.append(dict(index=index,QC=h.get('QC'),FLXSCALE=h.get('FLXSCALE')))
        row=dict(target=target,z=z,usable_exposures=len(ff),unusable_exposures=invalid)
        if len(ff)>=2:
            ff=np.asarray(ff);cc=np.asarray(cc);fit=common_ratio(ff,cc)
            row.update(ratio=fit['ratio'],chi2=fit['chi2'],dof=len(ff)-1,nominal_tail=float(chi2.sf(fit['chi2'],len(ff)-1)),grid_boundary=fit['grid_boundary'],nonpositive_amplitudes=int(np.sum(fit['amplitudes']==0)))
            file=folder/(target+'-exposures.npz');np.savez_compressed(file,flux=ff,covariance=cc,fit_amplitudes=fit['amplitudes'])
            row['exposure_metadata']=metadata;repeatfixtures[str(file.relative_to(ROOT))]=sha(file)
        rows.append(row)
        if number%100==0:print('exposure host',number,flush=True)
    rng=np.random.default_rng(927083);boot=[]
    selected=sorted([r for r in rows if r['usable_exposures']>=3],key=lambda r:hashlib.sha256(r['target'].encode()).hexdigest())[:40]
    for row in selected:
        f=np.load(folder/(row['target']+'-exposures.npz'));mean=f['fit_amplitudes'][:,None]*np.array([1.,row['ratio']]);cholesky=np.linalg.cholesky(f['covariance']);stats=[]
        for _ in range(200):
            draw=mean+np.einsum('nij,nj->ni',cholesky,rng.normal(size=mean.shape));stats.append(common_ratio(draw,f['covariance'])['chi2'])
        pvalue=(1+int(np.sum(np.asarray(stats)>=row['chi2'])))/201
        boot.append(dict(target=row['target'],replicates=200,observed_chi2=row['chi2'],tail_probability=pvalue,simulated_statistic_quantiles=np.quantile(stats,[.025,.5,.975]).tolist(),grid_boundary=row['grid_boundary']))
    (WORK/'ozdes-repeat-records.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
    (WORK/'ozdes-Hdelta-records.json').write_text(json.dumps(features,indent=2,allow_nan=False)+'\n')
    fitrows=[r for r in rows if r['usable_exposures']>=2];eligible=[r for r in features if r['eligible']]
    summary=dict(completed_utc=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),design_sha256=sha(design),
      dependencies_sha256={str(p.relative_to(ROOT)):sha(p) for p in [acquisition,Path(__file__).with_name('ozdes_bands.py')]},
      outputs_sha256={str((WORK/n).relative_to(ROOT)):sha(WORK/n) for n in ['ozdes-repeat-records.json','ozdes-Hdelta-records.json']},
      kernel_sha256=kernels,repeat_fixture_sha256=repeatfixtures,hosts=len(rows),hosts_with_repeated_usable_exposures=len(fitrows),
      total_usable_individual_exposures=sum(r['usable_exposures'] for r in rows),
      nominal_chi2_tail_below001=sum(r['nominal_tail']<.01 for r in fitrows),ratio_grid_boundary=sum(r['grid_boundary'] for r in fitrows),
      median_chi2_per_dof=float(np.median([r['chi2']/r['dof'] for r in fitrows])),
      bootstrap=boot,bootstrap_tail_below001=sum(r['tail_probability']<.01 for r in boot),
      Hdelta=dict(eligible=len(eligible),nonpositive_continuum=sum(r['continuum']<=0 for r in eligible),negative_absorption_contrast=sum(r['absorption_contrast']<0 for r in eligible),
        maximum_constant_continuum_residual=max(abs(r['linear_moment0'])for r in eligible),maximum_linear_continuum_residual_A=max(abs(r['linear_moment1_A'])for r in eligible)),
      scope='Observed differential exposure consistency, conditional on working pixel covariance. Does not calibrate shared response. Hdelta signed contrast cancels a locally linear count continuum; no physical age, dust or cosmological correction inferred. Nominal chi2 tails are diagnostic; selected40 parametric tests condition on fitted Gaussian model.')
    (OUT/'ozdes-repeatability.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['kernel_sha256','repeat_fixture_sha256','bootstrap']}))
if __name__=='__main__':main()
