"""Refit public YSE selected fluxes, preserving truncation and failed gates."""
from __future__ import annotations
import copy,gzip,json,os,sys
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd,sncosmo
from astropy.table import Table
from iminuit import Minuit
from scipy.special import log_ndtr
from scipy.stats import chi2
from acquire import ROOT,WORK,OUT,sha
sys.path.insert(0,str(ROOT))
from lib.flux_engine import Engine,ASSETS

BANDS={'g':'ps1::g','r':'ps1::r','i':'ps1::i','z':'ps1::z','X':'ztfg','Y':'ztfr'}
ENGINE=None

def initialize():
    global ENGINE
    ENGINE=Engine()
    for name in BANDS.values():sncosmo.get_bandpass(name)

def observations(path):
    rows=[]
    for line in path.read_text().splitlines():
        if line.startswith('OBS:'):
            a=line.split();rows.append(dict(time=float(a[1]),band=BANDS[a[2]],flux=float(a[4]),fluxerr=float(a[5]),zp=27.5,zpsys='ab',native_magerr=float(a[7])))
    return Table(rows=rows)

def native_fit(row):
    global ENGINE
    if ENGINE is None:initialize()
    name=row['name'];result=dict(name=name,status='attempted',fits={})
    path=ROOT/row['lightcurve'];data=observations(path)
    model=copy.deepcopy(ENGINE.snmodel);z=row['host_z_hel']
    model.set(z=z,mwebv=row['mwebv'])
    try:
        valid=np.isfinite(data['flux'])&np.isfinite(data['fluxerr'])&(data['fluxerr']>0)&(data['flux']>0)
        assert np.all(valid),'Unexpected released invalid/negative flux'
        candidates=data[np.isin(data['band'],['ps1::r','ztfr'])]
        if len(candidates)==0:candidates=data
        peak=float(candidates['time'][np.argmax(candidates['flux'])])
        initial=data[(data['time']>peak-20*(1+z))&(data['time']<peak+50*(1+z))]
        first,model=sncosmo.fit_lc(initial,model,['t0','x0','x1','c'],
          bounds={'t0':(peak-25,peak+25),'x1':(-4,4),'c':(-.5,.5)},modelcov=False,maxcall=5000,warn=False)
        if not first.success or model.get('x0')<=0:raise ValueError('Initial native fit did not converge with positive amplitude')
        seed=float(model.get('t0'));phase=(data['time']-seed)/(1+z)
        data=data[(phase>=-10)&(phase<=40)]
        if len(data)<10:raise ValueError('Fewer than10 epochs in frozen phase window')
        initial=np.array([np.log(model.get('x0')),model.get('x1'),model.get('c'),0.])
        band,time,flux,error=[np.asarray(data[key]) for key in ['band','time','flux','fluxerr']]
        def predicted(theta):
            model.set(x0=np.exp(theta[0]),x1=theta[1],c=theta[2],t0=seed+theta[3])
            return model.bandflux(band,time,zp=27.5,zpsys='ab')
        for mode in ['positive_truncated','gaussian_selected_sensitivity']:
            def objective(lnx0,x1,c,dt):
                f=predicted([lnx0,x1,c,dt]);nll=.5*np.sum(((flux-f)/error)**2+np.log(2*np.pi*error**2))
                if mode=='positive_truncated':nll+=np.sum(log_ndtr(f/error))
                return float(nll)
            fit=Minuit(objective,*initial,name=('lnx0','x1','c','dt'));fit.errordef=.5
            fit.limits=[(initial[0]-5,initial[0]+5),(-4,4),(-.5,.5),(-25,25)]
            fit.migrad(ncall=6000);fit.hesse()
            theta=np.array(fit.values);cov=np.array(fit.covariance) if fit.covariance is not None else None
            estimate=predicted(theta)
            phase=(time-seed-theta[3])/(1+z)
            phase_gate=bool(np.any(phase<0)&np.any(phase>10))
            goodbands=sum(np.sum((band==b)&(flux/error>=5))>=2 for b in np.unique(band))
            eigen=None if cov is None else np.linalg.eigvalsh(cov)
            boundary=any(abs(theta[i]-edge)<1e-3 for i,limits in enumerate(fit.limits) for edge in limits)
            numerical=bool(fit.valid and cov is not None and eigen.min()>0 and not boundary)
            quality=bool(numerical and abs(theta[1])<3 and abs(theta[2])<.3 and np.sqrt(cov[1,1])<1 and np.sqrt(cov[2,2])<.1 and phase_gate and goodbands>=2)
            jac=np.diag([-2.5/np.log(10),1,1,1])
            raw_chi2=float(np.sum(((flux-estimate)/error)**2))
            nominal_tail=float(chi2.sf(raw_chi2,len(data)-4))
            entry=dict(converged=bool(fit.valid),numerical_valid=numerical,sampling_and_numerical_gate=quality,
              nominal_chi2_tail=nominal_tail,nominal_residual_gate=bool(nominal_tail>=.01),
              analysis_gate=False,
              analysis_gate_reason='Model/template/calibration and censored-epoch likelihood not validated; nominal residual screening is exploratory and is not predictive validation.',
              parameters=dict(mB=float(-2.5/np.log(10)*theta[0]),x1=float(theta[1]),c=float(theta[2]),t0=float(seed+theta[3])),
              parameter_order=['mB_arbitrary_zero','x1','c','t0_MJD'],covariance=None if cov is None else (jac@cov@jac).tolist(),
              minimum_covariance_eigenvalue=None if cov is None else float(eigen.min()),boundary=bool(boundary),
              phase_gate=phase_gate,bands_with_two_SNR5=int(goodbands),epochs=len(data),nll=float(fit.fval),
              gaussian_chi2=raw_chi2,positive_truncation_logprob=float(np.sum(log_ndtr(estimate/error))),
              fluxerr_magerr_identity_max_difference=float(np.max(abs(data['native_magerr']-2.5/np.log(10)*error/flux))),
              min_observed_SNR=float(np.min(flux/error)),optimizer_edm=float(fit.fmin.edm),
              retained_phase_min=float(phase.min()),retained_phase_max=float(phase.max()))
            result['fits'][mode]=entry
        fixture=WORK/'yse-fit-fixtures';fixture.mkdir(exist_ok=True)
        np.savez_compressed(fixture/(name+'.npz'),band=band,time=time,flux=flux,error=error,z=z,mwebv=row['mwebv'],seed=seed)
        result['fixture_sha256']=sha(fixture/(name+'.npz'));result['status']='completed'
    except Exception as exc:result.update(status='failed',failure=type(exc).__name__+': '+str(exc))
    return result

def main():
    source=WORK/'yse-host-crosswalk.csv';table=pd.read_csv(source)
    chosen=table[table.primary_candidate]
    initialize()
    filters=WORK/'yse-passbands';filters.mkdir(exist_ok=True)
    for name in BANDS.values():
        bp=sncosmo.get_bandpass(name);np.savetxt(filters/(name.replace(':','_')+'.txt'),np.c_[bp.wave,bp.trans])
    results=[]
    with ProcessPoolExecutor(max_workers=4,initializer=initialize) as pool:
        for result in pool.map(native_fit,chosen.to_dict('records')):
            results.append(result);print(result['name'],result['status'],flush=True)
    record=dict(completed_utc=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),design_sha256=sha(Path(__file__).with_name('yse-design.json')),
      input_sha256={str(source.relative_to(ROOT)):sha(source)},
      dependencies_sha256={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'lib/flux_engine.py',*ASSETS.glob('*'),*filters.glob('*.txt')] if p.is_file()},
      candidates=len(chosen),completed=sum(r['status']=='completed' for r in results),records=results,
      scope='New conditional selected-flux fits, arbitrary common mB normalization. Only positive-flux truncation normalized; epochMAGERR/quality and event/hostselection, template/calibration covariance and actualredshift uncertainty remain unresolved. No certified bias-free distances.')
    (OUT/'yse-flux-fits.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    print('YSE fits complete',record['completed'],'/',record['candidates'])

if __name__=='__main__':main()
