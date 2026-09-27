"""Independent native integration, optimizer and observed residual diagnostics."""
from __future__ import annotations
import copy,json,sys
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,sncosmo
from scipy.optimize import minimize
from scipy.special import log_ndtr
from acquire import ROOT,WORK,OUT,sha
sys.path.insert(0,str(ROOT))
from lib.flux_engine import Engine

def main():
    resultfile=OUT/'yse-flux-fits.json';results=json.loads(resultfile.read_text());model=Engine().snmodel
    checks=[]
    for z in [.0257,.15]:
        model.set(z=z,mwebv=.024,x0=.003,x1=-.8,c=.05,t0=0)
        for band in ['ps1::g','ps1::r','ps1::i','ps1::z','ztfg','ztfr']:
            bp=sncosmo.get_bandpass(band);wave=np.linspace(bp.minwave(),bp.maxwave(),int(np.ceil(bp.maxwave()-bp.minwave()))*2+1)
            trans=bp(wave);ab=3.631e-20*2.99792458e18/wave**2
            for phase in [-5.,0.,15.,30.]:
                time=phase*(1+z);f=model.flux(time,wave)
                direct=np.trapezoid(f*wave*trans,wave)/np.trapezoid(ab*wave*trans,wave)*10**(.4*27.5)
                native=float(model.bandflux(band,time,zp=27.5,zpsys='ab'))
                difference=float(-2.5*np.log10(direct/native));assert abs(difference)<.0002
                checks.append(dict(z=z,band=band,phase=phase,difference_mag=difference))
    completed=[r for r in results['records']if r['status']=='completed'];residuals=[];objective_checks=[]
    for row in completed:
        name=row['name'];fixture=WORK/'yse-fit-fixtures'/(name+'.npz');assert sha(fixture)==row['fixture_sha256'];f=np.load(fixture)
        model.set(z=float(f['z']),mwebv=float(f['mwebv']));p=row['fits']['positive_truncated']['parameters'];seed=float(f['seed'])
        x=np.array([-p['mB']*np.log(10)/2.5,p['x1'],p['c'],p['t0']-seed])
        def predicted(x):
            model.set(x0=np.exp(x[0]),x1=x[1],c=x[2],t0=seed+x[3]);return model.bandflux(f['band'],f['time'],zp=27.5,zpsys='ab')
        pred=predicted(x);res=(f['flux']-pred)/f['error']
        recomputed_chi2=float(res@res);assert abs(recomputed_chi2-row['fits']['positive_truncated']['gaussian_chi2'])<1e-5
        rows=[dict(band=b,time=float(t),phase=float((t-p['t0'])/(1+f['z'])),flux=float(v),error=float(s),predicted=float(mu),residual_sigma=float(r))for b,t,v,s,mu,r in zip(f['band'],f['time'],f['flux'],f['error'],pred,res)]
        residuals.append(dict(name=name,chi2=recomputed_chi2,dof=len(pred)-4,maximum_absolute_residual=float(np.max(abs(res))),per_band={b:dict(epochs=int(np.sum(f['band']==b)),chi2=float(np.sum(res[f['band']==b]**2)))for b in np.unique(f['band'])},epochs=rows))
        if name in ['2019tvv','2021oat','2019zfv','2021oaw']:
            def nll(x):
                m=predicted(x);return float(.5*np.sum(((f['flux']-m)/f['error'])**2+np.log(2*np.pi*f['error']**2))+np.sum(log_ndtr(m/f['error'])))
            opt=minimize(nll,x,method='BFGS',options={'gtol':1e-6,'maxiter':500})
            delta=float(opt.fun-row['fits']['positive_truncated']['nll']);assert abs(delta)<1e-3
            objective_checks.append(dict(name=name,optimizer='scipy BFGS',success=bool(opt.success),message=str(opt.message),objective_difference=delta,maximum_parameter_change=float(np.max(abs(opt.x-x)))))
    crosspath=WORK/'yse-host-crosswalk.csv';cross=pd.read_csv(crosspath);meta=[]
    for row in cross[cross.primary_candidate].itertuples():
        for line in (ROOT/row.lightcurve).read_text().splitlines():
            if not line.startswith('OBS:'):continue
            a=line.split();flux,error,mag,magerr=map(float,a[4:8]);snr=flux/error
            meta.append(dict(name=row.name,MJD=float(a[1]),band=a[2],FLUXCAL=flux,FLUXCALERR=error,MAG=mag,MAGERR=magerr,FLAG=a[8],SNR=snr,
               magnitude_minus_logflux=float(mag-(27.5-2.5*np.log10(flux)))))
    (WORK/'yse-residual-records.json').write_text(json.dumps(residuals,indent=2)+'\n')
    pd.DataFrame(meta).to_csv(WORK/'yse-magnitude-metadata-audit.csv',index=False)
    summary=dict(completed_utc=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),
      inputs_sha256={str(p.relative_to(ROOT)):sha(p)for p in [resultfile,crosspath]},outputs_sha256={str((WORK/n).relative_to(ROOT)):sha(WORK/n)for n in ['yse-residual-records.json','yse-magnitude-metadata-audit.csv']},
      native_quadrature_checks=checks,maximum_integration_difference_mag=max(abs(q['difference_mag'])for q in checks),optimizer_checks=objective_checks,
      candidates=results['candidates'],completed=len(completed),sampling_and_numerical_pass=sum(r['fits']['positive_truncated']['sampling_and_numerical_gate']for r in completed),
      nominal_residual_screen_pass=sum(r['fits']['positive_truncated']['nominal_residual_gate']for r in completed),approved_host_brightness_analysis=0,
      median_chi2_per_dof=float(np.median([r['chi2']/r['dof']for r in residuals])),
      magnitude_metadata=dict(rows=len(meta),negative_MAGERR=sum(r['MAGERR']<0 for r in meta),
          max_absolute_MAG_minus_logFLUX=max(abs(r['magnitude_minus_logflux'])for r in meta),
          examples=sorted(meta,key=lambda r:abs(r['magnitude_minus_logflux']),reverse=True)[:5],
          highSNR_gt5_max_absolute_MAG_minus_logFLUX=max(abs(r['magnitude_minus_logflux'])for r in meta if r['SNR']>5)),
      worst_residual_examples=sorted(residuals,key=lambda r:r['chi2']/r['dof'],reverse=True)[:4],
      conclusion='Native filter integration and objective arithmetic agree; released fluxes show large within-band residual excursions. Ignored template/model covariance, calibration, classification and photometric extraction must be investigated before a physical host association. No error inflation, posterior host coefficient or empirical B(z)constraint is reported. Nominal residual screen uses exploratory chi-square tail>=.01 and is not selection-conditioned predictive certification.')
    (OUT/'yse-fit-audit.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:summary[k]for k in ['completed','sampling_and_numerical_pass','nominal_residual_screen_pass','median_chi2_per_dof','maximum_integration_difference_mag']}))
if __name__=='__main__':main()
