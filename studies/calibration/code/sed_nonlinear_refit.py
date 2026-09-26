"""Frozen six-object nonlinear gate; no cosmology or SED coefficient fit."""
from pathlib import Path
import hashlib
import importlib.util
import json
import platform
import sys

import numpy as np
import pandas as pd
import scipy
from numpy.polynomial.legendre import legvander
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares
import sncosmo
from sncosmo.constants import HC_ERG_AA
from sncosmo.utils import integration_grid

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'runs/research_2026_09_26/astra_design/validation1020'
OUT = ROOT/'runs/research_2026_09_26/sed_nonlinear_refit'
spec=importlib.util.spec_from_file_location('sed',ROOT/'scripts/research_2026_09_26/sed_identification.py')
sed=importlib.util.module_from_spec(spec);spec.loader.exec_module(sed)
K=.4*np.log(10.)
D=np.array([1.,.16087,-3.1178,0.])
STEPS=np.array([1e-4,1e-3,1e-4,.01])
STARTS=[np.zeros(4),np.array([.05,.5,.03,2.])]
LOWER=np.array([-1.,-3.,-.3,-10.]);UPPER=-LOWER
INPUTS={}


def record(path):
    path=Path(path)
    INPUTS[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
    return path


class Engine:
    def __init__(self,nominal,model,bands,offsets,coeff,observer):
        self.model=model
        self.t=nominal['MJD'];self.band=nominal['band']
        self.x0,self.x1,self.c,self.t0=nominal['parameters_x0_x1_c_t0']
        self.z=float(nominal['zHEL'][0]);self.ebv=float(nominal['MWEBV'][0])
        local=coeff[:26]+sed.redshift_coordinate(self.z)*coeff[26:]
        self.grids=[]
        for name in sorted(set(self.band)):
            idx=np.flatnonzero(self.band==name);bp=bands[name]
            wave,dw=integration_grid(bp.minwave(),bp.maxwave(),5.)
            wv=legvander(2*np.log(wave/(1+self.z)/2000)/np.log(5.5)-1,8)
            shapes=[]
            for q in range(3):
                take=[j for j,(n,qq) in enumerate(sed.PAIRS) if qq==q]
                shapes.append(wv[:,[sed.PAIRS[j][0] for j in take]]@local[take])
            norm=10**(.4*27.5)/sncosmo.get_magsystem('ab').zpbandflux(bp)
            norm*=10**(-.4*(.27+offsets[name]))
            weights=wave*bp(wave)*dw/HC_ERG_AA*norm
            self.grids.append((idx,wave,weights,np.array(shapes),float(observer['griz'.index(name)])))

    def flux(self,theta,mode='nominal',derivative=False):
        self.model.set(z=self.z,t0=self.t0+theta[3],x0=self.x0*np.exp(-K*theta[0]),
                       x1=self.x1+theta[1],c=self.c+theta[2],mwebv=self.ebv,mwrv=3.1,hostebv=0.,hostrv=3.1)
        f=np.empty(len(self.t))
        for idx,wave,weights,shapes,observer in self.grids:
            spectral=self.model.flux(self.t[idx],wave)
            if mode=='sed':
                phase=(self.t[idx]-self.t0-theta[3])/(1+self.z)
                dm=legvander(np.tanh(phase/20),2)@shapes
                factor=-K*dm if derivative else np.exp(-K*dm)
                spectral=spectral*factor
            elif mode=='observer':
                spectral=spectral*(-K*observer if derivative else np.exp(-K*observer))
            elif mode!='nominal':
                raise ValueError(mode)
            f[idx]=spectral@weights
        return f

    def jac(self,theta,mode='nominal',factor=1.):
        cols=[-K*self.flux(theta,mode)]
        for j in range(1,4):
            step=np.eye(4)[j]*STEPS[j]*factor
            cols.append((self.flux(theta+step,mode)-self.flux(theta-step,mode))/(2*STEPS[j]*factor))
        return np.column_stack(cols)


def fit_pair(engine,target,cov,keep,mode):
    chol=np.linalg.cholesky(cov[np.ix_(keep,keep)])
    whiten=lambda v:solve_triangular(chol,v,lower=True)
    objective=lambda theta:whiten(engine.flux(theta,mode)[keep]-target[keep])
    jac=lambda theta:whiten(engine.jac(theta,mode)[keep])
    fits=[]
    for start in STARTS:
        opt=least_squares(objective,start.copy(),jac=jac,bounds=(LOWER,UPPER),
                          x_scale=np.array([.1,1.,.1,3.]),xtol=1e-11,ftol=1e-11,gtol=1e-9,max_nfev=300)
        jw=jac(opt.x);u,s,v=np.linalg.svd(jw,full_matrices=False)
        fisher_gradient=float(np.linalg.norm(u.T@opt.fun))
        half=whiten(engine.jac(opt.x,mode,.5)[keep])
        jac_error=float(np.linalg.norm(jw-half)/np.linalg.norm(half))
        boundary=float(np.min(np.minimum(opt.x-LOWER,UPPER-opt.x)/(UPPER-LOWER)))
        gate=bool(opt.success and fisher_gradient<1e-4 and s[-1]>s[0]*1e-10 and boundary>1e-4 and jac_error<1e-4)
        fits.append(dict(theta=opt.x.tolist(),chi2=float(opt.fun@opt.fun),success=bool(opt.success),status=int(opt.status),message=str(opt.message),nfev=int(opt.nfev),fisher_gradient_norm=fisher_gradient,jacobian_halfstep_relative_error=jac_error,boundary_fraction=boundary,singular_values=s.tolist(),individual_gate=gate))
    best=min(fits,key=lambda f:f['chi2'])
    diff=np.array(fits[0]['theta'])-np.array(fits[1]['theta'])
    agreement_chi=abs(fits[0]['chi2']-fits[1]['chi2'])
    agreement_mag=abs(D@diff)
    gate=all(f['individual_gate'] for f in fits) and agreement_chi<1e-6 and agreement_mag<1e-4
    return np.array(best['theta']),dict(gate=bool(gate),starts=fits,start_chi2_difference=float(agreement_chi),start_standardized_difference_mag=float(agreement_mag),best_chi2=best['chi2'])


def main():
    assert (OUT/'freeze.json').exists() and not (OUT/'result.json').exists()
    protocol=ROOT/'docs/research-2026-09-26/sed-nonlinear-refit.md'
    (OUT/'protocol.md').write_bytes(protocol.read_bytes())
    (OUT/'executed_source.py').write_bytes(Path(__file__).read_bytes())
    for path in [Path(__file__),OUT/'freeze.json',OUT/'cohort.csv',OUT/'protocol.md',ROOT/'scripts/research_2026_09_26/sed_identification.py',ROOT/'scripts/salt_dust_audit/flux_response.py']:
        record(path)
    frozen=json.loads((OUT/'freeze.json').read_text())
    for name,digest in frozen['input_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    coefficients=ROOT/'runs/research_2026_09_26/sed_identification/amplitude-constrained-coefficients.npz'
    with np.load(record(coefficients)) as data: coeff=data['broad_drift_discovery_rms0.05'][:,3]
    with np.load(record(BASE/'frozen-discovery-coefficients.npz')) as data: observer=data['gauge_griz']@data['basis_mean']
    model,bands,paths,zp=sed.fr.build_model()
    for path in paths:record(path)
    offsets={str(r['Filter Name'])[-1]:float(r['Primary Mag']) for r in zp}
    chosen=pd.read_csv(OUT/'cohort.csv',dtype={'CID':str})
    fit_details=[];responses=[];mean_checks=[];negative_rows=[];closure=[]
    prediction_arrays={}
    for row in chosen.itertuples():
        cid=row.CID
        with np.load(record(BASE/'objectives'/f'objective_{cid}.npz')) as data:
            names=['MJD','band','model_flux','data_flux','data_fluxerr','zHEL','MWEBV','parameters_x0_x1_c_t0','frozen_flux_covariance']
            nominal={k:data[k] for k in names}
        engine=Engine(nominal,model,bands,offsets,coeff,observer)
        native=engine.flux(np.zeros(4));official=nominal['model_flux'];cov=nominal['frozen_flux_covariance']
        full=np.arange(len(native));chol=np.linalg.cholesky(cov)
        cached=np.load(record(BASE/'analysis/objects'/f'{cid}.npz'))
        order=pd.DataFrame({'t':nominal['MJD'],'b':nominal['band']}).sort_values(['t','b'],kind='stable').index.to_numpy()
        mean_error=float(np.max(abs(native[order]-cached['native_flux_model'])/np.maximum(abs(cached['native_flux_model']),1e-30)))
        assert mean_error<1e-10,(cid,mean_error)
        mean_checks.append(dict(CID=cid,zHEL=engine.z,epochs=len(native),native_archive_relative_error=mean_error,native_minus_official_exactC_norm=float(np.linalg.norm(solve_triangular(chol,native-official,lower=True))),native_minus_official_max_quoted_sigma=float(np.max(abs(native-official)/nominal['data_fluxerr'])),native_nonpositive_epochs=int(np.sum(native<=0)),official_nonpositive_epochs=int(np.sum(official<=0))))
        bad=(native<=0)|(official<=0)
        for i in np.flatnonzero(bad):
            negative_rows.append(dict(CID=cid,band=nominal['band'][i],MJD=nominal['MJD'][i],phase=(nominal['MJD'][i]-engine.t0)/(1+engine.z),native_mean=native[i],official_mean=official[i],quoted_error=nominal['data_fluxerr'][i]))
        masks={'accepted_full':full}
        if np.any(bad):masks['positive_reference_subset']=np.flatnonzero(~bad)
        targets={'native_noiseless':native,'official_mean_sensitivity':official,'observed_flux':nominal['data_flux']}
        for mask,keep in masks.items():
            local_chol=np.linalg.cholesky(cov[np.ix_(keep,keep)])
            whiten=lambda v:solve_triangular(local_chol,v,lower=True)
            for target_name,target in targets.items():
                baseline,gate=fit_pair(engine,target,cov,keep,'nominal')
                fit_details.append(dict(CID=cid,target=target_name,mask=mask,mode='nominal',**gate))
                if target_name=='native_noiseless':
                    closure_gate=bool(gate['gate'] and gate['best_chi2']<1e-10 and np.max(abs(baseline))<1e-6)
                    closure.append(dict(CID=cid,mask=mask,gate=closure_gate,theta=baseline.tolist(),chi2=gate['best_chi2']))
                baseflux=engine.flux(baseline)
                jw=whiten(engine.jac(baseline)[keep]);inverse=np.linalg.pinv(jw,rcond=1e-10)
                for mode in ['observer','sed']:
                    actual,altgate=fit_pair(engine,target,cov,keep,mode)
                    fit_details.append(dict(CID=cid,target=target_name,mask=mask,mode=mode,**altgate))
                    shift=actual-baseline
                    delta=engine.flux(baseline,mode)-baseflux
                    infinitesimal=engine.flux(baseline,mode,derivative=True)
                    finite_linear=-inverse@whiten(delta[keep])
                    small_linear=-inverse@whiten(infinitesimal[keep])
                    record_response=dict(CID=cid,zHEL=engine.z,target=target_name,mask=mask,mode=mode,epochs=len(keep),paired_gate=bool(gate['gate'] and altgate['gate']),baseline_theta=baseline.tolist(),alternative_theta=actual.tolist(),actual_shift=shift.tolist(),finite_mean_tangent_shift=finite_linear.tolist(),infinitesimal_tangent_shift=small_linear.tolist(),nonlinear_standardized_mag=float(D@shift),finite_mean_tangent_standardized_mag=float(D@finite_linear),infinitesimal_tangent_standardized_mag=float(D@small_linear),nonlinear_minus_finite_tangent_mag=float(D@(shift-finite_linear)),nonlinear_minus_infinitesimal_tangent_mag=float(D@(shift-small_linear)),baseline_chi2=gate['best_chi2'],alternative_chi2=altgate['best_chi2'])
                    responses.append(record_response)
                    prediction_arrays[f'{cid}__{mask}__{target_name}__{mode}']=engine.flux(actual,mode)
                print(cid,mask,target_name,'done',flush=True)
    frame=pd.DataFrame(responses)
    frame.to_csv(OUT/'paired-responses.csv',index=False)
    pd.DataFrame(mean_checks).to_csv(OUT/'native-mean-checks.csv',index=False)
    pd.DataFrame(negative_rows,columns=['CID','band','MJD','phase','native_mean','official_mean','quoted_error']).to_csv(OUT/'nonpositive-reference-epochs.csv',index=False)
    np.savez_compressed(OUT/'fitted-model-means.npz',**prediction_arrays)
    contrasts={}
    low=set(chosen.CID.iloc[:2]);high=set(chosen.CID.iloc[-2:])
    for (target,mode),data in frame[frame['mask']=='accepted_full'].groupby(['target','mode']):
        contrasts[target+'__'+mode]={key:float(data.loc[data.CID.isin(high),key].mean()-data.loc[data.CID.isin(low),key].mean()) for key in ['nonlinear_standardized_mag','finite_mean_tangent_standardized_mag','infinitesimal_tangent_standardized_mag']}
    result={'scope':'Constructed six-object independent sncosmo nonlinear gate; exact exported SNANA covariance and accepted masks. No physical correction, SED population fit, retraining, selection/BBC or cosmology inference.',
            'cohort':chosen.to_dict('records'),'sed_coefficient_source':str(coefficients.relative_to(ROOT)),'sed_coefficient_key':'broad_drift_discovery_rms0.05 column3','frozen_observer_griz_mag':observer.tolist(),
            'all_fit_gates_pass':all(x['gate'] for x in fit_details),'all_native_noiseless_closure_gates_pass':all(x['gate'] for x in closure),
            'fit_details':fit_details,'closure':closure,'responses':responses,'mean_checks':mean_checks,
            'descriptive_six_object_highest2_minus_lowest2':contrasts,'environment':dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,sncosmo=sncosmo.__version__),
            'limits':['Two starts do not establish global optimization or calibrated interval coverage.','Fixed covariance/mask and finite family; phase follows fitted t0 but spectral training and covariance are not refitted.','Noiseless generator and fitter use the same independent sncosmo mean, not exact SNANA spectral code.','Observed-data branch is a perturbation diagnostic, not a fitted physical correction.','Six central-support objects cannot validate full1020 tails or establish full high255/low255 response.']}
    (OUT/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    outputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'}
    (OUT/'manifest.json').write_text(json.dumps({'inputs_sha256':INPUTS,'outputs_sha256':outputs,'argv':sys.argv},indent=2)+'\n')
    print(json.dumps({'all_fit_gates_pass':result['all_fit_gates_pass'],'all_native_noiseless_closure_gates_pass':result['all_native_noiseless_closure_gates_pass'],'contrasts':contrasts},indent=2))


if __name__=='__main__':main()
