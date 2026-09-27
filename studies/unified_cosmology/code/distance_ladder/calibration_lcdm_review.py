"""Independent joint-GLS and scalar-quadrature review of anchored SN-only LCDM.

No production projection or cosmological-integrator function is imported.
Analysis uses general LU solves of the full symmetric released covariance.
The established ReleasedCalibration class is imported only for fixed-point
density comparisons after the independent GLS calculation.
"""
from __future__ import annotations
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
from numpy.polynomial.legendre import leggauss
from scipy import linalg, special, stats
from scipy.integrate import quad
from scipy.optimize import toms748

ROOT = Path(__file__).resolve().parents[4]
WORK = ROOT/'.work/unified-cosmology/calibration-interface'
K = np.log(10.)/5
HREF = 70.
PROBS = np.array([.025, .16, .5, .84, .975])
PINS = {'Pantheon+SH0ES.dat': '1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8',
        'Pantheon+SH0ES_STAT+SYS.cov': 'abf806d966485e64afdb359c87bffc0ecc00d05eff0a31ced66f247385df0fdc'}
RULES = {'distance_order': 80, 'Omega_GL_order': 384,
         'direct_Omega_cases': [.1, .3, .5, .8], 'direct_H0_cases': [50., 65., 75., 90.],
         'H0_mean_sd_quantile_max_abs': 1e-5, 'Omega_mean_sd_quantile_max_abs': 1e-7,
         'H0_Omega_covariance_max_abs': 1e-7, 'loglike_max_abs': 1e-7,
         'distance_max_relative': 1e-10, 'adaptive_normalization_relative': 1e-8,
         'adaptive_CDF_max_abs': 1e-8, 'conditional_moment_relative': 1e-9,
         'initial_design': 'General full-C LU, X=[1,-I_noncal], eta=5log10(H0/70), flat-H0 exponential tilt, direct H0 scalar quadrature. No CMB or MCMC.'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


class JointGLS:
    def __init__(self, order=80, h_bounds=(50.,90.), omega_bounds=(.01,.99)):
        for name, digest in PINS.items():
            assert sha(WORK/name)==digest
        data=pd.read_csv(WORK/'Pantheon+SH0ES.dat',sep=r'\s+')
        mask=(data.zHD.to_numpy()>.01)|data.IS_CALIBRATOR.to_numpy().astype(bool)
        self.data=data.loc[mask].reset_index(drop=True)
        self.cal=self.data.IS_CALIBRATOR.to_numpy().astype(bool)
        self.z=self.data.zHD.to_numpy()[~self.cal]
        self.zhel=self.data.zHEL.to_numpy()[~self.cal]
        with (WORK/'Pantheon+SH0ES_STAT+SYS.cov').open() as stream:
            n=int(stream.readline());raw=np.loadtxt(stream).reshape(n,n)
        C=(raw+raw.T)/2
        self.C=C[np.ix_(mask,mask)]
        self.LU=linalg.lu_factor(self.C)
        self.logdetC=float(np.log(abs(np.diag(self.LU[0]))).sum())
        self.X=np.column_stack((np.ones(len(self.data)),-(~self.cal).astype(float)))
        self.PX=linalg.lu_solve(self.LU,self.X)
        self.F=self.X.T@self.PX
        self.V=linalg.solve(self.F,np.eye(2),assume_a='gen')
        assert np.linalg.eigvalsh((self.F+self.F.T)/2).min()>0
        self.logdetF=float(np.linalg.slogdet(self.F)[1])
        self.s=float(np.sqrt(self.V[1,1]))
        x,w=leggauss(order);self.x=(x+1)/2;self.w=w/2
        self.h_bounds=np.array(h_bounds);self.omega_bounds=np.array(omega_bounds)
        self.eta_bounds=np.log(self.h_bounds/HREF)/K

    def distance(self, omega, H0=HREF):
        # Independent 80-node rule, distinct from producer's 64/128-node rules.
        iz=self.z*np.sum(self.w/np.sqrt(omega*(1+self.z[:,None]*self.x)**3+1-omega),axis=1)
        return 299792.458/H0*iz/(1+self.z)

    def response(self, omega):
        theory=self.data.CEPH_DIST.to_numpy().copy()
        theory[~self.cal]=5*np.log10((1+self.z)*(1+self.zhel)*self.distance(omega))+25
        return self.data.m_b_corr.to_numpy()-theory

    @lru_cache(maxsize=20000)
    def joint(self, omega):
        y=self.response(omega)
        # General joint solve; no one-dimensional projected precision is used.
        fit=linalg.solve(self.F,self.PX.T@y,assume_a='gen')
        residual=y-self.X@fit
        chi2=float(residual@linalg.lu_solve(self.LU,residual))
        # Multiplying by dH0/deta=K*70*exp(K*eta) tilts BOTH coordinates.
        tilted=fit+K*self.V[:,1]
        mass=self.interval_mass(tilted[1])
        assert mass>0 and np.isfinite(mass)
        log_integral=(-.5*(chi2+self.logdetC+self.logdetF+(len(y)-2)*np.log(2*np.pi))
                      +np.log(K*HREF/np.diff(self.h_bounds)[0])+K*fit[1]
                      +.5*K*K*self.V[1,1]+np.log(mass))
        return {'mean_gls':fit,'mean_flat_H0':tilted,'minimum_chi2':chi2,
                'mass':float(mass),'log_integral':float(log_integral)}

    def interval_mass(self, eta_mean, upper=None):
        hi=self.eta_bounds[1] if upper is None else upper
        lo=self.eta_bounds[0]
        # In this declared bounded application both endpoints surround the
        # conditional means. Fail instead of silently cancelling tiny tails.
        return float(stats.norm.cdf((hi-eta_mean)/self.s)-stats.norm.cdf((lo-eta_mean)/self.s))

    def h_moment(self, omega, power):
        r=self.joint(float(omega));m=r['mean_flat_H0'][1];t=power*K
        return float(HREF**power*np.exp(t*m+.5*t*t*self.V[1,1])
                     *self.interval_mass(m+t*self.V[1,1])/r['mass'])

    def h_cdf(self, omega, value):
        if value<=self.h_bounds[0]:return 0.
        if value>=self.h_bounds[1]:return 1.
        r=self.joint(float(omega));m=r['mean_flat_H0'][1]
        return self.interval_mass(m,np.log(value/HREF)/K)/r['mass']

    def direct_flat_M_loglike(self, omega, H0):
        """Conditional M-only GLS at a fixed H0, independent of the joint integral."""
        y=self.response(omega)-self.X[:,1]*(np.log(H0/HREF)/K)
        a=self.F[0,0];M=float(self.PX[:,0]@y/a)
        r=y-self.X[:,0]*M
        q=float(r@linalg.lu_solve(self.LU,r))
        return -.5*(q+self.logdetC+np.log(a)+(len(y)-1)*np.log(2*np.pi))


def posterior(model):
    x,w=leggauss(RULES['Omega_GL_order']);lower,upper=model.omega_bounds
    omegas=(x+1)*(upper-lower)/2+lower
    logs=np.array([model.joint(float(o))['log_integral'] for o in omegas])+np.log(w/2)
    norm=float(special.logsumexp(logs));weights=np.exp(logs-norm)
    mh=np.array([model.h_moment(o,1) for o in omegas]);sh=np.array([model.h_moment(o,2) for o in omegas])
    mean_o=float(weights@omegas);mean_h=float(weights@mh)
    sd_o=float(np.sqrt(weights@((omegas-mean_o)**2)))
    sd_h=float(np.sqrt(weights@sh-mean_h**2))
    # Independent adaptive normalization and Omega CDF; fixed uniform breakpoints
    # are numerical integration support, not data-selected cuts or posterior priors.
    offset=max(model.joint(float(o))['log_integral'] for o in omegas)
    split=np.linspace(lower,upper,9)[1:-1]
    def integral(fn=lambda o:1., a=lower, b=upper):
        return quad(lambda o: np.exp(model.joint(float(o))['log_integral']-offset)*fn(o),
                    a,b,epsabs=1e-12,epsrel=3e-11,points=[s for s in split if a<s<b],limit=160)
    ad_norm,ad_err=integral()
    omega_q=[toms748(lambda v:integral(b=v)[0]/ad_norm-p,lower,upper,xtol=1e-11) for p in PROBS]
    h_q=[toms748(lambda v:float(weights@np.array([model.h_cdf(o,v) for o in omegas]))-p,
                *model.h_bounds,xtol=1e-9) for p in PROBS]
    adaptive={'log_normalization_abs_error':abs(np.log(ad_norm)+offset-np.log(upper-lower)-norm),
              'relative_error_estimate':ad_err/ad_norm,
              'Omega_mean_abs_error':abs(integral(lambda o:o)[0]/ad_norm-mean_o),
              'H0_mean_abs_error':abs(integral(lambda o:model.h_moment(o,1))[0]/ad_norm-mean_h),
              'H0_CDF_abs_errors':[abs(integral(lambda o:model.h_cdf(o,q))[0]/ad_norm-p) for q,p in zip(h_q,PROBS)]}
    return {'Omega_m':{'mean':mean_o,'sd':sd_o,'quantiles':omega_q},
            'H0_km_s_Mpc':{'mean':mean_h,'sd':sd_h,'quantiles':h_q},
            'H0_Omega_m_covariance':float(weights@((omegas-mean_o)*(mh-mean_h))),
            'quantile_probabilities':PROBS.tolist()},adaptive


def numerical_checks(model):
    from calibration_interface import ReleasedCalibration
    reference=ReleasedCalibration()
    rows=[];d_err=[];moment_errors=[]
    for om in RULES['direct_Omega_cases']:
        rr=model.joint(om)
        for z in [.01,.1,.5,1.,float(model.z.max())]:
            a=quad(lambda zz:1/np.sqrt(om*(1+zz)**3+1-om),0,z,epsabs=2e-13,epsrel=2e-13)[0]
            b=z*np.sum(model.w/np.sqrt(om*(1+z*model.x)**3+1-om))
            d_err.append(abs(a-b)/a)
        for H0 in RULES['direct_H0_cases']:
            own=model.direct_flat_M_loglike(om,H0)
            official=reference.evaluate(model.distance(om,H0))['loglike']
            rows.append({'Omega_m':om,'H0':H0,'independent_loglike':own,'interface_loglike':official,'abs_error':abs(own-official)})
        # Integrate the actual M-profiled density directly in H0. No transformed
        # eta Jacobian or analytic lognormal moment is used in this comparator.
        offset=model.direct_flat_M_loglike(om,HREF)
        def integ(power):
            return quad(lambda H: np.exp(model.direct_flat_M_loglike(om,H)-offset)*H**power,
                        *model.h_bounds,epsabs=1e-10,epsrel=1e-10,points=[60.,70.,80.],limit=100)[0]
        I=integ(0)
        moment_errors.append({'Omega_m':om,'log_integral_abs_error':abs(np.log(I/np.diff(model.h_bounds)[0])+offset-rr['log_integral']),
                              'first_moment_relative_error':abs(integ(1)/I/model.h_moment(om,1)-1),
                              'second_moment_relative_error':abs(integ(2)/I/model.h_moment(om,2)-1)})
    return {'fixed_point_loglike':rows,'distance_relative_max_error':max(d_err),
            'direct_H0_integration':moment_errors,'joint_covariance_M_eta':model.V.tolist(),
            'joint_GLS_information':model.F.tolist(),
            'normal_equation_identity_max_error':float(abs(model.F@model.V-np.eye(2)).max())}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--producer',type=Path,required=True)
    ap.add_argument('--output',type=Path,default=ROOT/'studies/unified_cosmology/results/distance_ladder/calibration-lcdm-review.json')
    ap.add_argument('--execution-design',type=Path,help='Optional fresh or identical immutable design under .work; default is keyed by producer and review source identities.')
    args=ap.parse_args();start=time.monotonic()
    producer=json.loads(args.producer.read_text())
    assert producer['status']=='qualified_conditional_SN_only_LCDM_quadrature'
    assert producer['priors']['Omega_m']==[.01,.99] and producer['priors']['H0_km_s_Mpc']==[50.,90.]
    bindings={**producer['source_sha256'],**producer['input_sha256'],relative(args.producer):sha(args.producer),relative(Path(__file__)):sha(__file__)}
    for path,digest in bindings.items():assert sha(ROOT/path)==digest
    design={'rules':RULES,'producer_sha256':sha(args.producer),'source_sha256':sha(__file__),'no_CMB_or_MCMC':True}
    design_key=hashlib.sha256(json.dumps(design,sort_keys=True,separators=(',',':')).encode()).hexdigest()[:20]
    frozen=(args.execution_design or ROOT/f'.work/unified-cosmology/calibration-lcdm-review/execution-design-{design_key}.json').resolve()
    assert frozen.is_relative_to(ROOT/'.work'), 'Execution designs belong in ignored .work.'
    if frozen.exists():assert json.loads(frozen.read_text())==design
    else:write(frozen,design)
    model=JointGLS(order=RULES['distance_order'])
    checks=numerical_checks(model)
    result,adaptive=posterior(model)
    differences={}
    for name in ['H0_km_s_Mpc','Omega_m']:
        a=result[name];b=producer['posterior'][name]
        differences[name]={k:abs(a[k]-b[k]) for k in ['mean','sd']}
        differences[name]['quantiles_abs_difference']=abs(np.array(a['quantiles'])-b['quantiles']).tolist()
    differences['H0_Omega_m_covariance']=abs(result['H0_Omega_m_covariance']-producer['posterior']['H0_Omega_m_covariance'])
    failed=[]
    for name,key in [('H0_km_s_Mpc','H0_mean_sd_quantile_max_abs'),('Omega_m','Omega_mean_sd_quantile_max_abs')]:
        d=differences[name]
        if max(d['mean'],d['sd'],max(d['quantiles_abs_difference']))>RULES[key]:failed.append(key)
    if differences['H0_Omega_m_covariance']>RULES['H0_Omega_covariance_max_abs']:failed.append('joint_covariance')
    if max(r['abs_error'] for r in checks['fixed_point_loglike'])>RULES['loglike_max_abs']:failed.append('density_normalization')
    if checks['distance_relative_max_error']>RULES['distance_max_relative']:failed.append('distance_quadrature')
    if adaptive['log_normalization_abs_error']>RULES['adaptive_normalization_relative']:failed.append('adaptive_normalization')
    if max(adaptive['H0_CDF_abs_errors'])>RULES['adaptive_CDF_max_abs']:failed.append('adaptive_CDF')
    if max(r['log_integral_abs_error'] for r in checks['direct_H0_integration'])>RULES['loglike_max_abs']:failed.append('direct_H0_integral')
    if max(r[k] for r in checks['direct_H0_integration'] for k in ['first_moment_relative_error','second_moment_relative_error'])>RULES['conditional_moment_relative']:failed.append('direct_H0_moments')
    for path,digest in bindings.items():assert sha(ROOT/path)==digest
    report={'status':'passed_independent_joint_GLS_review' if not failed else 'failed_independent_joint_GLS_review',
            'failed_gates':failed,'independent_posterior':result,'producer_comparison':differences,
            'numerical_checks':checks,'adaptive_integrals':adaptive,'rules':RULES,
            'input_source_sha256':bindings,'execution_design_path':relative(frozen),'execution_design_sha256':sha(frozen),
            'seconds':time.monotonic()-start,'CMB_calls':0,'MCMC_steps':0,
            'interpretation':'Independent numerical solution of the same released observations and priors; not an independent observational measurement or a verification of the original covariance construction.'}
    write(args.output,report)
    print(json.dumps({k:report[k] for k in ['status','failed_gates','independent_posterior','producer_comparison','seconds']},indent=2))
    if failed:raise SystemExit(1)


if __name__=='__main__':
    main()
