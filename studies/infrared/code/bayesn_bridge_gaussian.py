"""Conjugate-only test of a direct optical-to-joint normalizing-constant bridge.

This file imports no BayeSN, JAX, NumPyro or cosmological model. Its Gaussian
draws are independent analytic samples, not a new physical posterior fit.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit, logsumexp
from scipy.stats import multivariate_normal

ROOT=Path(__file__).resolve().parents[3]


def log_bridge(log_likelihood_optical,log_likelihood_joint):
    """Return log Zjoint/Zoptical, retaining the full likelihood normalization.

    The bridge is asymptotically optimal for independent draws at fixed counts.
    Its consistency does not turn MCMC draws into independent observations.
    """
    x=np.asarray(log_likelihood_optical,dtype=float).ravel()
    y=np.asarray(log_likelihood_joint,dtype=float).ravel()
    assert len(x)>1 and len(y)>1 and np.isfinite(x).all() and np.isfinite(y).all()
    s0=len(x)/(len(x)+len(y));s1=1-s0
    # Center before solving, to preserve precision under an arbitrary density
    # constant. Adding that constant to both arrays must add it to log Z ratio.
    center=float(np.median(np.concatenate([x,y])))
    x=x-center;y=y-center;odds=np.log(s1/s0)
    def terms(t):
        a=expit(x-t+odds)/s1
        b=expit(t-y-odds)/s0
        return a,b
    def equation(t):
        a,b=terms(t);return float(a.mean()-b.mean())
    lo=float(min(x.min(),y.min())-50);hi=float(max(x.max(),y.max())+50)
    t=brentq(equation,lo,hi,xtol=1e-12,rtol=1e-14)
    a,b=terms(t)
    # Derivative of E0[a]-E1[b] with respect to log-ratio. This factor is
    # necessary when propagating uncertainty through the implicit estimator.
    b0=expit(t-x-odds)/s0;a1=expit(y-t+odds)/s1
    derivative=float(s0*np.mean(a*b0)+s1*np.mean(a1*b))
    assert derivative>0 and np.isfinite(derivative), 'No numerically resolved overlap.'
    iid_se=float(np.sqrt(np.var(a,ddof=1)/len(a)+np.var(b,ddof=1)/len(b))/derivative)
    return dict(log_ratio=float(t+center),iid_delta_standard_error=iid_se,
                estimating_equation_residual=equation(t),overlap_derivative=derivative)


def gaussian_case(rng,noise_scale,shift,repetitions,n=4000):
    # Six informed linear coordinates of a47-dimensional Gaussian optical
    # posterior; the other41 coordinates cancel exactly. Correlated noise is
    # included to test a genuinely joint vector normalization.
    loading=np.diag([1.,.8,.6,.4,.2,.1])
    covariance=noise_scale**2*(.8*np.eye(6)+.2*np.ones((6,6)))
    predictive=loading@loading.T+covariance
    observation=np.sqrt(np.diag(predictive))*shift
    mean=loading.T@np.linalg.solve(predictive,observation)
    postcov=np.eye(6)-loading.T@np.linalg.solve(predictive,loading)
    chol=np.linalg.cholesky(postcov)
    exact=float(multivariate_normal.logpdf(observation,mean=np.zeros(6),cov=predictive))
    def likelihood(points):
        return multivariate_normal.logpdf(observation-points@loading.T,mean=np.zeros(6),cov=covariance)
    direct=[];bridge=[];brute=[];reported=[];overlap=[]
    for _ in range(repetitions):
        p0=rng.normal(size=(16*n,6))
        p1=rng.normal(size=(n,6))@chol.T+mean
        ll0=likelihood(p0);ll1=likelihood(p1)
        estimate=log_bridge(ll0[:n],ll1)
        direct.append(float(logsumexp(ll0[:n])-np.log(n)))
        brute.append(float(logsumexp(ll0)-np.log(16*n)))
        bridge.append(estimate['log_ratio']);reported.append(estimate['iid_delta_standard_error'])
        overlap.append(estimate['overlap_derivative'])
    methods={}
    for name,values in [('direct_4000',direct),('bridge_4000_plus_4000',bridge),('direct_64000',brute)]:
        v=np.array(values)
        methods[name]=dict(mean_error=float(np.mean(v-exact)),RMSE=float(np.sqrt(np.mean((v-exact)**2))),
                           empirical_standard_deviation=float(v.std(ddof=1)))
    return dict(noise_scale=noise_scale,shift=shift,dimension=47,likelihood_rank=6,repetitions=repetitions,
                exact_log_predictive_density=exact,methods=methods,
                bridge_mean_reported_iid_error=float(np.mean(reported)),minimum_overlap_derivative=float(min(overlap)),
                replicate_values=dict(direct_4000=direct,bridge_4000_plus_4000=bridge,direct_64000=brute))


def run():
    rng=np.random.default_rng(8273992)
    constant=log_bridge(np.full(100,-17.),np.full(120,-17.))
    assert abs(constant['log_ratio']+17)<1e-12
    x=rng.normal(-10,2,4000);y=rng.normal(-7,2,4000)
    base=log_bridge(x,y);reverse=log_bridge(-y,-x);offset=log_bridge(x-100000,y-100000)
    assert abs(base['log_ratio']+reverse['log_ratio'])<1e-12
    assert abs(offset['log_ratio']-(base['log_ratio']-100000))<2e-11
    cases=[gaussian_case(rng,noise,shift,64) for noise,shift in [(.6,.5),(.3,.5),(.6,2.)]]
    # Deterministic conjugate quadrature controls the sign and normalization
    # independently of simulation error in the repeated Monte Carlo tests.
    from scipy.integrate import quad
    from scipy.stats import norm
    quadrature=[]
    for y,sigma in [(0.,.4),(2.,.7),(-3.,1.2)]:
        z=quad(lambda t:norm.pdf(t)*norm.pdf(y,loc=t,scale=sigma),-12,12,epsabs=1e-13,epsrel=1e-13)[0]
        exact=norm.logpdf(y,scale=np.hypot(1,sigma))
        assert abs(np.log(z)-exact)<1e-12
        quadrature.append(float(abs(np.log(z)-exact)))
    return dict(status='passed_conjugate_only_bridge_identity_checks',seed=8273992,cases=cases,
                constant_likelihood_error=abs(constant['log_ratio']+17),
                reciprocal_identity_error=abs(base['log_ratio']+reverse['log_ratio']),
                large_density_offset_error=abs(offset['log_ratio']-(base['log_ratio']-100000)),
                scalar_quadrature_log_errors=quadrature,
                physical_model_or_chain_calls=0,observed_or_synthetic_SN_payloads_read=False,
                limits='IID conjugateGaussian validation only. Auxiliary draws are exact and cheap here; no BayeSN speedup, posterior convergence, or real score precision is established.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();assert not args.output.exists()
    result=run()
    here=Path(__file__).resolve();design=here.with_name('bayesn-bridge-integration-design.json')
    result['source_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [here,design]}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='cases'}))
