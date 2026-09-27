"""Independent scalar-distance integration and synthetic isolation checks.

No SED fitting, observed flux, or posterior sampler is invoked.
"""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import sys
from model import GRAY,HERE,ROOT,broad_logpdf,configure_cpu,dump,sha,source_hashes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'.work/infrared/heldout-synthetic-validation.json')
    args = parser.parse_args()
    configure_cpu()
    import numpy as np
    from scipy.integrate import quad,simpson
    from scipy.special import log_ndtr,logsumexp
    from scipy.stats import norm,qmc
    from synthetic import worker_spec
    from score import monte_carlo_score
    source = source_hashes()
    # The broad-D density is independently integrated over the original mu.
    prior_checks = []
    for d in np.r_[np.linspace(19.8,20.2,11),35.,np.linspace(49.8,50.2,11)]:
        expected = quad(lambda mu:norm.pdf(d,loc=mu,scale=GRAY)/30,20,50,
                        points=[max(20,min(50,d))],epsabs=1e-14)[0]
        actual = float(broad_logpdf(d))
        prior_checks.append(abs(actual-np.log(expected)))
    assert max(prior_checks)<1e-6
    K = .4*np.log(10)
    f,y,s = np.array([2.,4.,1.]),np.array([1.8,3.7,-.2]),np.array([.5,.8,1.])
    fn,yn,sn = np.array([1.5,2.]),np.array([1.6,1.7]),np.array([.4,.5])
    def llike(d,flux,obs,err):
        amplitude = np.exp(-K*(np.asarray(d)-35))
        return np.sum(norm.logpdf(obs,amplitude[...,None]*flux,err),axis=-1)
    def scipy_broad(d):
        d = np.asarray(d)
        # Reflect upper tail before taking differences; not the JAX implementation.
        a = np.minimum(d-20,50-d)/GRAY
        b = a-30/GRAY
        la,lb = log_ndtr(a),log_ndtr(b)
        return la+np.log(-np.expm1(lb-la))-np.log(30)
    integrals = []
    for arm in ('LCDM','broad'):
        logprior = (lambda d:norm.logpdf(d,35,np.hypot(.03,GRAY))) if arm=='LCDM' else scipy_broad
        lo,hi = (35-12*np.hypot(.03,GRAY),35+12*np.hypot(.03,GRAY)) if arm=='LCDM' else (20-12*GRAY,50+12*GRAY)
        points = sorted(set([lo,hi]+[x for x in (20,30,34,34.5,35,35.5,36,40,50) if lo<x<hi]))
        def integrand(d,with_nir):
            return np.exp(llike(d,f,y,s)+logprior(d)+(llike(d,fn,yn,sn) if with_nir else 0))
        opt = quad(lambda d:integrand(d,False),lo,hi,points=points[1:-1],epsabs=1e-13,epsrel=1e-10)[0]
        joint = quad(lambda d:integrand(d,True),lo,hi,points=points[1:-1],epsabs=1e-13,epsrel=1e-10)[0]
        exact = float(np.log(joint/opt))
        grid = np.linspace(lo,hi,240001)
        dense = float(np.log(simpson(integrand(grid,True),x=grid)/simpson(integrand(grid,False),x=grid)))
        # Independent normalized importance sampling of the conditional D posterior.
        # Scrambled Sobol uniforms -> Gaussian + exact-prior mixture. The prior
        # component explicitly covers the nonzero zero-amplitude likelihood tail.
        sampled = []
        for seed in range(273501,273509):
            u = qmc.Sobol(d=2,scramble=True,seed=seed).random_base2(18)
            local = 35+.5*norm.ppf(u[:,0])
            prior_draw = (35+np.hypot(.03,GRAY)*norm.ppf(u[:,0]) if arm=='LCDM'
                          else 20+30*u[:,0]+GRAY*norm.ppf(u[:,1]))
            # Integrate each mixture component separately with its normalized
            # probability, avoiding a discontinuous random component indicator.
            numerator,denominator = [],[]
            for probability,d in ((.9,local),(.1,prior_draw)):
                logq = np.logaddexp(np.log(.9)+norm.logpdf(d,35,.5),np.log(.1)+logprior(d))
                lw = llike(d,f,y,s)+logprior(d)-logq
                numerator.append(np.log(probability)+logsumexp(lw+llike(d,fn,yn,sn))-np.log(len(d)))
                denominator.append(np.log(probability)+logsumexp(lw)-np.log(len(d)))
            sampled.append(float(logsumexp(numerator)-logsumexp(denominator)))
        error = abs(float(np.mean(sampled))-exact)
        mcse = float(np.std(sampled,ddof=1)/np.sqrt(len(sampled)))
        assert abs(dense-exact)<1e-6
        assert error < 1e-6 and error < max(1e-7,5*mcse)
        # A product of separately marginalized epoch likelihoods is a different object.
        wrong = 0.
        for k in range(2):
            value = quad(lambda d:np.exp(llike(d,f,y,s)+logprior(d)+llike(d,fn[k:k+1],yn[k:k+1],sn[k:k+1])),lo,hi,points=points[1:-1],epsabs=1e-13)[0]
            wrong += np.log(value/opt)
        assert abs(wrong-exact)>1e-4
        integrals.append(dict(arm=arm,log_joint_over_optical=exact,dense_grid_difference=dense-exact,
                              conditional_D_sample_difference=float(np.mean(sampled))-exact,scramble_MCSE=mcse,
                              incorrectly_independent_NIR_epoch_score_difference=float(wrong-exact),
                              grids=240001,QMC_draws=2*8*2**18,
                              QMC_proposal='Stratified 0.9 Normal(35,.5^2)+0.1 exact proper D prior; no tail truncation'))
    # A held-out payload cannot enter the worker specification/initialization.
    prepared = dict(cases=[dict(optical_sha256='opt',nir_sha256='nir')],archive='archive',source_sha256=source,
                    design_sha256='design',filter_config_sha256='filters')
    changed = copy.deepcopy(prepared); changed['cases'][0]['nir_sha256']='arbitrary replacement'
    assert worker_spec(Path('test'),0,'LCDM',0,prepared)==worker_spec(Path('test'),0,'LCDM',0,changed)
    # Check shared-vector score against direct normalized summation, constant and dependent cases.
    score_checks = []
    rng = np.random.default_rng(273510)
    for index in range(12):
        values = np.full((4,1000),-2.) if index==0 else rng.normal(-5,1,(4,1000))
        actual = monte_carlo_score(values)
        expected = float(logsumexp(values)-np.log(values.size))
        assert abs(actual['log_predictive_density']-expected)<1e-12
        score_checks.append(abs(actual['log_predictive_density']-expected))
    extreme = np.full((4,1000),-2.); extreme[0] = -10000.
    stressed = monte_carlo_score(extreme)
    assert not stressed['score_support_passed'] and stressed['chain_scaled_mean_underflow']==[True,False,False,False]
    assert np.isfinite(stressed['chain_log_predictive_density']).all()
    json.dumps(stressed,allow_nan=False)
    assert source_hashes()==source
    result = dict(status='passed',scope='Synthetic algebra and worker-isolation checks only; no posterior, native SED call or observed outcome.',
                  source_sha256=source,broad_prior_max_log_density_error=max(prior_checks),scalar_D_checks=integrals,
                  NIR_swap_worker_spec_invariant=True,score_checks=12,max_score_sum_error=max(score_checks),
                  whole_chain_underflow_retained_as_support_failure=True,
                  tail_bounds='D integration encloses 12 Gaussian-gray standard deviations outside proper mu support; omitted prior tail below 4e-33.',
                  limitation='These conditional one-dimensional checks establish the integration/sampling identity; they do not establish high-dimensional NUTS mixing.')
    dump(args.output,result)
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    main()
