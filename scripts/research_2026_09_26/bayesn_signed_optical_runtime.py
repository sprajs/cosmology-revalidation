"""Capped synthetic-only NUTS runtime forecast for both normalized D priors."""
from __future__ import annotations
import os
os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:2])
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
os.environ.setdefault('XLA_FLAGS','--xla_cpu_multi_thread_eigen=false')
import json,time
import numpy as np
from scipy.integrate import quad
from scipy.stats import norm
import jax
import jax.numpy as jnp
from jax.scipy.special import log_ndtr
import numpyro
import numpyro.distributions as dist
from numpyro.distributions import constraints
from numpyro.infer import MCMC,NUTS
from bayesn_signed_optical_gate import OUT,make_refined,sha,save

GRAY=.088

def broad_logpdf(d,lo=20.,hi=50.,s=GRAY):
    """Log of Uniform(lo,hi) convolved with Normal(0,s), reflected tails."""
    a=jnp.where(d<(lo+hi)/2,(d-lo)/s,(hi-d)/s)
    b=jnp.where(d<(lo+hi)/2,(d-hi)/s,(lo-d)/s)
    la=log_ndtr(a);lb=log_ndtr(b)
    return la+jnp.log(-jnp.expm1(lb-la))-jnp.log(hi-lo)

class BroadD(dist.Distribution):
    support=constraints.real
    arg_constraints={}
    def __init__(self,batch_shape=(2,),validate_args=None):
        super().__init__(batch_shape=batch_shape,event_shape=(),validate_args=validate_args)
    def sample(self,key,sample_shape=()):
        k1,k2=jax.random.split(key)
        shape=tuple(sample_shape)+self.batch_shape
        return jax.random.uniform(k1,shape,minval=20.,maxval=50.)+GRAY*jax.random.normal(k2,shape)
    def log_prob(self,value):return broad_logpdf(value)

def check_prior():
    grid=np.r_[np.linspace(19.7,20.3,21),np.linspace(34,36,9),np.linspace(49.7,50.3,21),[15.,55.]]
    got=np.asarray(broad_logpdf(jnp.array(grid)))
    ref=[]
    for d in grid:
        # Independent adaptive integration over mu, with Gaussian gray kernel.
        q=quad(lambda mu:norm.pdf(d,mu,GRAY)/30,20,50,epsabs=1e-12,points=[max(20,min(50,d))])[0]
        ref.append(np.log(q) if q>0 else -np.inf)
    ref=np.asarray(ref)
    central=np.isfinite(ref)
    err=float(np.max(np.abs(got[central]-ref[central])))
    Z=quad(lambda d:float(np.exp(broad_logpdf(jnp.array(d)))),15,55,epsabs=1e-12,points=[20,35,50])[0]
    return dict(max_logpdf_vs_independent_mu_quadrature=err,normalization_15_55=float(Z),tested_points=int(central.sum()),
                extreme_tail_points=int((~central).sum()),tail_caveat='Independent direct mu quadrature underflows beyond support; stable reflected logCDF retained.')

def main():
    start=time.monotonic()
    protocol=json.loads((OUT/'runtime-protocol.json').read_text())
    assert sha(__file__)==protocol['source_sha256']
    assert sha(OUT/'forward-algebra-result.json')==protocol['forward_result_sha256']
    prior=check_prior();save('runtime-prior-check.json',prior)
    if prior['max_logpdf_vs_independent_mu_quadrature']>=1e-6 or abs(prior['normalization_15_55']-1)>=1e-9:
        raise AssertionError('Broad D prior algebra gate failed')
    optical=json.loads((OUT/'optical-payload.json').read_text());objs=optical['objects'];rows=optical['rows'];n=2;nobs=15
    t=np.zeros((nobs,n));errors=np.ones((nobs,n));mask=np.zeros((nobs,n));bands=np.full((nobs,n),'NULL_BAND',object)
    for j,o in enumerate(objs):
        rr=[r for r in rows if r['cid']==o['cid']]
        for i,r in enumerate(rr):t[i,j]=r['trigger_rest_time'];errors[i,j]=r['error'];mask[i,j]=1;bands[i,j]='RAISIN_DES_'+r['band']
    Refined,_=make_refined();Refined.audit_spectrum_bins=2400
    model=Refined(load_model='M20_model',num_devices=1,filter_yaml=str(OUT/'release-filter-config.yaml'))
    weights=model._calculate_band_weights(jnp.array([o['zHEL'] for o in objs]),jnp.array([o['MWEBV'] for o in objs]))
    bi=jnp.array([[model.band_dict[b] for b in row] for row in bands]);tj=jnp.array(t);maskj=jnp.array(mask)
    sig=jnp.array(errors);mu=jnp.array([o['mu_LCDM'] for o in objs]);sigma=jnp.array([np.hypot(o['sigma_external'],GRAY) for o in objs])
    def flux(p):
        phase=tj-p[:,46][None,:]
        jt=model.J_t_map(phase.flatten(order='F'),model.tau_knots,model.KD_t).reshape((nobs,n,6),order='F').transpose(1,2,0)
        hs=jnp.array([19+jnp.floor(phase),19+jnp.ceil(phase),jnp.remainder(phase,1)])
        eps=(model.L_Sigma@p[:,4:46].T).T.reshape((n,7,6),order='F')
        eps=jnp.zeros((n,9,6)).at[:,1:-1,:].set(eps)
        return model.get_flux_batch(model.M0,p[:,3],p[:,1],model.W0,model.W1,eps,p[:,0],p[:,2],bi,maskj,jt,hs,weights)
    fid=np.zeros((n,47));fid[:,0]=np.asarray(mu);fid[:,1]=.3;fid[:,2]=3.1;fid[:,46]=5.
    synthetic_y=jax.lax.stop_gradient(flux(jnp.array(fid)))
    def make_model(arm):
        def toy_model():
            D=numpyro.sample('D',dist.Normal(mu,sigma).to_event(1) if arm=='LCDM' else BroadD().to_event(1))
            AV=numpyro.sample('AV',dist.Exponential(jnp.ones(n)/.329).to_event(1))
            RV=numpyro.sample('RV',dist.Uniform(jnp.full(n,1.2),jnp.full(n,6.)).to_event(1))
            theta=numpyro.sample('theta',dist.Normal(jnp.zeros(n),jnp.ones(n)).to_event(1))
            eps=numpyro.sample('epsilon_white',dist.Normal(jnp.zeros((n,42)),jnp.ones((n,42))).to_event(2))
            tau=numpyro.sample('tau',dist.Uniform(jnp.full(n,-10.),jnp.full(n,20.)).to_event(1))
            p=jnp.concatenate([D[:,None],AV[:,None],RV[:,None],theta[:,None],eps,tau[:,None]],axis=1)
            logp=dist.Normal(flux(p),sig).log_prob(synthetic_y)
            numpyro.factor('synthetic_optical_normalized',jnp.sum(jnp.where(maskj>0,logp,0.)))
        return toy_model
    results=[]
    for i,arm in enumerate(('LCDM','broad')):
        ts=time.monotonic()
        mcmc=MCMC(NUTS(make_model(arm)),num_warmup=80,num_samples=20,num_chains=1,progress_bar=False)
        mcmc.run(jax.random.PRNGKey(20260926+i),extra_fields=('num_steps','diverging','accept_prob'))
        fields=mcmc.get_extra_fields();steps=np.asarray(fields['num_steps']);acc=np.asarray(fields['accept_prob'])
        row=dict(arm=arm,warmup=80,discarded_draws=20,seconds=time.monotonic()-ts,
                 postwarmup_num_steps_min=int(steps.min()),postwarmup_num_steps_median=float(np.median(steps)),
                 postwarmup_num_steps_p90=float(np.quantile(steps,.9)),postwarmup_num_steps_max=int(steps.max()),
                 postwarmup_num_steps_sum=int(steps.sum()),divergences=int(np.asarray(fields['diverging']).sum()),
                 accept_prob_mean=float(acc.mean()),per_transition_seconds=(time.monotonic()-ts)/100)
        results.append(row);save('runtime-partial.json',dict(prior_check=prior,arms=results,elapsed_seconds=time.monotonic()-start))
        print(json.dumps(row),flush=True)
        if time.monotonic()-start>180:raise TimeoutError('Cumulative 180-second runtime benchmark exceeded')
    result=dict(scope='Synthetic-only runtime forecast, not chain convergence or observed likelihood.',prior_check=prior,arms=results,
                elapsed_seconds=time.monotonic()-start,forecast_1000warmup_1000draws_single_chain_seconds={
                    x['arm']:2000*x['per_transition_seconds'] for x in results},
                caveat='Timings include compilation and short adaptation; full 47-latent chain cost and mixing remain uncertain; no full chain authorized.')
    save('runtime-result.json',result);print(json.dumps(result),flush=True)
if __name__=='__main__':main()
