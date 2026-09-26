"""Capped synthetic-only BayeSN gradient and tiny warmup timing benchmark."""
from __future__ import annotations
import os
os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:2])
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
os.environ.setdefault('XLA_FLAGS','--xla_cpu_multi_thread_eigen=false')
import hashlib,json,time
from pathlib import Path
import numpy as np
import jax
import jax.numpy as jnp
import numpyro
import numpyro.distributions as dist
from numpyro.infer import MCMC,NUTS
from bayesn_signed_optical_gate import ROOT,OUT,make_refined,sha,save

def main():
    start=time.monotonic()
    protocol=json.loads((OUT/'benchmark-protocol.json').read_text())
    assert protocol['source_sha256']==sha(__file__)
    assert protocol['forward_result_sha256']==sha(OUT/'forward-algebra-result.json')
    optical=json.loads((OUT/'optical-payload.json').read_text())
    objs=optical['objects'];rows=optical['rows'];n=2;nobs=15
    t=np.zeros((nobs,n));errors=np.ones((nobs,n));mask=np.zeros((nobs,n));bands=np.full((nobs,n),'NULL_BAND',object)
    for j,o in enumerate(objs):
        rr=[r for r in rows if r['cid']==o['cid']]
        for i,r in enumerate(rr):t[i,j]=r['trigger_rest_time'];errors[i,j]=r['error'];mask[i,j]=1;bands[i,j]='RAISIN_DES_'+r['band']
    Refined,_=make_refined();Refined.audit_spectrum_bins=2400
    model=Refined(load_model='M20_model',num_devices=1,filter_yaml=str(OUT/'release-filter-config.yaml'))
    weights=model._calculate_band_weights(jnp.array([o['zHEL'] for o in objs]),jnp.array([o['MWEBV'] for o in objs]))
    bi=jnp.array([[model.band_dict[b] for b in row] for row in bands]);tj=jnp.array(t);maskj=jnp.array(mask)
    sig=jnp.array(errors);mu=jnp.array([o['mu_LCDM'] for o in objs]);sigma=jnp.array([np.hypot(o['sigma_external'],.088) for o in objs])
    def flux(p):
        phase=tj-p[:,46][None,:]
        jt=model.J_t_map(phase.flatten(order='F'),model.tau_knots,model.KD_t).reshape((nobs,n,6),order='F').transpose(1,2,0)
        hs=jnp.array([19+jnp.floor(phase),19+jnp.ceil(phase),jnp.remainder(phase,1)])
        eps=(model.L_Sigma@p[:,4:46].T).T.reshape((n,7,6),order='F')
        eps=jnp.zeros((n,9,6)).at[:,1:-1,:].set(eps)
        return model.get_flux_batch(model.M0,p[:,3],p[:,1],model.W0,model.W1,eps,p[:,0],p[:,2],bi,maskj,jt,hs,weights)
    fid=np.zeros((n,47));fid[:,0]=np.asarray(mu);fid[:,1]=.3;fid[:,2]=3.1;fid[:,46]=5.
    synthetic_y=jax.lax.stop_gradient(flux(jnp.array(fid)))
    def synthetic_logp(p):return -.5*jnp.sum((((synthetic_y-flux(p))/sig)*maskj)**2)
    g=jax.jit(jax.grad(synthetic_logp))
    p=jnp.array(fid)
    g(p).block_until_ready();warm_compile=time.monotonic()-start
    ts=time.monotonic()
    for _ in range(100):g(p).block_until_ready()
    gradient_seconds=time.monotonic()-ts
    result=dict(scope='Synthetic-only, no observed likelihood or NIR payload',compile_seconds=warm_compile,
                gradients=100,gradient_seconds=gradient_seconds,seconds_per_gradient=gradient_seconds/100)
    save('benchmark-gradient-partial.json',result)
    print(json.dumps(result),flush=True)
    # Five warmup transitions and one discarded synthetic draw time the full latent model.
    def toy_model():
        D=numpyro.sample('D',dist.Normal(mu,sigma).to_event(1))
        AV=numpyro.sample('AV',dist.Exponential(jnp.ones(n)/.329).to_event(1))
        RV=numpyro.sample('RV',dist.Uniform(jnp.full(n,1.2),jnp.full(n,6.)).to_event(1))
        theta=numpyro.sample('theta',dist.Normal(jnp.zeros(n),jnp.ones(n)).to_event(1))
        eps=numpyro.sample('epsilon_white',dist.Normal(jnp.zeros((n,42)),jnp.ones((n,42))).to_event(2))
        tau=numpyro.sample('tau',dist.Uniform(jnp.full(n,-10.),jnp.full(n,20.)).to_event(1))
        x=jnp.concatenate([D[:,None],AV[:,None],RV[:,None],theta[:,None],eps,tau[:,None]],axis=1)
        numpyro.factor('synthetic_optical',-.5*jnp.sum((((synthetic_y-flux(x))/sig)*maskj)**2))
    ts=time.monotonic();mcmc=MCMC(NUTS(toy_model),num_warmup=5,num_samples=1,num_chains=1,progress_bar=False)
    mcmc.run(jax.random.PRNGKey(20260926))
    result.update(warmup_steps=5,discarded_draws=1,warmup_and_draw_seconds=time.monotonic()-ts,total_seconds=time.monotonic()-start)
    save('benchmark-result.json',result)
    print(json.dumps(result),flush=True)
if __name__=='__main__':main()
