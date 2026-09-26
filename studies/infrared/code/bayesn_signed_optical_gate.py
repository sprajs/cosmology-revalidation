"""Forward/algebra gate only; this module never opens the sealed NIR payload."""
from __future__ import annotations

import os
os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:2])
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key] = '1'
os.environ.setdefault('XLA_FLAGS','--xla_cpu_multi_thread_eigen=false')

import hashlib
import inspect
import json
import sys
import textwrap
import time
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr, logsumexp
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'runs/research_2026_09_26/bayesn_signed_optical_pilot'
OLD = ROOT/'runs/research_2026_09_26/bayesn_distance_identification'
sys.path.insert(0,str(OLD/'official-code'))
import bayesn.bayesn_model as module
from bayesn import SEDmodel
import jax
import jax.numpy as jnp

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(name,obj):(OUT/name).write_text(json.dumps(obj,indent=2)+'\n')

def merged_reference(T,S):
    x=np.unique(np.r_[T[:,0],S[(S[:,0]>T[0,0])&(S[:,0]<T[-1,0]),0]])
    u,w=np.polynomial.legendre.leggauss(4)
    h=np.diff(x)/2;m=(x[1:]+x[:-1])/2;l=m[:,None]+h[:,None]*u
    r=np.interp(l,T[:,0],T[:,1]);spec=np.interp(l,S[:,0],S[:,1])
    return float(np.sum(h[:,None]*w*l*r*spec))

def normalized_d_integral():
    # Independent toy with negative synthetic optical flux; no observed values.
    K=.4*np.log(10);gray=.088;lo,hi=20.,50.;dref=35.
    fo=np.array([2.,4.,1.]);yo=np.array([1.8,3.7,-.2]);so=np.array([.5,.8,1.])
    fn=np.array([1.5,2.]);yn=np.array([1.6,1.7]);sn=np.array([.4,.5])
    def prior(d):return (ndtr((d-lo)/gray)-ndtr((d-hi)/gray))/(hi-lo)
    def lp(d,f,y,s):return norm.logpdf(y,np.exp(-K*(d-dref))*f,s).sum()
    def integrand(d,nir):return np.exp(lp(d,fo,yo,so)+(lp(d,fn,yn,sn) if nir else 0))*prior(d)
    aq=quad(lambda d:integrand(d,False),19,51,epsabs=1e-13,points=[20,30,34,35,36,40,50])[0]
    bq=quad(lambda d:integrand(d,True),19,51,epsabs=1e-13,points=[20,30,34,35,36,40,50])[0]
    grid=np.linspace(19,51,320001)
    amplitude=np.exp(-K*(grid-dref))
    pop=np.exp(norm.logpdf(yo[None,:],amplitude[:,None]*fo[None,:],so).sum(1))
    pnir=np.exp(norm.logpdf(yn[None,:],amplitude[:,None]*fn[None,:],sn).sum(1))
    pgrid=prior(grid)
    ag=np.trapezoid(pop*pgrid,grid)
    bg=np.trapezoid(pop*pnir*pgrid,grid)
    # Proper broad D prior integrates to one, including Gaussian gray tails.
    pn=quad(prior,15,55,epsabs=1e-12)[0]
    sd=np.sqrt(.088**2+.003**2);cdf=quad(lambda d:norm.pdf(d,41.4,sd),35,48)[0]
    return dict(log_joint_over_optical_quad=float(np.log(bq/aq)),log_grid_error=float(abs(np.log(bg/ag)-np.log(bq/aq))),
                broad_prior_normalization_15_55=float(pn),lcdm_normalization_35_48=float(cdf),negative_synthetic_optical=True)

def make_refined():
    source=textwrap.dedent(inspect.getsource(SEDmodel._setup_band_weights))
    assert source.count('self.spectrum_bins = 300')==1
    source=source.replace('self.spectrum_bins = 300','self.spectrum_bins = self.audit_spectrum_bins')
    scope={};exec(compile(source,'<resolution_override>','exec'),module.__dict__,scope)
    class Refined(SEDmodel):_setup_band_weights=scope['_setup_band_weights']
    return Refined,source

def main():
    tic=time.monotonic()
    protocol=json.loads((OUT/'execution-protocol.json').read_text())
    assert sha(OUT/'execution-protocol.json')==protocol['self_sha256_excluding_this_field'] if 'self_sha256_excluding_this_field' in protocol else True
    adapter=json.loads((OUT/'adapter-manifest.json').read_text())
    assert sha(OUT/'optical-payload.json')==adapter['outputs']['optical-payload.json']
    # Deliberately no access to nir-sealed.json, including hash, from this process.
    payload=json.loads((OUT/'optical-payload.json').read_text())
    objs=payload['objects'];rows=payload['rows'];n=len(objs);nobs=max(sum(r['cid']==o['cid'] for r in rows) for o in objs)
    assert n==2 and nobs==15
    t=np.zeros((nobs,n));errors=np.ones((nobs,n));mask=np.zeros((nobs,n));bands=np.full((nobs,n),'NULL_BAND',object)
    for j,o in enumerate(objs):
        selected=[r for r in rows if r['cid']==o['cid']]
        for i,r in enumerate(selected):
            t[i,j]=r['trigger_rest_time']; errors[i,j]=r['error'];bands[i,j]='RAISIN_DES_'+r['band'];mask[i,j]=1
    assert np.all((t[mask>0]>10)&(t[mask>0]<30))
    assert np.all((t[mask>0]-20>-10)&(t[mask>0]+10<40))
    pars=np.zeros((n,47));pars[:,0]=[o['mu_LCDM'] for o in objs];pars[:,1]=.3;pars[:,2]=3.1;pars[:,46]=5.
    Refined,override=make_refined();(OUT/'resolution-method.py').write_text(override)
    z=jnp.array([o['zHEL'] for o in objs]);ebv=jnp.array([o['MWEBV'] for o in objs]);maskj=jnp.array(mask);tj=jnp.array(t)
    reference=[];resolution=[];prev=None;last=None
    T=np.loadtxt(OUT/'RAISIN_DES_J.dat');S=np.loadtxt(OLD/'release-filters/DES_AB_primary.dat')
    exact_j=2.5*np.log10(merged_reference(T,S)/merged_reference(T,np.column_stack([T[:,0],np.ones(len(T))])))
    for bins in [300,600,1200,2400]:
        Refined.audit_spectrum_bins=bins
        start=time.monotonic();model=Refined(load_model='M20_model',num_devices=1,filter_yaml=str(OUT/'release-filter-config.yaml'))
        weights=model._calculate_band_weights(z,ebv)
        bi=jnp.array([[model.band_dict[b] for b in row] for row in bands])
        def forward(p):
            phase=tj-p[:,46][None,:]
            jt=model.J_t_map(phase.flatten(order='F'),model.tau_knots,model.KD_t).reshape((nobs,n,6),order='F').transpose(1,2,0)
            hs=jnp.array([19+jnp.floor(phase),19+jnp.ceil(phase),jnp.remainder(phase,1)])
            eps=(model.L_Sigma@p[:,4:46].T).T.reshape((n,7,6),order='F')
            eps=jnp.zeros((n,9,6)).at[:,1:-1,:].set(eps)
            return model.get_flux_batch(model.M0,p[:,3],p[:,1],model.W0,model.W1,eps,p[:,0],p[:,2],bi,maskj,jt,hs,weights)
        compiled=jax.jit(forward)
        pred=np.asarray(compiled(jnp.array(pars)));eager=np.asarray(forward(jnp.array(pars)))
        assert pred.shape==(nobs,n) and np.isfinite(pred).all()
        gap=float(np.max(abs(pred-eager)/errors*mask));delta=None if prev is None else float(np.max(abs(pred-prev)/errors*mask))
        resolution.append(dict(bins=bins,eager_jit_error_sigma=gap,last_delta_sigma=delta,seconds=time.monotonic()-start))
        reference.append(dict(bins=bins,exact_J_mag=float(exact_j),model_J_mag=float(model.zp_dict['RAISIN_DES_J']),
                              delta_mag=float(model.zp_dict['RAISIN_DES_J']-exact_j)))
        np.savez_compressed(OUT/f'forward-{bins}.npz',prediction=pred,errors=errors,mask=mask,trigger_rest_time=t,pars=pars,
                            CID=np.array([o['cid'] for o in objs]))
        print(json.dumps(dict(resolution=resolution[-1],reference=reference[-1])),flush=True)
        prev=pred;last=(model,forward,compiled,pred)
        if time.monotonic()-tic>900:raise TimeoutError('Forward/algebra 15-minute cap')
    model,forward,compiled,pred=last
    p=jnp.array(pars)
    # Endpoint and phase-knot crossing checks are independent of observed flux.
    phase_checks=[]
    for tau in [-10.,20.]:
        q=pars.copy();q[:,46]=tau
        a=np.asarray(compiled(jnp.array(q)));b=np.asarray(forward(jnp.array(q)))
        phase_checks.append(dict(tau=tau,max_eager_jit_sigma=float(np.max(abs(a-b)/errors*mask)),finite=bool(np.isfinite(a).all())))
    knots=np.asarray(model.tau_knots)
    crossing=[]
    for j,o in enumerate(objs):
        i=0;ph=t[i,j];candidate=ph-knots
        valid=candidate[(candidate>-10)&(candidate<20)]
        for tau in [float(valid[0]),float(valid[len(valid)//2]),float(valid[-1])]:
            q=pars.copy();q[j,46]=tau
            at=np.asarray(compiled(jnp.array(q)));q[j,46]=tau-1e-6;left=np.asarray(compiled(jnp.array(q)))
            q[j,46]=tau+1e-6;right=np.asarray(compiled(jnp.array(q)))
            crossing.append(dict(cid=o['cid'],tau=tau,finite=bool(np.isfinite(at).all()),
                                 symmetric_jump_sigma=float(np.max(abs(left+right-2*at)/errors*mask))))
    # All 47 parameter finite derivatives, with half-step consistency; tau is dynamic.
    steps=np.r_[1e-4,1e-4,1e-3,1e-4,np.repeat(1e-4,42),1e-4]
    def differences(scale):
        cols=[]
        for k,h in enumerate(steps*scale):
            plus=pars.copy();minus=pars.copy();plus[:,k]+=h;minus[:,k]-=h
            cols.append((np.asarray(compiled(jnp.array(plus)))-np.asarray(compiled(jnp.array(minus))))/(2*h))
        return np.stack(cols,-1)
    J=differences(1);Jhalf=differences(.5)
    relative=float(np.linalg.norm((J-Jhalf)/errors[:,:,None]*mask[:,:,None])/max(1.,np.linalg.norm(Jhalf/errors[:,:,None]*mask[:,:,None])))
    q=pars.copy();q[:,1]=0;q[:,2]=3.2;up=np.asarray(compiled(jnp.array(q)));q[:,2]=3.;down=np.asarray(compiled(jnp.array(q)))
    av0=float(np.max(abs(up-down)/errors*mask))
    q=pars.copy();q[:,0]+=.15;gray=np.asarray(compiled(jnp.array(q)));expected=pred*10**(-.4*.15)
    gray_error=float(np.max(abs(gray-expected)/np.maximum(abs(expected),1e-200)*mask))
    # Synthetic linear functional gradient checks all coordinates, no measured y.
    coeff=jnp.array(np.where(mask>0,1/errors,0.))
    def scalar(x):return jnp.sum(compiled(x)*coeff)
    grad=np.asarray(jax.grad(scalar)(p));fd=np.sum(Jhalf*np.asarray(coeff)[:,:,None],axis=0)
    graderr=float(np.linalg.norm(grad-fd)/max(1.,np.linalg.norm(fd)))
    algebra=normalized_d_integral()
    result=dict(scope='No observed likelihood, posterior, or NIR flux access',adapter_manifest_sha256=sha(OUT/'adapter-manifest.json'),
      optical_payload_sha256=sha(OUT/'optical-payload.json'),resolution=resolution,reference=reference,
      phase_endpoints=phase_checks,knot_crossings=crossing,derivative_halfstep_relative=relative,
      gradient_scaled_error=graderr,AV0_RV_error_sigma=av0,gray_scaling_relative=gray_error,
      D_algebra=algebra,elapsed_seconds=time.monotonic()-tic)
    result['gates_pass']=bool(all(r['eager_jit_error_sigma']<1e-8 for r in resolution) and resolution[-1]['last_delta_sigma']<.1 and
      abs(reference[-1]['delta_mag'])<.001 and relative<1e-4 and graderr<1e-5 and gray_error<1e-12 and av0<1e-12 and
      algebra['log_grid_error']<1e-6 and all(x['finite'] for x in phase_checks+crossing))
    save('forward-algebra-result.json',result)
    np.savez_compressed(OUT/'derivative-gate.npz',J=J,Jhalf=Jhalf,parameters=pars,prediction=pred,errors=errors,mask=mask)
    print(json.dumps(dict(gates_pass=result['gates_pass'],elapsed_seconds=result['elapsed_seconds'],halfstep=relative,gradient=graderr)),flush=True)

if __name__=='__main__':main()
