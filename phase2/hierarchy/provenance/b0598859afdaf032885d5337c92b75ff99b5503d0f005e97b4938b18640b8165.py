"""Repeated synthetic recovery and prespecified misspecification challenges.

Frequentist MLE/Hessian coverage is a diagnostic, not Bayesian SBC. Generation
uses NumPy draws and independent SciPy distance quadrature. No DES data enter.
"""
import os
os.environ.setdefault('JAX_PLATFORMS','cpu')
os.environ.setdefault('XLA_FLAGS','--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1')
import argparse,datetime,hashlib,json,time
from pathlib import Path
import numpy as np
import jax,jax.numpy as jnp
from scipy.integrate import quad
from scipy.optimize import minimize
from scipy.special import ndtr
from scipy.stats import binomtest
from core import population_logpdf,probit_selection_logz,distance_kinematic

ROOT=Path(__file__).resolve().parents[3]
NAMES=['M','alpha','beta','mx','mc','sx','sc','sm','mx_h','rb','tau','q0','q1']
TRUTH=np.array([-19.3,.15,2.2,0.,-.06,1.,.05,.1,-.35,3.4,.085,-.5,1.2])
SCALE=np.array([1.,.1,1.,1.,.1,1.,1.,1.,1.,1.,1.,1.,1.])
LOGIDX=[5,6,7,10]

def encode(p):
    u=np.array(p,copy=True);u[LOGIDX]=np.log(u[LOGIDX]);return u/SCALE

def decode(u):
    v=u*jnp.asarray(SCALE[:len(u)])
    v=v.at[jnp.array(LOGIDX)].set(jnp.exp(v[jnp.array(LOGIDX)]))
    return {k:v[i] for i,k in enumerate(NAMES[:len(u)])}

def generate(n,seed,stress):
    rng=np.random.default_rng(seed);z=rng.uniform(.02,.95,n);h=rng.integers(0,2,n)
    # Deliberately independent of the inference distance implementation.
    inv_e=lambda v:np.exp(-(1+TRUTH[11]+TRUTH[12])*np.log1p(v)+TRUTH[12]*v/(1+v))
    mu=np.array([5*np.log10((1+v)*299792.458/70*quad(inv_e,0,v,epsabs=1e-11)[0])+25 for v in z])
    err=np.column_stack([.035+.05*z/.95,.12+.25*z/.95,.012+.025*z/.95]);rho=np.array([[1,.15,-.12],[.15,1,.05],[-.12,.05,1.]])
    cov=err[:,:,None]*err[:,None,:]*rho;chol=np.linalg.cholesky(cov);y=np.zeros((n,3));pending=np.ones(n,bool);proposals=0
    while pending.any():
        ids=np.flatnonzero(pending);nn=len(ids);proposals+=nn
        x=rng.normal(TRUTH[3]+TRUTH[8]*h[ids],TRUTH[5]);ci=rng.normal(TRUTH[4],TRUTH[6],nn)
        e=rng.exponential(TRUTH[10],nn)
        if stress=='dust_mixture':e*=np.where(rng.random(nn)<.25,3.,1./3.) # Same mean E, non-exponential shape.
        scatter=rng.normal(0.,TRUTH[7],nn)
        if stress=='heavy_tails':scatter=rng.standard_t(4,nn)*TRUTH[7]/np.sqrt(2.) # Same variance.
        drift=.30*z[ids]/(1+z[ids]) if stress=='luminosity_drift' else 0.
        m=mu[ids]+TRUTH[0]-TRUTH[1]*x+TRUTH[2]*ci+TRUTH[9]*e+scatter+drift
        cand=np.column_stack([m,x,ci+e])+np.einsum('nij,nj->ni',chol[ids],rng.normal(size=(nn,3)))
        accept=rng.random(nn)<ndtr((23.4-cand[:,0])/.3)
        y[ids[accept]]=cand[accept];pending[ids[accept]]=False
        if proposals>n*100000:raise RuntimeError('Unexpectedly low selection acceptance')
    return dict(y=y,cov=cov,z=z,zhel=z,host_prob=h.astype(float),distance=mu,limit=np.ones(n)*23.4,width=np.ones(n)*.3),proposals

def objective(u,data,free_q,selected):
    p=decode(u);mu=distance_kinematic(data['z'],data['z'],p['q0'],p['q1']) if free_q else data['distance']
    lp=population_logpdf(data['y'],data['cov'],mu,data['host_prob'],data['z'],p,'dust')
    if selected:lp-=probit_selection_logz(data['cov'],mu,data['host_prob'],data['z'],p,data['limit'],data['width'],'dust')
    # log epsilon(y) is constant in parameters and omitted for optimization only.
    return -lp.sum()

def main():
    a=argparse.ArgumentParser();a.add_argument('--replicates',type=int,default=30);a.add_argument('--n',type=int,default=1000);a.add_argument('--seed',type=int,default=2026092120);args=a.parse_args()
    out=ROOT/'phase2/hierarchy/coverage';out.mkdir(parents=True,exist_ok=True)
    config={'arguments':vars(args),'truth':dict(zip(NAMES,TRUTH)),'design':'Repeated independent selected draws; MLE plus observed-Hessian errors, not posterior SBC. Correct-law fixed-distance recovery first, then free-q tests with omitted selection, same-variance t4 scatter, equal-mean two-exponential dust, and 0.30*z/(1+z) luminosity drift. At most four seeds per stress; descriptive falsification, not precision coverage.','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
    (out/'configuration.json').write_text(json.dumps(config,indent=2)+'\n')
    code={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [Path(__file__),Path(__file__).with_name('core.py')]}
    for f in [Path(__file__),Path(__file__).with_name('core.py')]:
        dst=ROOT/'phase2/hierarchy/provenance'/(code[str(f.relative_to(ROOT))]+'.py');dst.write_bytes(f.read_bytes())
    bounds=[(-21,-17),(.2,3.5),(.1,6),(-3,3),(-3,3),(np.log(.2),np.log(3)),(np.log(.008),np.log(.3)),(np.log(.015),np.log(.5)),(-3,3),(1,6),(np.log(.008),np.log(.4)),(-2,1.5),(-6,8)]
    results=[];began=time.time()
    # Compile once per distance/selection mode, with data arrays passed as arguments.
    engines={}
    for fq,sel in [(False,True),(True,True),(True,False)]:
        fn=lambda u,d,fq=fq,sel=sel:objective(u,d,fq,sel)
        engines[(fq,sel)]=(jax.jit(jax.value_and_grad(fn)),jax.jit(jax.hessian(fn)))
    cases=[('correct',False,True,args.replicates),('correct',True,True,4),('correct',True,False,4),('heavy_tails',True,True,4),('dust_mixture',True,True,4),('luminosity_drift',True,True,4)]
    for stress,free_q,selected,count in cases:
        for rep in range(min(count,args.replicates)):
            seed=args.seed+rep;data,proposals=generate(args.n,seed,stress);datafile=out/f'data-{stress}-{seed}.npz'
            if not datafile.exists():np.savez_compressed(datafile,**data)
            d={k:jnp.asarray(v) for k,v in data.items()};vg,hes=engines[(free_q,selected)];dim=13 if free_q else 11
            def fg(u):v,g=vg(u,d);return float(v),np.asarray(g,dtype=float)
            fit=minimize(fg,encode(TRUTH)[:dim],jac=True,method='L-BFGS-B',bounds=bounds[:dim],options={'maxiter':1200,'ftol':1e-12,'gtol':1e-5,'maxls':40})
            # An alternate start avoids presenting a single initialized basin as unique.
            alt=encode(TRUTH)[:dim].copy();alt[2]+=.7;alt[9]-=.4;alt[10]+=.4
            second=minimize(fg,alt,jac=True,method='L-BFGS-B',bounds=bounds[:dim],options={'maxiter':1200,'ftol':1e-12,'gtol':1e-5,'maxls':40})
            gap=abs(fit.fun-second.fun)
            if second.fun<fit.fun:fit=second
            hessian=np.asarray(hes(fit.x,d));eig=np.linalg.eigvalsh(hessian);phys=decode(jnp.asarray(fit.x));estimate=np.array([float(phys[k]) for k in NAMES[:dim]])
            jac=SCALE[:dim].copy();jac[LOGIDX]=estimate[LOGIDX]
            sd=jac*np.sqrt(np.maximum(0,np.diag(np.linalg.inv(hessian)))) if eig.min()>0 else np.full(dim,np.nan)
            boundary=[NAMES[k] for k in range(dim) if min(fit.x[k]-bounds[k][0],bounds[k][1]-fit.x[k])<1e-4]
            record={'stress':stress,'free_q':free_q,'selection_normalized':selected,'seed':seed,'proposals':proposals,'success':bool(fit.success),'message':str(fit.message),'objective':float(fit.fun),'alternate_start_objective_gap':float(gap),'max_abs_gradient':float(np.abs(fit.jac).max()),'minimum_hessian_eigenvalue':float(eig.min()),'boundary_parameters':boundary,'parameters':{k:{'estimate':float(estimate[i]),'standard_error':float(sd[i]),'truth':float(TRUTH[i]),'z_score':float((estimate[i]-TRUTH[i])/sd[i])} for i,k in enumerate(NAMES[:dim])}}
            results.append(record);(out/'replicates.json').write_text(json.dumps(results,indent=2)+'\n')
            print(json.dumps({'case':stress,'q':free_q,'sel':selected,'seed':seed,'ok':record['success'],'gradient':record['max_abs_gradient'],'gap':gap,'tau_z':record['parameters']['tau']['z_score'],'q0':record['parameters'].get('q0')}),flush=True)
    good=[r for r in results if r['stress']=='correct' and not r['free_q'] and r['success'] and r['minimum_hessian_eigenvalue']>0 and not r['boundary_parameters']]
    coverage={}
    for key in NAMES[:11]:
        zz=np.array([r['parameters'][key]['z_score'] for r in good]);n95=int((np.abs(zz)<1.96).sum());ci=binomtest(n95,len(zz)).proportion_ci()
        coverage[key]={'mean_z':float(zz.mean()),'sd_z':float(zz.std(ddof=1)),'covered_95':n95,'valid_replicates':len(zz),'coverage_binomial_95_interval':[ci.low,ci.high]}
    (out/'summary.json').write_text(json.dumps({'coverage':coverage,'runtime_seconds':time.time()-began,'interpretation':'Monte Carlo uncertainty is large for 30 seeds. MLE-Hessian checks are neither SBC nor proof of empirical population correctness. Misspecification cases are explicit synthetic counterexamples.'},indent=2)+'\n')
    (out/'manifest.json').write_text(json.dumps({'code_sha256':code,'outputs_sha256':{str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in out.iterdir() if f.name!='manifest.json'}},indent=2)+'\n')
if __name__=='__main__':main()
