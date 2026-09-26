#!/usr/bin/env python3
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.special import erfcx,ndtr
from scipy.stats import norm
from scipy.integrate import quad
R=Path(__file__).resolve().parents[3];O=R/'phase2/inference_audit';sys.path.insert(0,str(R/'scripts/phase2/hierarchy'))
import core,jax,jax.numpy as jnp
tail=[]
for x in [-1e8,-1e5,-1000.,-100.,-30.,-20.000001,-20.,-19.999999,-10.,0.,10.]:
 actual=float(core.scaled_log_normal_cdf(x));reference=float(np.log(erfcx(-x/np.sqrt(2.)))-np.log(2.)) if x<=0 else .5*x*x+norm.logcdf(x)
 grad=float(jax.grad(core.scaled_log_normal_cdf)(x))
 if x < -19.:
  t=-x;den=quad(lambda v:np.exp(-v-.5*(v/t)**2),0,np.inf,epsabs=1e-12)[0];num=quad(lambda v:v*np.exp(-v-.5*(v/t)**2),0,np.inf,epsabs=1e-12)[0];gref=num/(t*den)
 else:gref=x+np.exp(norm.logpdf(x)-norm.logcdf(x))
 tail.append({'x':x,'actual':actual,'reference_erfcx':reference,'difference':actual-reference,'gradient':grad,'gradient_reference_integral':gref,'gradient_difference':grad-gref})
 assert abs(actual-reference)<2e-12 and np.isfinite(grad) and abs(grad-gref)<2e-11
old=json.load(open(O/'independent_counterexamples.json'))['small_dust_scale_numerical_stability'];checks=[]
for row in old:
 s=row['sigma'];tau=row['tau'];r=s*row['residual_sigma'];actual=float(core.exp_gaussian_logpdf(jnp.array([r]),jnp.array([[s*s]]),jnp.array([3.2]),tau));checks.append({**row,'updated_logpdf':actual,'updated_difference':actual-row['reference_logpdf']})
assert max(abs(c['updated_difference']) for c in checks)<2e-10
limits=[]
for tau in [1e-6,1e-8,1e-10,1e-12]:
 for sigma in [.03,.1,.3]:
  for k in [-10,-3,0,3,10]:
   x=k*sigma;actual=float(core.exp_gaussian_logpdf(jnp.array([x]),jnp.array([[sigma*sigma]]),jnp.array([3.2]),tau));gaussian=norm.logpdf(x,scale=sigma)
   cdf=float(core.log_exgauss_cdf(x,sigma,3.2*tau))
   # Independent integral over unit exponential u; rescale by Gaussian CDF
   # to retain relative precision in small lower-tail probabilities.
   base=norm.logcdf(k);value=quad(lambda u:np.exp(norm.logcdf((x-3.2*tau*u)/sigma)-base-u),0,50,epsabs=1e-12,epsrel=1e-11)[0];ref=base+np.log(value)
   limits.append({'tau':tau,'sigma':sigma,'x_sigma':k,'logpdf_minus_gaussian':actual-gaussian,'cdf_log_error':cdf-ref})
assert max(abs(c['cdf_log_error']) for c in limits)<2e-10
assert max(abs(c['logpdf_minus_gaussian']) for c in limits if c['tau']==1e-12)<2e-8
result={'tail_function':tail,'original_grid_recheck':checks,'gaussian_limits_and_cdf':limits,'max_updated_grid_error':max(abs(c['updated_difference']) for c in checks),'max_tail_error':max(abs(c['difference']) for c in tail),'max_cdf_error':max(abs(c['cdf_log_error']) for c in limits),'core_sha256':hashlib.file_digest((R/'scripts/phase2/hierarchy/core.py').open('rb'),'sha256').hexdigest(),'passed':True,'scientific_scope':'Numerical robustness at tiny dust scales; no empirical population inference or correction ranking.'}
(O/'stable_tail_validation.json').write_text(json.dumps(result,indent=2)+'\n');print({k:v for k,v in result.items() if not isinstance(v,list)})
