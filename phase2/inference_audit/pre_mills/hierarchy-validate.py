"""Independent numerical-integral and analytic-limit tests before DES outcomes."""
from pathlib import Path
import hashlib,json,sys,datetime
import numpy as np
from scipy.integrate import quad
from scipy.stats import multivariate_normal,exponnorm,norm
import jax,jax.numpy as jnp
from core import exp_gaussian_logpdf,log_exgauss_cdf,distance_kinematic,distance_qbins,population_logpdf,probit_selection_logz

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'phase2/hierarchy/validation';OUT.mkdir(parents=True,exist_ok=True)
rng=np.random.default_rng(2026092101);records=[]
for dim in [1,2,3]:
    for i in range(12):
        l=rng.normal(size=(dim,dim));cov=l@l.T*.02+np.eye(dim)*.015
        direction=rng.uniform(.1,4,size=dim);r=rng.normal(size=dim)*.15;tau=rng.uniform(.02,.2)
        reference=quad(lambda ee: multivariate_normal.pdf(r,mean=direction*ee,cov=cov)*np.exp(-ee/tau)/tau,0,np.inf,epsabs=1e-11,epsrel=1e-10)[0]
        actual=float(exp_gaussian_logpdf(jnp.array(r),jnp.array(cov),jnp.array(direction),tau))
        records.append({'dimension':dim,'actual':actual,'reference':float(np.log(reference)),'log_difference':actual-np.log(reference)})
assert max(abs(r['log_difference']) for r in records)<2e-8
cdf=[]
for sigma in [.08,.2,.7]:
 for scale in [.01,.08,.4]:
  for x in [-1,-.2,0,.3,1.5]:
   a=float(log_exgauss_cdf(x,sigma,scale));b=float(exponnorm.logcdf(x,scale/sigma,scale=sigma))
   cdf.append({'sigma':sigma,'exp_mean':scale,'x':x,'difference':a-b})
assert max(abs(r['difference']) for r in cdf if np.isfinite(r['difference']))<2e-8
z=np.geomspace(.001,1.2,200);exact={};distance_checks={}
for name,q in [('deSitter',-1.),('coasting',0.),('EdS',.5)]:
 chi=z if q==-1 else np.log1p(z) if q==0 else 2*(1-1/np.sqrt(1+z))
 expected=5*np.log10((1+z)*chi*299792.458/70)+25
 v=np.array(distance_kinematic(jnp.array(z),jnp.array(z),q,0.))
 vb=np.array(distance_qbins(jnp.array(z),jnp.array(z),jnp.ones(4)*q))
 distance_checks[name]={'smooth_max_mag_error':float(np.max(abs(v-expected))),'bin_max_mag_error':float(np.max(abs(vb-expected)))}
 assert max(distance_checks[name].values())<1e-11
deriv=jax.jacrev(lambda q:distance_qbins(jnp.array([.08,.2,.5,1.]),jnp.array([.08,.2,.5,1.]),q))(jnp.zeros(4))
assert np.isfinite(deriv).all()

# Integrate the selected brightness distribution independently of the 3D code.
pars={'M':-19.3,'mx':0.,'mc':-.06,'sx':1.,'sc':.05,'sm':.1,'alpha':.15,'beta':2.2,'rb':3.2,'tau':.09}
cov=np.diag([.05**2,.2**2,.02**2])[None];mu=np.array([42.]);host=np.array([0.]);zz=np.array([.4]);limit=23.;width=.2
sigma=np.sqrt(cov[0,0,0]+pars['sm']**2+(pars['alpha']*pars['sx'])**2+(pars['beta']*pars['sc'])**2)
mean=mu[0]+pars['M']-pars['alpha']*pars['mx']+pars['beta']*pars['mc']
integral=quad(lambda m:exponnorm.pdf(m,pars['rb']*pars['tau']/sigma,loc=mean,scale=sigma)*norm.cdf((limit-m)/width),mean-4,mean+5,epsabs=1e-11)[0]
analytic=float(probit_selection_logz(jnp.array(cov),jnp.array(mu),jnp.array(host),jnp.array(zz),pars,limit,width,'dust')[0])
assert abs(analytic-np.log(integral))<1e-9
result={'seed':2026092101,'exponential_convolution_checks':records,'cdf_checks':cdf,'distance_limits':distance_checks,'q0_zero_gradient_finite':True,'selection_normalization_log_error':analytic-float(np.log(integral)),'qualification':'Independent quadrature and analytic tests of the declared effective-SALT and synthetic-selection engine; no DES model claim.'}
(OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
files=[Path(__file__),Path(__file__).with_name('core.py'),ROOT/'phase2/hierarchy/uv.lock',ROOT/'docs/phase2/plan.md']
(OUT/'manifest.json').write_text(json.dumps({'purpose':'Pre-data numerical validation','time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'seed':2026092101,'inputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files},'outputs_sha256':{'phase2/hierarchy/validation/results.json':sha(OUT/'results.json')}},indent=2)+'\n')
print(json.dumps({'max_convolution_log_error':max(abs(r['log_difference']) for r in records),'distance':distance_checks,'selection_log_error':result['selection_normalization_log_error']}))
