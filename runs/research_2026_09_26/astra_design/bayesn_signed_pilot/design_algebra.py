"""Pure model/prior and synthetic probability algebra; no observed flux values."""
from pathlib import Path
import ast,json
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
from scipy.special import ndtr,logsumexp
from scipy.stats import norm
P=Path(__file__).resolve().parent;S=P.parent/'next_experiment_review/bayesn_dust_support.py'
tree=ast.parse(S.read_text());ns={'np':np,'CubicSpline':CubicSpline,'xk':np.array([0,1e4/26500,1e4/12200,1e4/6000,1e4/5470,1e4/4670,1e4/4110,1e4/2700,1e4/2600])}
exec(compile(ast.Module(body=[x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name in ('knots','min_law')],type_ignores=[]),str(S),'exec'),ns)
rvs=np.linspace(1.2,6,4097);mins=np.array([ns['min_law'](x)[0] for x in rvs]);assert mins.min()>0
# Synthetic optical posterior on gray amplitude with proper prior; compare ratio of joint and optical integrals to conditional posterior mixture.
s=.088;lo,hi=20,50;Dref=35.;K=.4*np.log(10);fO=np.array([2.,4.,1.]);yO=np.array([1.8,3.7,-.2]);sigO=np.array([.5,.8,1.]);fN=np.array([1.5,2.]);yN=np.array([1.6,1.7]);sigN=np.array([.4,.5])
def prior(d):return (ndtr((d-lo)/s)-ndtr((d-hi)/s))/(hi-lo)
def lO(d):return norm.pdf(yO,np.exp(-K*(d-Dref))*fO,sigO).prod()
def lN(d):return norm.pdf(yN,np.exp(-K*(d-Dref))*fN,sigN).prod()
Z=quad(lambda d:lO(d)*prior(d),19,51,epsabs=1e-13,points=[30,34,35,36,40,50])[0]
ZN=quad(lambda d:lO(d)*lN(d)*prior(d),19,51,epsabs=1e-13,points=[30,34,35,36,40,50])[0]
d=np.linspace(19,51,320001);a=np.exp(-K*(d-Dref));L1=norm.pdf(yO[None,:],a[:,None]*fO[None,:],sigO).prod(1);L2=norm.pdf(yN[None,:],a[:,None]*fN[None,:],sigN).prod(1);p=(ndtr((d-lo)/s)-ndtr((d-hi)/s))/(hi-lo);pred=np.trapezoid(L1*L2*p,d)/np.trapezoid(L1*p,d)
err=abs(np.log(pred)-np.log(ZN/Z));assert err<1e-9
r={'status':'Pure source-law and toy probability algebra, no real observed likelihoods','F99':{'RV_domain':[1.2,6.],'RV_grid_N':4097,'wavelength_domain_A':[3000,18500],'each_RV_min':'all cubic stationary points and wavelength endpoints','minimum_A_over_AV':float(mins.min()),'minimum_RV_grid':float(rvs[np.argmin(mins)]),'caveat':'Dense RV verification, not empirical validation of the law or a universal physical RV bound.'},'proper_prior_predictive':{'joint_over_optical_logdensity':float(np.log(ZN/Z)),'independent_grid_vs_quad_log_error':err,'noise_includes_negative_optical_flux':True},'amplitude':{'D':'mu+gray, grayNormal(0,.088²)','broad':'muUniform20..50, convolved Normal; integrate normalized density in D','prediction':'p(NIR|opt)=integral p(NIR|latent,D)p(latent,D|opt); preserve entire common latent draw across all NIR rows'}}
(P/'design-algebra-result.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
