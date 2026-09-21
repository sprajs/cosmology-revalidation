#!/usr/bin/env python3
"""Independent analytic/quad challenges; no DES outcome or fit scores read."""
from pathlib import Path
import importlib.util,sys,json,math
import numpy as np
from scipy.integrate import quad
from scipy.stats import norm,t
R=Path(__file__).resolve().parents[3];O=R/'phase2/inference_audit'
manifest=json.load(open(O/'initial_source_manifest.json'))
corepath=R/next(f['snapshot'] for f in manifest['files'] if f['path'].endswith('/core.py'))
spec=importlib.util.spec_from_file_location('core',corepath);core=importlib.util.module_from_spec(spec);sys.modules['core']=core;spec.loader.exec_module(core)
import jax,jax.numpy as jnp
out={}
# Exact Gaussian-conditioning calculation from joint covariance, not root code.
linear=[]
for sx in [.1,.5,.8,1.2]:
 for rho in [-.5,0,.5]:
  latent_var=1.;slope=.3;sm=.1;sy=.05;cross=rho*sx*sy
  varx=latent_var+sx*sx;covxy=slope*latent_var+cross;true_slope=covxy/varx
  true_var=slope*slope*latent_var+sm*sm+sy*sy-covxy*covxy/varx
  projected=sy*sy+slope*slope*sx*sx-2*slope*cross+sm*sm
  mean_error_variance=(slope-true_slope)**2*varx
  kl=.5*(np.log(projected/true_var)+(true_var+mean_error_variance)/projected-1)
  linear.append({'latent_x_sd':1.,'observed_x_error_sd':sx,'measurement_rho':rho,'physical_slope':slope,'exact_observed_predictor_conditional_slope':true_slope,'exact_conditional_variance':true_var,'root_residual_propagation_variance_at_physical_slope':projected,'expected_logscore_loss_nats_per_object_at_physical_slope':kl})
out['finite_predictor_population']=linear
# At noisy predictor zero, finite Gaussian latent prior gives shrunk posterior variance.
hinges=[]
for obs in [-.05,0,.05]:
 tau=.08;err=.05;base=3.;hinge=2.;measurement_m=.03;scatter=.1
 postvar=1/(1/tau**2+1/err**2);postsd=np.sqrt(postvar);postmean=postvar*obs/err**2
 f=lambda c:base*c+hinge*max(c,0)
 em=quad(lambda u:f(postmean+postsd*u)*norm.pdf(u),-10,10,points=[-postmean/postsd],epsabs=1e-12)[0]
 e2=quad(lambda u:f(postmean+postsd*u)**2*norm.pdf(u),-10,10,points=[-postmean/postsd],epsabs=1e-12)[0]
 approxmean=f(obs);approxvar=(base+hinge*(obs>0))**2*err**2+measurement_m**2+scatter**2
 hinges.append({'observed_colour':obs,'colour_population_sigma':tau,'colour_measurement_sigma':err,'exact_latent_colour_posterior_mean':postmean,'exact_latent_colour_posterior_sigma':postsd,'exact_predictive_mean_offset':em,'plug_in_mean_offset':approxmean,'mean_error_mag':approxmean-em,'exact_predictive_variance':e2-em**2+measurement_m**2+scatter**2,'projected_variance':approxvar})
out['colour_hinge']=hinges
quadratic=[]
for obs in [0.,1.,2.]:
 tau=1.;err=.5;gamma=.1;scatter=.1;me=.05;v=1/(1/tau**2+1/err**2);m=v*obs/err**2
 quadratic.append({'x_observed':obs,'gamma':gamma,'latent_x_posterior_mean':m,'latent_x_posterior_var':v,'exact_predictive_mean_offset':gamma*(m*m+v),'plug_in_mean_offset':gamma*obs*obs,'exact_predictive_variance':gamma**2*(2*v*v+4*m*m*v)+scatter**2+me**2,'projected_variance':(2*gamma*obs)**2*err**2+scatter**2+me**2})
out['quadratic_stretch']=quadratic
# Student intrinsic scatter + independent Gaussian measurement: convolution CDF.
student=[]
for intrinsic_sd,measurement_sd in [(.1,.03),(.1,.1),(.1,.2)]:
 total=np.hypot(intrinsic_sd,measurement_sd);iscl=intrinsic_sd/np.sqrt(2);tscl=total/np.sqrt(2)
 for nsigma in [1.,2.,3.,5.]:
  bound=nsigma*total
  exact=quad(lambda u:(t.cdf((bound-measurement_sd*u)/iscl,4)-t.cdf((-bound-measurement_sd*u)/iscl,4))*norm.pdf(u),-12,12,epsabs=1e-11)[0]
  approx=t.cdf(bound/tscl,4)-t.cdf(-bound/tscl,4)
  student.append({'intrinsic_sd':intrinsic_sd,'measurement_sd':measurement_sd,'bound_total_sigma':nsigma,'exact_gaussian_plus_student_coverage':exact,'student_total_variance_coverage':approx,'difference':approx-exact})
out['intrinsic_student_measurement_convolution']=student
# Normalized host mixtures differ from fractional-host substitution even when means agree.
ph=.5;step=.3;sigma=.1;y=0.;exact=.5*norm.pdf(y,-step/2,sigma)+.5*norm.pdf(y,step/2,sigma)
out['host_mixture']={'p_high':ph,'host_step':step,'within_host_sigma':sigma,'density_at_global_mean_exact':exact,'density_fractional_host_within_sigma':norm.pdf(y,0,sigma),'mixture_variance':sigma**2+ph*(1-ph)*step**2,'fractional_host_variance':sigma**2,'assessment':'Root correctly uses a mixture, not fractional-host substitution. But host weights still require their stated conditional mass model.'}
# E= tau*u substitution makes the independent integral well-behaved near tau=0.
cases=[]
for sigma in [.03,.1,.3]:
 for tau in [1e-8,1e-6,1e-4,.001,.01,.08,.5]:
  for rsigma in [-10,-3,0,3,10]:
   r=rsigma*sigma;direction=3.2;a=direction*tau/sigma;rr=r/sigma;um=max(0.,(a*rr-1)/(a*a));f=lambda u:-u-.5*(rr-a*u)**2-np.log(sigma)-.5*np.log(2*np.pi);peak=f(um)
   if um==0:
    # Scale a boundary peak by its local exponential/Gaussian width; integrating
    # directly over a huge u range can miss a narrow boundary mass.
    B=1-a*rr;scale=1/(B+a)
    val=quad(lambda t:np.exp(-B*scale*t-.5*(a*scale*t)**2),0,np.inf,epsabs=1e-12,epsrel=1e-11,limit=300)[0]
    reference=f(0)+np.log(scale)+np.log(val)
   else:
    val=quad(lambda v:np.exp(-.5*v*v),max(-um*a,-40),np.inf,epsabs=1e-12,epsrel=1e-11,limit=300)[0]
    reference=peak-np.log(a)+np.log(val)
   actual=float(core.exp_gaussian_logpdf(jnp.array([r]),jnp.array([[sigma*sigma]]),jnp.array([direction]),tau))
   cases.append({'sigma':sigma,'tau':tau,'residual_sigma':rsigma,'reference_logpdf':reference,'root_logpdf':actual,'difference':actual-reference})
out['small_dust_scale_numerical_stability']=cases
out['interpretation']='Analytical counterexamples to physical/measurement interpretation of normalized descriptive predictors; numerical root-core challenges are separate. No real heldout outcomes read.'
(O/'independent_counterexamples.json').write_text(json.dumps(out,indent=2)+'\n')
print('Maximum finite-predictor KL',max(x['expected_logscore_loss_nats_per_object_at_physical_slope'] for x in linear));print('Hinge',hinges);print('Maximum root-core logpdf error',max(abs(x['difference']) for x in cases))
