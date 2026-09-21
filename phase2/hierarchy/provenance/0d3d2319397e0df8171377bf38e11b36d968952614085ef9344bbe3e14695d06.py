"""Independent marginal population likelihood for effective SALT observables.

This is a stated approximation to the calibrated-flux likelihood, not a claim
that SALT colour equals physical reddening or that R_B equals extinction R_V.
Selection normalization is deliberately separate and never silently omitted.
"""
import os
os.environ.setdefault('JAX_ENABLE_X64','true')
import jax
jax.config.update('jax_enable_x64',True)
import jax.numpy as jnp
from jax.scipy.special import log_ndtr
import numpy as np

LOG2PI=np.log(2*np.pi)
C_KMS=299792.458
GL_X,GL_W=np.polynomial.legendre.leggauss(64)
GL_X=jnp.asarray(GL_X);GL_W=jnp.asarray(GL_W)

def mvn_logpdf(residual,cov):
    """Batched normalized Gaussian log density; covariance must be positive."""
    chol=jnp.linalg.cholesky(cov)
    w=jax.scipy.linalg.solve_triangular(chol,residual[...,None],lower=True)[...,0]
    return -.5*(jnp.sum(w*w,axis=-1)+2*jnp.log(jnp.diagonal(chol,axis1=-2,axis2=-1)).sum(axis=-1)+residual.shape[-1]*LOG2PI)

def exp_gaussian_logpdf(residual,cov,direction,tau):
    """Integrate E~Exp(tau), residual|E~N(direction*E,cov), E>=0 exactly."""
    invr=jnp.linalg.solve(cov,residual[...,None])[...,0]
    invu=jnp.linalg.solve(cov,direction[...,None])[...,0]
    aa=jnp.sum(direction*invu,axis=-1)
    bb=jnp.sum(direction*invr,axis=-1)-1/tau
    x=bb/jnp.sqrt(aa)
    integral=.5*LOG2PI-.5*jnp.log(aa)-jnp.log(tau)+.5*x*x+log_ndtr(x)
    return mvn_logpdf(residual,cov)+integral

def effective_population(measurement_cov,distance,host,z,params,family='gaussian'):
    """Return marginal base mean/covariance before exponential dust integration.

    host is a binary latent category; uncertain membership is mixed by caller.
    Optional drifts use g=z/(1+z), zero at z=0. Parameters are effective SALT
    population quantities; physical SED dust is tested by forward simulation.
    """
    g=z/(1+z)
    mx=params['mx']+params.get('mx_z',0.)*g+params.get('mx_h',0.)*host
    mc=params['mc']+params.get('mc_z',0.)*g+params.get('mc_h',0.)*host
    alpha=params['alpha']+params.get('alpha_z',0.)*g
    beta=params['beta']+params.get('beta_z',0.)*g
    lum=params['M']+params.get('gamma',0.)*host+params.get('lum_z',0.)*g
    mean=jnp.stack([distance+lum-alpha*mx+beta*mc,mx,mc],axis=-1)
    one=jnp.ones_like(z);zero=jnp.zeros_like(z)
    bx=jnp.stack([-alpha*one,one,zero],axis=-1)
    bc=jnp.stack([beta*one,zero,one],axis=-1)
    v=measurement_cov+params['sx']**2*bx[...,None]*bx[...,None,:]+params['sc']**2*bc[...,None]*bc[...,None,:]
    v=v.at[:,0,0].add(params['sm']**2)
    return mean,v

def population_logpdf(observed,measurement_cov,distance,host,z,params,family='gaussian'):
    mean,cov=effective_population(measurement_cov,distance,host,z,params,family)
    residual=observed-mean
    if family=='gaussian':return mvn_logpdf(residual,cov)
    if family=='dust':
        rb=params['rb']+params.get('rb_h',0.)*host+params.get('rb_z',0.)*z/(1+z)
        tau=params['tau']*jnp.exp(params.get('tau_h',0.)*host+params.get('tau_z',0.)*z/(1+z))
        direction=jnp.stack([jnp.ones_like(z)*rb,jnp.zeros_like(z),jnp.ones_like(z)],axis=-1)
        return exp_gaussian_logpdf(residual,cov,direction,tau)
    raise ValueError(family)

def marginal_xc_logpdf(observed,measurement_cov,distance,host,z,params,family='gaussian'):
    mean,cov=effective_population(measurement_cov,distance,host,z,params,family)
    residual=observed[:,1:]-mean[:,1:];v=cov[:,1:,1:]
    if family=='gaussian':return mvn_logpdf(residual,v)
    tau=params['tau']*jnp.exp(params.get('tau_h',0.)*host+params.get('tau_z',0.)*z/(1+z))
    direction=jnp.stack([jnp.zeros_like(z),jnp.ones_like(z)],axis=-1)
    return exp_gaussian_logpdf(residual,v,direction,tau)

def log_exgauss_cdf(x,sigma,exp_mean):
    """CDF of N(0,sigma)+Exp(exp_mean), stable logarithmic subtraction."""
    a=x/sigma
    first=log_ndtr(a)
    second=-x/exp_mean+.5*(sigma/exp_mean)**2+log_ndtr(a-sigma/exp_mean)
    return first+jnp.log(-jnp.expm1(jnp.minimum(second-first,-1e-14)))

def probit_selection_logz(measurement_cov,distance,host,z,params,limit,width,family='gaussian'):
    """Known synthetic selection P(S|m_obs)=Phi((limit-m_obs)/width).

    This analytic toy kernel is not asserted to be the DES survey selection.
    Integrates population+measurement scatter exactly for fixed R_B.
    """
    mean,cov=effective_population(measurement_cov,distance,host,z,params,family)
    sigma=jnp.sqrt(cov[:,0,0]+width**2)
    if family=='gaussian':return log_ndtr((limit-mean[:,0])/sigma)
    rb=params['rb']+params.get('rb_h',0.)*host+params.get('rb_z',0.)*z/(1+z)
    tau=params['tau']*jnp.exp(params.get('tau_h',0.)*host+params.get('tau_z',0.)*z/(1+z))
    return log_exgauss_cdf(limit-mean[:,0],sigma,rb*tau)

def selected_logpdf(observed,measurement_cov,distance,host,z,params,selection,family='gaussian'):
    lp=population_logpdf(observed,measurement_cov,distance,host,z,params,family)
    if selection['kind']=='complete':return lp
    if selection['kind']=='synthetic_probit':
        ls=log_ndtr((selection['limit']-observed[:,0])/selection['width'])
        lz=probit_selection_logz(measurement_cov,distance,host,z,params,selection['limit'],selection['width'],family)
        return lp+ls-lz
    raise ValueError('Unknown/unvalidated selection; cannot silently omit it: '+selection['kind'])

def distance_kinematic(z,zhel,q0,q1,H0=70.):
    """Flat FLRW q(z)=q0+q1*z/(1+z), independent of matter/DE split."""
    zz=z[:,None]*(GL_X+1)/2
    loge=(1+q0+q1)*jnp.log1p(zz)-q1*zz/(1+zz)
    chi=(z/2)*jnp.sum(GL_W*jnp.exp(-loge),axis=-1)
    return 5*jnp.log10((1+zhel)*chi*C_KMS/H0)+25

def distance_qbins(z,zhel,q,edges=(0.,.15,.35,.6,1.3),H0=70.):
    """Exact segment integrals for constant q, with a stable q=0 limit."""
    edges=jnp.asarray(edges);loge=jnp.zeros_like(z);chi=jnp.zeros_like(z)
    for i in range(len(edges)-1):
        lo=edges[i];hi=jnp.maximum(lo,jnp.minimum(z,edges[i+1]))
        width=jnp.log1p(hi)-jnp.log1p(lo)
        a=-q[i]*width
        # Taylor branch prevents singular gradient at exactly q=0.
        safe=jnp.where(jnp.abs(a)<1e-5,1.,a)
        exprel=jnp.where(jnp.abs(a)<1e-5,1+a/2+a*a/6+a*a*a/24,jnp.expm1(a)/safe)
        chi+=(1+lo)*jnp.exp(-loge)*width*exprel
        loge+=(1+q[i])*(jnp.log1p(edges[i+1])-jnp.log1p(lo))
    return 5*jnp.log10((1+zhel)*chi*C_KMS/H0)+25
