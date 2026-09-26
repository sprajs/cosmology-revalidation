"""NumPyro implementation of the independent, marginalized population engine."""
import jax
import jax.numpy as jnp
from jax.scipy.special import log_ndtr,logsumexp
import numpyro
import numpyro.distributions as dist
from core import population_logpdf,probit_selection_logz,distance_kinematic,distance_qbins

MODELS={
 'G0':{'family':'gaussian','population_drift':False,'coefficient_drift':False,'host_dust':False,'luminosity_host':True},
 'G1':{'family':'gaussian','population_drift':True,'coefficient_drift':False,'host_dust':False,'luminosity_host':True},
 'G2':{'family':'gaussian','population_drift':True,'coefficient_drift':True,'host_dust':False,'luminosity_host':True},
 'D0':{'family':'dust','population_drift':False,'coefficient_drift':False,'host_dust':False,'luminosity_host':False},
 'Dhost':{'family':'dust','population_drift':True,'coefficient_drift':False,'host_dust':True,'luminosity_host':False},
 'DhostM':{'family':'dust','population_drift':True,'coefficient_drift':False,'host_dust':True,'luminosity_host':True},
 'Ddrift':{'family':'dust','population_drift':True,'coefficient_drift':False,'host_dust':True,'luminosity_host':True,'dust_drift':True},
}

def parameters(spec):
    p={}
    def draw(k,d):p[k]=numpyro.sample(k,d)
    draw('M',dist.Normal(-19.3,.6));draw('alpha',dist.Uniform(.02,.35));draw('beta',dist.Uniform(.1,6.))
    draw('mx',dist.Normal(0.,1.));draw('mc',dist.Normal(-.04,.15))
    draw('sx',dist.LogNormal(0.,.5));draw('sc',dist.LogNormal(jnp.log(.07),.7));draw('sm',dist.HalfNormal(.2))
    draw('mx_h',dist.Normal(0.,.6))
    if spec['family']=='gaussian':draw('mc_h',dist.Normal(0.,.08))
    if spec['luminosity_host']:draw('gamma',dist.Normal(0.,.15))
    if spec['population_drift']:
        draw('mx_z',dist.Normal(0.,2.));
        if spec['family']=='gaussian':draw('mc_z',dist.Normal(0.,.25))
    if spec['coefficient_drift']:
        draw('alpha_z',dist.Normal(0.,.15));draw('beta_z',dist.Normal(0.,1.5))
    if spec['family']=='dust':
        draw('rb',dist.Uniform(1.,6.));draw('tau',dist.LogNormal(jnp.log(.08),.7))
        if spec['host_dust']:
            draw('tau_h',dist.Normal(0.,.7));draw('rb_high',dist.Uniform(1.,6.));p['rb_h']=p['rb_high']-p['rb']
        if spec.get('dust_drift'):draw('tau_z',dist.Normal(0.,1.))
    return p

def distance(data,mode):
    z=data['z'];zh=data['zhel']
    if mode=='fixed':return data['distance']
    if mode=='kinematic':
        q0=numpyro.sample('q0',dist.Uniform(-2.,1.5));q1=numpyro.sample('q1',dist.Normal(1.,2.))
        return distance_kinematic(z,zh,q0,q1)
    if mode=='qbins':
        q=numpyro.sample('q',dist.Uniform(-2.,2.).expand([4]).to_event(1))
        numpyro.factor('q_smoothing',dist.Normal(0.,data.get('q_smoothing',1.)).log_prob(jnp.diff(q)).sum())
        return distance_qbins(z,zh,q)
    if mode=='offsetbins':
        knots=jnp.array([.01,.08,.2,.35,.5,.7,.9,1.2])
        offset=numpyro.sample('distance_offsets',dist.Normal(0.,.5).expand([len(knots)-1]).to_event(1))
        return distance_kinematic(z,zh,0.,0.)+jnp.interp(z,knots,jnp.concatenate([jnp.zeros(1),offset]))
    raise ValueError(mode)

def pointwise_logpdf(params,mu,data,spec):
    n=len(data['z']);host0=jnp.zeros(n);host1=jnp.ones(n)
    args=(data['y'],data['cov'],mu)
    a=population_logpdf(*args,host0,data['z'],params,spec['family'])
    b=population_logpdf(*args,host1,data['z'],params,spec['family'])
    ph=jnp.clip(data['host_prob'],1e-8,1-1e-8)
    ll=jnp.logaddexp(jnp.log1p(-ph)+a,jnp.log(ph)+b)
    if data['selection']=='complete':return ll
    if data['selection']=='synthetic_probit':
        za=probit_selection_logz(data['cov'],mu,host0,data['z'],params,data['limit'],data['width'],spec['family'])
        zb=probit_selection_logz(data['cov'],mu,host1,data['z'],params,data['limit'],data['width'],spec['family'])
        z=jnp.logaddexp(jnp.log1p(-ph)+za,jnp.log(ph)+zb)
        return ll+log_ndtr((data['limit']-data['y'][:,0])/data['width'])-z
    raise ValueError('Real survey selection is not yet implemented/validated; refusing '+data['selection'])

def model(data,name='D0',distance_mode='fixed'):
    spec=MODELS[name];p=parameters(spec);mu=distance(data,distance_mode)
    ll=pointwise_logpdf(p,mu,data,spec)
    numpyro.factor('observation_loglikelihood',ll.sum())

