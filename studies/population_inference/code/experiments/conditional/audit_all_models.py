"""Audit all conditional feature maps against independent NumPy algebra."""
import sys,json,hashlib,os
from pathlib import Path
import numpy as np
import jax,jax.numpy as jnp
from scipy.special import logsumexp
R=Path(__file__).resolve().parents[3];sys.path.insert(0,str(R/'scripts/phase2/hierarchy'))
from conditional import MODELS,FEATURES,PRIORS,KNOTS,components,predictive,features
from core import distance_kinematic
raw=dict(np.load(R/'phase2/hierarchy/data/conditioned-multistart-best/data.npz',allow_pickle=False));sel=(raw['survey']==10)&(raw['pIa']>.999)&(raw['fold']!=0)
d={k:jnp.asarray(v[sel]) for k,v in raw.items()};d['distance_reference']=distance_kinematic(d['z'],d['zhel'],0.,0.)
n={k:np.asarray(v) for k,v in d.items()}; cohort_sha=hashlib.sha256((R/'phase2/hierarchy/conditional-cohort.json').read_bytes()).hexdigest()
def legacy_components(params,data,name):
    indices=jnp.asarray([FEATURES.index(x) for x in MODELS[name]],dtype=int)
    baseline=data['distance_reference']+params['M']+jnp.interp(data['z'],KNOTS,jnp.concatenate([jnp.zeros(1),params['offsets']]))
    coef=params.get('coefficients',jnp.empty(0));means=[];variances=[]
    for h in [0,1]:
        f,dx,dc=features(data,h);mu=baseline+f[:,indices]@coef
        v=jnp.stack([jnp.ones(len(mu)),-dx[:,indices]@coef,-dc[:,indices]@coef],axis=-1)
        variance=jnp.einsum('ni,nij,nj->n',v,data['cov'],v)+params['scatter'][h]**2
        means.append(mu);variances.append(variance)
    return jnp.stack(means,axis=-1),jnp.stack(variances,axis=-1)
def np_components(params,name):
    x=n['y'][:,1];c=n['y'][:,2];g=n['z']/(1+n['z'])-.2
    baseline=n['distance_reference']+params['M']+np.interp(n['z'],np.asarray(KNOTS),np.r_[0.,params['offsets']]);coef=params.get('coefficients',np.empty(0));idx=[FEATURES.index(k) for k in MODELS[name]]
    means=[];variances=[];gradvars=[]
    for h in [0.,1.]:
        hc=h-.5;one=np.ones_like(x);zero=np.zeros_like(x)
        f=np.stack([x,c,hc*one,hc*c,np.maximum(c,0),x*g,c*g,x*x],axis=-1)[:,idx]
        dx=np.stack([one,zero,zero,zero,zero,g,zero,2*x],axis=-1)[:,idx]
        dc=np.stack([zero,one,zero,hc*one,(c>0).astype(float),zero,g,zero],axis=-1)[:,idx]
        mu=baseline+f@coef;v=np.stack([one,-dx@coef,-dc@coef],axis=-1)
        variance=np.einsum('ni,nij,nj->n',v,n['cov'],v)+params['scatter'][int(h)]**2
        means.append(mu);variances.append(variance)
        dv=np.stack([np.zeros_like(dx),-dx,-dc],axis=-1)
        gradvars.append(2*np.einsum('nki,nij,nj->k',dv,n['cov'],v))
    return np.stack(means,axis=-1),np.stack(variances,axis=-1),np.stack(gradvars,axis=0)
def json_float(v):return float(np.asarray(v))
results=[]
for name,features_names in MODELS.items():
    prior=np.asarray([PRIORS[k] for k in features_names]).reshape(-1,2)
    for draw in range(3):
        coef=prior[:,0]+(.25*draw-.25)*prior[:,1]
        params={'M':-19.3,'offsets':np.linspace(-.05,.1,7),'scatter':np.array([.14,.18])}
        if len(coef):params['coefficients']=coef
        jp={k:jnp.asarray(v) for k,v in params.items()}
        em,ev,ng=np_components(params,name)
        oldm,oldv=jax.jit(lambda pp:legacy_components(pp,d,name))(jp)
        newm,newv=jax.jit(lambda pp:components(pp,d,name))(jp)
        # Whole predictive function is compiled; compare its Gaussian mixture to NumPy.
        lp,mean,var=jax.jit(lambda pp:predictive(pp,d,name,'gaussian'))(jp)
        ph=np.clip(n['host_prob'],1e-10,1-1e-10)
        lps=-.5*(np.log(2*np.pi*ev)+(n['y'][:,0,None]-em)**2/ev)
        nlp=logsumexp(np.stack([np.log1p(-ph)+lps[:,0],np.log(ph)+lps[:,1]],axis=-1),axis=-1)
        nmean=(1-ph)*em[:,0]+ph*em[:,1]
        nvar=(1-ph)*ev[:,0]+ph*ev[:,1]+ph*(1-ph)*(em[:,1]-em[:,0])**2
        # Analytic NumPy gradient of sum of both host variances.
        if len(coef):
            grad=np.asarray(jax.jit(jax.grad(lambda cc:components({**jp,'coefficients':cc},d,name)[1].sum()))(jnp.asarray(coef)))
            grad_err=float(np.max(np.abs(grad-ng.sum(axis=0))))
        else:grad_err=0.
        results.append({'model':name,'draw':draw,'legacy_max_mean_error':float(np.max(np.abs(np.asarray(oldm)-em))),'legacy_max_variance_error':float(np.max(np.abs(np.asarray(oldv)-ev))),'fixed_max_mean_error':float(np.max(np.abs(np.asarray(newm)-em))),'fixed_max_variance_error':float(np.max(np.abs(np.asarray(newv)-ev))),'fixed_max_logpdf_error':float(np.max(np.abs(np.asarray(lp)-nlp))),'fixed_max_predictive_mean_error':float(np.max(np.abs(np.asarray(mean)-nmean))),'fixed_max_predictive_variance_error':float(np.max(np.abs(np.asarray(var)-nvar))),'fixed_max_variance_gradient_error':grad_err})
        print(name,draw,'legacy variance max',results[-1]['legacy_max_variance_error'],'fixed variance max',results[-1]['fixed_max_variance_error'],flush=True)
assert len(results)==27
for row in results:
    for key,value in row.items():
        if key.startswith('fixed_'):assert value<1e-8,(row['model'],row['draw'],key,value)
output={'n_train':int(sel.sum()),'cohort_sha256':cohort_sha,'checks':results}
(R/'runs/research_2026_09_26/conditional/all_models_audit.json').write_text(json.dumps(output,indent=2)+'\n')
