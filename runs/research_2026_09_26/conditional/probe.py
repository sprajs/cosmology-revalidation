import sys,time
from pathlib import Path
import numpy as np
import jax,jax.numpy as jnp
from numpyro.infer import MCMC,NUTS,init_to_median
from numpyro.infer.util import initialize_model
R=Path(__file__).resolve().parents[3];sys.path.insert(0,str(R/'scripts/phase2/hierarchy'))
from conditional import model
from core import distance_kinematic
a=dict(np.load(R/'phase2/hierarchy/data/conditioned-multistart-best/data.npz',allow_pickle=False));sel=(a['survey']==10)&(a['pIa']>.999)&(a['fold']!=0)
d={k:jnp.asarray(v[sel]) for k,v in a.items()};d['distance_reference']=distance_kinematic(d['z'],d['zhel'],0.,0.)
key=jax.random.PRNGKey(2026092160)
info=initialize_model(key,model,init_strategy=init_to_median(num_samples=20),model_args=(d,'tripp','gaussian'))
p=info.param_info.z;f=info.potential_fn;g=jax.grad(f)(p)
for k in p:
 x=np.asarray(p[k]); idx=(0,)*x.ndim
 for h in (1e-4,1e-6):
  delta=np.zeros_like(x);delta[idx]=h
  fp=float(f({**p,k:jnp.asarray(x+delta)}));fm=float(f({**p,k:jnp.asarray(x-delta)})); fd=(fp-fm)/(2*h)
  print('GRAD',k,h,float(np.asarray(g[k])[idx]),fd,flush=True)
for size in (0.001,0.00001):
 t=time.time();m=MCMC(NUTS(model,step_size=size,adapt_step_size=False,adapt_mass_matrix=False,dense_mass=False,max_tree_depth=5,init_strategy=init_to_median(num_samples=20)),num_warmup=1,num_samples=10,num_chains=1,progress_bar=False)
 m.run(key,d,'tripp','gaussian',extra_fields=('diverging','num_steps','accept_prob'))
 e=m.get_extra_fields();s=m.get_samples()
 print('FIXED',size,'time',time.time()-t,'accept',np.asarray(e['accept_prob']),'steps',np.asarray(e['num_steps']),'div',np.asarray(e['diverging']),'M',np.asarray(s['M']),flush=True)
for size in (0.001,1e-6):
 n=NUTS(model,step_size=size,adapt_step_size=False,adapt_mass_matrix=False,dense_mass=False,max_tree_depth=3,init_strategy=init_to_median(num_samples=20))
 state=n.init(key,num_warmup=1,model_args=(d,'tripp','gaussian'))
 print('STATE_INIT',size,'pe',float(state.potential_energy),'energy',float(state.energy),'zM',float(state.z['M']),'inv_mass',state.adapt_state.inverse_mass_matrix,flush=True)
 st2=n.sample(state,(d,'tripp','gaussian'),{})
 print('STATE_NEXT',size,'pe',float(st2.potential_energy),'energy',float(st2.energy),'accept',float(st2.accept_prob),'steps',int(st2.num_steps),'zM',float(st2.z['M']),'div',bool(st2.diverging),flush=True)
n=NUTS(model,step_size=.001,adapt_step_size=False,adapt_mass_matrix=False,dense_mass=False,max_tree_depth=5,init_strategy=init_to_median(num_samples=20))
st=n.init(key,num_warmup=1,model_args=(d,'tripp','gaussian'))
for i in range(5):
 st=n.sample(st,(d,'tripp','gaussian'),{})
 print('TRANSITION',i+1,'pe',float(st.potential_energy),'accept',float(st.accept_prob),'M',float(st.z['M']),'step',float(st.adapt_state.step_size),'inv_mass',np.asarray(st.adapt_state.inverse_mass_matrix[('M','coefficients','offsets','scatter')])[:3],flush=True)
n=NUTS(model,step_size=.001,adapt_step_size=False,adapt_mass_matrix=False,dense_mass=False,max_tree_depth=5,init_strategy=init_to_median(num_samples=20))
st=n.init(key,num_warmup=1,model_args=(d,'tripp','gaussian'))
advance=jax.jit(lambda ss:n.sample(ss,(d,'tripp','gaussian'),{}))
for i in range(5):
 st=advance(st)
 print('JIT_TRANSITION',i+1,'pe',float(st.potential_energy),'accept',float(st.accept_prob),'M',float(st.z['M']),'steps',int(st.num_steps),flush=True)
info=initialize_model(key,model,init_strategy=init_to_median(num_samples=20),model_args=(d,'tripp','gaussian'))
p=info.param_info.z;f=info.potential_fn
jf=jax.jit(f);jg=jax.jit(jax.grad(f));ge=jax.grad(f)(p);gj=jg(p)
print('COMPARE_P',float(f(p)),float(jf(p)),flush=True)
for k in p:print('COMPARE_GRAD',k,float(jnp.max(jnp.abs(ge[k]-gj[k]))),np.asarray(ge[k]).reshape(-1)[:3],np.asarray(gj[k]).reshape(-1)[:3],flush=True)
params=info.postprocess_fn(p)
from conditional import predictive,components
pe=predictive(params,d,'tripp','gaussian')[0];pj=jax.jit(lambda pp:predictive(pp,d,'tripp','gaussian')[0])(params)
print('PRED_COMPARE',float(pe.sum()),float(pj.sum()),float(jnp.max(jnp.abs(pe-pj))),flush=True)
for k in params:print('POST',k,np.asarray(params[k]),flush=True)
print('REPEAT',float(f(p)),float(f(p)),float(jf(p)),float(jf(p)),flush=True)
ce,ve=components(params,d,'tripp');cj,vj=jax.jit(lambda pp:components(pp,d,'tripp'))(params)
for label,e,j in [('mu',ce,cj),('variance',ve,vj)]:
 ee=np.asarray(e);jj=np.asarray(j);idx=np.argmax(np.abs(ee-jj)); print('COMPONENT',label,'max',float(np.max(np.abs(ee-jj))),'row',idx//2,'e',ee.reshape(-1)[idx],'j',jj.reshape(-1)[idx],flush=True)
from conditional import KNOTS
base=lambda pp: d['distance_reference']+pp['M']+jnp.interp(d['z'],KNOTS,jnp.concatenate([jnp.zeros(1),pp['offsets']]))
be=base(params);bj=jax.jit(base)(params)
print('BASE_DIFF',float(jnp.max(jnp.abs(be-bj))),flush=True)
from conditional import features,FEATURES,MODELS
idx=jnp.array([FEATURES.index(x) for x in MODELS['tripp']]);_,dx,dc=features(d,0);v=jnp.stack([jnp.ones(len(d['z'])),-dx[:,idx]@params['coefficients'],-dc[:,idx]@params['coefficients']],axis=-1)
ve2=jnp.einsum('ni,nij,nj->n',v,d['cov'],v);vj2=jax.jit(lambda vv:jnp.einsum('ni,nij,nj->n',vv,d['cov'],vv))(v)
row=322
print('V',np.asarray(v[row]),'COV',np.asarray(d['cov'][row]),'PROJ_E',float(ve2[row]),'PROJ_J',float(vj2[row]),'MANUAL',float(v[row]@d['cov'][row]@v[row]),flush=True)
print('COVFINITE',bool(jnp.all(jnp.isfinite(d['cov']))),'vfinite',bool(jnp.all(jnp.isfinite(v))),flush=True)
def get_v(pp):
 _,dx,dc=features(d,0); ii=jnp.asarray([FEATURES.index(x) for x in MODELS['tripp']],dtype=int)
 return jnp.stack([jnp.ones(len(d['z'])),-dx[:,ii]@pp['coefficients'],-dc[:,ii]@pp['coefficients']],axis=-1)
vv_e=get_v(params);vv_j=jax.jit(get_v)(params)
print('V_DIFF',float(jnp.max(jnp.abs(vv_e-vv_j))),'V_ROW',np.asarray(vv_e[161]),np.asarray(vv_j[161]),flush=True)
for h in (0,1):
 ff,ddx,ddc=features(d,h)
 for label,arr in [('f',ff),('dx',ddx),('dc',ddc)]:
  fn=lambda aa:aa[:,jnp.asarray([0,1],dtype=int)]
  e=fn(arr);j=jax.jit(fn)(arr)
  print('GATHER',h,label,'max',float(jnp.max(jnp.abs(e-j))),'e',np.asarray(e[161]),'j',np.asarray(j[161]),flush=True)
for label,arr in [('dx',ddx),('dc',ddc),('f',ff)]:
 fn=lambda aa,cc:aa[:,jnp.asarray([0,1],dtype=int)]@cc
 e=fn(arr,params['coefficients']);j=jax.jit(fn)(arr,params['coefficients'])
 print('MATMUL',label,'max',float(jnp.max(jnp.abs(e-j))),'e',float(e[161]),'j',float(j[161]),flush=True)
def get_parts(pp):
 _,ax,ac=features(d,0); ii=jnp.asarray([0,1],dtype=int)
 return ax[:,ii]@pp['coefficients'],ac[:,ii]@pp['coefficients']
for label,e,j in zip(('gx','gc'),get_parts(params),jax.jit(get_parts)(params)):
 print('PART',label,'max',float(jnp.max(jnp.abs(e-j))),'e',float(e[161]),'j',float(j[161]),flush=True)
def get_v2(pp):
 gx,gc=get_parts(pp)
 return jnp.stack([jnp.ones(len(d['z'])),-gx,-gc],axis=-1)
qe=get_v2(params);qj=jax.jit(get_v2)(params)
print('V2_DIFF',float(jnp.max(jnp.abs(qe-qj))),np.asarray(qe[161]),np.asarray(qj[161]),flush=True)
def get_v3(pp):
 _,ax,ac=features(d,0); ii=jnp.asarray([0,1],dtype=int)
 gx=ax[:,ii]@pp['coefficients']; gc=ac[:,ii]@pp['coefficients']
 return jnp.stack([jnp.ones(len(d['z'])),-gx,-gc],axis=-1)
qe=get_v3(params);qj=jax.jit(get_v3)(params)
print('V3_DIFF',float(jnp.max(jnp.abs(qe-qj))),np.asarray(qe[161]),np.asarray(qj[161]),flush=True)
