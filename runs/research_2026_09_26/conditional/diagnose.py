import sys,json,time
from pathlib import Path
import numpy as np,pandas as pd
import jax,jax.numpy as jnp
from numpyro.infer import MCMC,NUTS,init_to_median
from numpyro.infer.util import initialize_model
R=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(R/'scripts/phase2/hierarchy'))
from conditional import model,components,ROOT
from core import distance_kinematic
folder=R/'phase2/hierarchy/data/conditioned-multistart-best'
a=dict(np.load(folder/'data.npz',allow_pickle=False));sel=(a['survey']==10)&(a['pIa']>.999)&(a['fold']!=0)
d={k:jnp.asarray(v[sel]) for k,v in a.items()};d['distance_reference']=distance_kinematic(d['z'],d['zhel'],0.,0.)
keys=jax.random.split(jax.random.PRNGKey(2026092160),2)
res=[]
for key in keys:
 info=initialize_model(key,model,init_strategy=init_to_median(num_samples=20),model_args=(d,'tripp','gaussian'))
 pos=info.param_info.z;potential=info.potential_fn
 grad=jax.grad(potential)(pos)
 res.append({'potential':float(potential(pos)),'gradient_l2':float(jnp.sqrt(sum(jnp.sum(v*v) for v in grad.values()))),'gradient_max':{k:float(jnp.max(jnp.abs(v))) for k,v in grad.items()},'position':{k:np.asarray(v).tolist() for k,v in pos.items()},'finite_gradient':bool(all(jnp.all(jnp.isfinite(v)) for v in grad.values()))})
 print('INIT',res[-1],flush=True)
start=time.time()
m=MCMC(NUTS(model,target_accept_prob=.9,max_tree_depth=10,dense_mass=True,init_strategy=init_to_median(num_samples=20)),num_warmup=30,num_samples=20,num_chains=2,chain_method='sequential',progress_bar=True)
m.run(jax.random.PRNGKey(2026092160),d,'tripp','gaussian',extra_fields=('diverging','num_steps','accept_prob'))
s=m.get_samples(group_by_chain=True);ex=m.get_extra_fields(group_by_chain=True)
state=m.last_state
out={'initial':res,'elapsed':time.time()-start,'accept':np.asarray(ex['accept_prob']).tolist(),'steps':np.asarray(ex['num_steps']).tolist(),'divergences':np.asarray(ex['diverging']).tolist(),'chain_std':{k:np.asarray(v).std(axis=1).tolist() for k,v in s.items()},'adapt_step_size':np.asarray(state.adapt_state.step_size).tolist(),'adapt_inverse_mass_matrix':np.asarray(state.adapt_state.inverse_mass_matrix).tolist()}
(R/'runs/research_2026_09_26/conditional/diagnostic_30w20d.json').write_text(json.dumps(out,indent=2)+'\n')
print('RESULT',json.dumps({k:out[k] for k in ('elapsed','chain_std','adapt_step_size')}),flush=True)
