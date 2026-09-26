"""Reproduce eager/compiled projection discrepancy on frozen training inputs."""
from pathlib import Path
import json,os,platform
import numpy as np
import jax,jax.numpy as jnp,jaxlib
jax.config.update('jax_enable_x64',True)
R=Path(__file__).resolve().parents[3]
a=dict(np.load(R/'phase2/hierarchy/data/conditioned-multistart-best/data.npz',allow_pickle=False))
mask=(a['survey']==10)&(a['pIa']>.999)&(a['fold']!=0)
y=jnp.asarray(a['y'][mask]); cov=jnp.asarray(a['cov'][mask]);z=jnp.asarray(a['z'][mask]);h=0
x=y[:,1];c=y[:,2];hc=h-.5;g=z/(1+z)-.2;one=jnp.ones_like(x);zero=jnp.zeros_like(x)
# These arrays are copied verbatim in structure from conditional.features.
dx=jnp.stack([one,zero,zero,zero,zero,g,zero,2*x],axis=-1)
dc=jnp.stack([zero,one,zero,hc*one,(c>0).astype(float),zero,g,zero],axis=-1)
idx=jnp.asarray([0,1],dtype=int)
coef=jnp.asarray([-.20623263862481683,3.4957056719902764])
def fused(cc):return jnp.stack([jnp.ones(len(x)),-dx[:,idx]@cc,-dc[:,idx]@cc],axis=-1)
def split(cc):
    gx=dx[:,idx]@cc;gc=dc[:,idx]@cc
    return jnp.stack([jnp.ones(len(x)),-gx,-gc],axis=-1)
def projection(fn,cc):
    v=fn(cc)
    return jnp.einsum('ni,nij,nj->n',v,cov,v)
expected_v=np.stack([np.ones(len(x)),-np.asarray(dx)[:,[0,1]]@np.asarray(coef),-np.asarray(dc)[:,[0,1]]@np.asarray(coef)],axis=-1)
expected_var=np.einsum('ni,nij,nj->n',expected_v,np.asarray(cov),expected_v)
res={'environment':{'python':platform.python_version(),'jax':jax.__version__,'jaxlib':jaxlib.__version__,'backend':jax.default_backend(),'x64':bool(jax.config.jax_enable_x64),'XLA_FLAGS':os.environ.get('XLA_FLAGS'),'OPENBLAS_NUM_THREADS':os.environ.get('OPENBLAS_NUM_THREADS'),'OMP_NUM_THREADS':os.environ.get('OMP_NUM_THREADS')},'n_train':len(x),'numpy':{'v_row161':expected_v[161].tolist(),'variance_row161':float(expected_var[161]),'variance_sum':float(expected_var.sum())}}
for label,fn in [('fused_eager',fused),('fused_jit',jax.jit(fused)),('split_jit',jax.jit(split))]:
    vv=np.asarray(fn(coef)); pp=np.asarray(projection(fn,coef)); grad=np.asarray(jax.grad(lambda cc:projection(fn,cc).sum())(coef))
    res[label]={'v_row161':vv[161].tolist(),'variance_row161':float(pp[161]),'variance_sum':float(pp.sum()),'max_abs_v_error':float(np.max(np.abs(vv-expected_v))),'max_abs_variance_error':float(np.max(np.abs(pp-expected_var))),'grad_variance_sum':grad.tolist()}
print(json.dumps(res,indent=2))
