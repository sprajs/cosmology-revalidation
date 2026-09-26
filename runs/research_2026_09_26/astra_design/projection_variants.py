"""Independent expression variants for the pinned-JAX projection failure."""
from pathlib import Path
import json, os, platform, sys, hashlib
import numpy as np
import jax, jax.numpy as jnp, jaxlib
jax.config.update('jax_enable_x64', True)
ROOT = Path(__file__).resolve().parents[3]
a = dict(np.load(ROOT/'phase2/hierarchy/data/conditioned-multistart-best/data.npz'))
mask=(a['survey']==10)&(a['pIa']>.999)&(a['fold']!=0)
x=a['y'][mask,1];c=a['y'][mask,2];z=a['z'][mask]
o=np.ones_like(x);zz=np.zeros_like(x);g=z/(1+z)-.2
dx=np.column_stack([o,zz,zz,zz,zz,g,zz,2*x])
dc=np.column_stack([zz,o,zz,-.5*o,(c>0).astype(float),zz,g,zz])
jdx,jdc=jnp.array(dx),jnp.array(dc)
ids=jnp.array([0,1]);coef=np.array([-.20623263862481683,3.4957056719902764])
expected=np.column_stack([o,-dx[:,[0,1]]@coef,-dc[:,[0,1]]@coef])
variants={
 'original_closure': lambda b:jnp.stack([jnp.ones(len(x)),-jdx[:,ids]@b,-jdc[:,ids]@b],axis=-1),
 'parenthesized_closure': lambda b:jnp.stack([jnp.ones(len(x)),-(jdx[:,ids]@b),-(jdc[:,ids]@b)],axis=-1),
 'multiply_negative_closure': lambda b:jnp.stack([jnp.ones(len(x)),-1*(jdx[:,ids]@b),-1*(jdc[:,ids]@b)],axis=-1),
 'no_negation_closure': lambda b:jnp.stack([jnp.ones(len(x)),jdx[:,ids]@b,jdc[:,ids]@b],axis=-1)*jnp.array([1.,-1.,-1.]),
 'matmul_preselected': lambda b:jnp.stack([jnp.ones(len(x)), -(jnp.array(dx[:,[0,1]])@b), -(jnp.array(dc[:,[0,1]])@b)],axis=-1),
 'explicit_scalar_tripp': lambda b:jnp.broadcast_to(jnp.array([1.,-b[0],-b[1]]),(len(x),3)),
}
result={'environment':{'python':platform.python_version(),'jax':jax.__version__,'jaxlib':jaxlib.__version__,'XLA_FLAGS':os.environ.get('XLA_FLAGS')},'variants':{}}
cov=a['cov'][mask]
expected_grad=np.array([-88.54629360885976,13.197675581788143])
for name,f in variants.items():
    jf=jax.jit(f);got=np.asarray(jf(jnp.array(coef)))
    grad=np.asarray(jax.grad(lambda b:jnp.einsum('ni,nij,nj->',jf(b),cov,jf(b)))(jnp.array(coef)))
    result['variants'][name]={'max_abs_error':float(abs(got-expected).max()),'row161':got[161].tolist(),'variance_sum':float(np.einsum('ni,nij,nj->',got,cov,got)),'gradient':grad.tolist(),'max_gradient_error':float(abs(grad-expected_grad).max())}
    if name in ('original_closure','parenthesized_closure'):
        (Path(__file__).parent/f'{name}.stablehlo.txt').write_text(str(jf.lower(jnp.array(coef)).compiler_ir('stablehlo')))
def dynamic(b,xx,cc):
    return jnp.stack([jnp.ones(xx.shape[0]),-xx[:,ids]@b,-cc[:,ids]@b],axis=-1)
got=np.asarray(jax.jit(dynamic)(jnp.array(coef),jdx,jdc))
result['variants']['original_dynamic_arguments']={'max_abs_error':float(abs(got-expected).max()),'row161':got[161].tolist()}
print(json.dumps(result,indent=2))
