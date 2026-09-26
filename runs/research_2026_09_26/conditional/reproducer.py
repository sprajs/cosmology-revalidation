"""Run with pinned Python/JAX to compare covariance projection implementations."""
import ast,json,os,platform
import numpy as np
import jax,jax.numpy as jnp,jaxlib
jax.config.update('jax_enable_x64',True)

# Keep the two selected feature derivatives and their matrix product visible.
# The padded columns mimic the conditional model's eight-feature stack.
x=np.array([1.0,1.0,1.0]); z=np.zeros_like(x)
dx=np.stack([x,z,z,z,z,z,z,z],axis=-1)
dc=np.stack([z,x,z,z,z,z,z,z],axis=-1)
coef=np.array([-0.20623263862481683,3.4957056719902764])
indices=np.array([0,1]); cov=np.array([[.0033094,.04752451,.00219889],[.04752451,1.57899487,.01996853],[.00219889,.01996853,.00270431]])

# Preserve exactly the original expression as the source of the investigation.
def fused(c):
    one=jnp.ones(3)
    ax=jnp.asarray(dx);ac=jnp.asarray(dc);ii=jnp.asarray(indices)
    return jnp.stack([one,-ax[:,ii]@c,-ac[:,ii]@c],axis=-1)
def split(c):
    one=jnp.ones(3)
    ax=jnp.asarray(dx);ac=jnp.asarray(dc);ii=jnp.asarray(indices)
    gx=ax[:,ii]@c;gc=ac[:,ii]@c
    return jnp.stack([one,-gx,-gc],axis=-1)

expected=np.stack([np.ones(3),-(dx[:,indices]@coef),-(dc[:,indices]@coef)],axis=-1)
result={'environment':{'python':platform.python_version(),'jax':jax.__version__,'jaxlib':jaxlib.__version__,'backend':jax.default_backend(),'x64':bool(jax.config.jax_enable_x64),'XLA_FLAGS':os.environ.get('XLA_FLAGS'),'OPENBLAS_NUM_THREADS':os.environ.get('OPENBLAS_NUM_THREADS'),'OMP_NUM_THREADS':os.environ.get('OMP_NUM_THREADS')},'ast':ast.dump(ast.parse('-dx[:,indices]@coef',mode='eval'),include_attributes=False),'expected_v':expected[0].tolist()}
for name,fn in [('fused_eager',fused),('fused_jit',jax.jit(fused)),('split_jit',jax.jit(split))]:
    v=np.asarray(fn(jnp.asarray(coef)))[0]
    variance=float(v@cov@v)
    grad=np.asarray(jax.grad(lambda c:jnp.einsum('i,ij,j->',fn(c)[0],jnp.asarray(cov),fn(c)[0]))(jnp.asarray(coef)))
    result[name]={'v':v.tolist(),'variance':variance,'variance_grad':grad.tolist()}
# Independent NumPy gradient for v=[1,-a,-b], a=coefficient[0], b=coefficient[1].
v=expected[0];dv=np.array([[0.,-1.,0.],[0.,0.,-1.]])
result['numpy']={'variance':float(v@cov@v),'variance_grad':(2*dv@cov@v).tolist()}
print(json.dumps(result,indent=2))
