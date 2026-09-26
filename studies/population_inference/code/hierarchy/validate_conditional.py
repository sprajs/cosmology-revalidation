"""Independent density integration and finite-difference covariance checks."""
import hashlib,json
from pathlib import Path
import numpy as np
import jax,jax.numpy as jnp
from scipy.integrate import quad
from conditional import MODELS,PRIORS,predictive,features

R=Path(__file__).resolve().parents[3];out=R/'phase2/hierarchy/conditional-validation';out.mkdir(parents=True,exist_ok=True)
rng=np.random.default_rng(2026092159);err=np.array([.07,.3,.04]);cov=np.array([[1,.25,-.18],[.25,1,.22],[-.18,.22,1.]])*err[:,None]*err[None,:]
d={'y':jnp.array([[23.4,.7,.09]]),'cov':jnp.asarray(cov[None]),'host_prob':jnp.array([.37]),'z':jnp.array([.5]),'distance_reference':jnp.array([42.5])};results=[]
for name,names in MODELS.items():
    coef=np.array([PRIORS[k][0]+.25*PRIORS[k][1] for k in names]);p={'M':jnp.array(-19.3),'offsets':jnp.linspace(0.,.2,7),'scatter':jnp.array([.11,.16])}
    if len(coef):p['coefficients']=jnp.asarray(coef)
    for noise in ['gaussian','student4']:
        f=jax.jit(lambda m:predictive(p,{**d,'y':d['y'].at[0,0].set(m)},name,noise)[0][0])
        # Integrate the whole normalized conditional density independently.
        integral,error=quad(lambda v:np.exp(float(f(v))),-np.inf,np.inf,epsabs=2e-9,points=None,limit=500)
        # Infinite-range quadrature can miss a narrow density near magnitude23;
        # shift the integration coordinate to the separately returned mean.
        centre=float(predictive(p,d,name,noise)[1][0])
        integral,error=quad(lambda v:np.exp(float(f(v+centre))),-np.inf,np.inf,epsabs=2e-10,limit=500)
        assert abs(integral-1)<2e-8,(name,noise,integral)
        results.append({'model':name,'noise':noise,'density_integral':integral,'quadrature_error_estimate':error})
    # Independent central differences of conditional mean at known host state;
    # compare Gaussian measurement propagation with Monte Carlo perturbations.
    dh={**d,'host_prob':jnp.array([0.])};grad=[]
    for col in [1,2]:
        ep=1e-5
        plus={**dh,'y':dh['y'].at[0,col].add(ep)};minus={**dh,'y':dh['y'].at[0,col].add(-ep)}
        grad.append(float((predictive(p,plus,name,'gaussian')[1][0]-predictive(p,minus,name,'gaussian')[1][0])/(2*ep)))
    vec=np.array([1.,-grad[0],-grad[1]]);expected=vec@cov@vec
    total=float(predictive(p,dh,name,'gaussian')[2][0]);assert abs(total-.11**2-expected)<1e-9,(name,total,expected)
    draws=rng.multivariate_normal(np.zeros(3),cov,size=50000);observed=np.var(draws@vec,ddof=1)
    assert abs(observed/expected-1)<.04
    results.append({'model':name,'projected_measurement_variance':expected,'monte_carlo_variance':observed,'relative_difference':observed/expected-1})
(out/'results.json').write_text(json.dumps({'checks':results,'qualification':'Numerical tests of a declared empirical conditional predictor. This does not establish astrophysical correctness or selection-free inference.'},indent=2)+'\n')
files=[Path(__file__),Path(__file__).with_name('conditional.py'),Path(__file__).with_name('core.py')];sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
(out/'manifest.json').write_text(json.dumps({'inputs_sha256':{str(f.relative_to(R)):sha(f) for f in files},'outputs_sha256':{str((out/'results.json').relative_to(R)):sha(out/'results.json')}},indent=2)+'\n');print('All18density integrals and9covariance MonteCarlo checks passed')
