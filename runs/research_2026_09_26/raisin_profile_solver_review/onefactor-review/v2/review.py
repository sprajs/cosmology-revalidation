"""Prespecified adaptive-quadrature comparison; no observed means or fitting."""
from pathlib import Path
import hashlib, importlib.util, itertools, json, warnings
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import log_ndtr, erfcx
from scipy.stats import multivariate_normal

ROOT=Path(__file__).resolve().parents[5];OUT=Path(__file__).resolve().parent
SOURCE=ROOT/'runs/research_2026_09_26/onefactor_sign_validation/v2/onefactor_sign_likelihood.py'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
save=lambda name,o:(OUT/name).write_text(json.dumps(o,indent=2,allow_nan=False)+'\n')
spec=importlib.util.spec_from_file_location('kernel_v2',SOURCE);k=importlib.util.module_from_spec(spec);spec.loader.exec_module(k)
protocol=dict(scope='Synthetic means and covariance-only checks, no observed outcomes.',kernel_sha256=sha(SOURCE),
    synthetic_design=dict(dimensions=[2,4,16,74],factor_precisions=[.01,.25,.5],diagonal_variance='logspace0.5to2',SNR=[-3.,0.,3.],signs=['allpositive','alternating']),
    near_mode_boundary=dict(dimensions=[2,16],factor_precision=[.25,.5],target_modes=[-4.01,-3.99,3.99,4.01]),
    native_design='Native C_B only; standardized mean -10,-3,0,3,10; allpositive,alternating,first9negative patterns.',
    reference='Independent stable erfcx inverse-Mills derivative; bracket exact latent mode; mode-scaled scipy adaptive integration on mode +/-12; strong log concavity bounds omitted tails.',
    acceptance_tolerance_log_probability=2e-10,orders=[64,128,256],tail_guard_stresses=[-1e6,-1e12,-1e150],
    analytic_checks='Strong scalar; empty; independent; explicit sparse covariance refusal; synthetic conditional posterior latent against adaptive integration.')
save('protocol.json',protocol)

def mills(t):
    t=np.asarray(t);out=np.empty_like(t);negative=t<0
    out[negative]=np.sqrt(2/np.pi)/erfcx(-t[negative]/np.sqrt(2))
    out[~negative]=np.exp(-.5*t[~negative]**2-.5*np.log(2*np.pi)-log_ndtr(t[~negative]))
    return out
def reference(mu,d,v,s):
    a=s*mu/np.sqrt(d);b=s*v/np.sqrt(d)
    score=lambda z:float(-z+np.sum(b*mills(a+b*z)))
    radius=4.
    while score(-radius)<0 or score(radius)>0:radius*=2
    mode=brentq(score,-radius,radius,xtol=1e-13)
    ell=lambda z:float(-.5*z*z-.5*np.log(2*np.pi)+np.sum(log_ndtr(a+b*z)))
    maximum=ell(mode)
    integral,error=quad(lambda u:np.exp(ell(mode+u)-maximum),-12,12,epsabs=2e-13,epsrel=2e-13,limit=150,points=[0])
    return maximum+np.log(integral),mode,error/integral

cases=[]
for n,p,snr,pattern in itertools.product([2,4,16,74],[.01,.25,.5],[-3.,0.,3.],['allpositive','alternating']):
    d=np.geomspace(.5,2.,n);b=np.full(n,np.sqrt(p/n));v=b*np.sqrt(d);mu=snr*np.sqrt(d);s=np.ones(n) if pattern=='allpositive' else np.where(np.arange(n)%2,1.,-1.)
    cases.append((f'design_n{n}_p{p}_snr{snr}_{pattern}',mu,d,v,s))
for n,p,target in itertools.product([2,16],[.25,.5],[-4.01,-3.99,3.99,4.01]):
    d=np.geomspace(.5,2.,n);b=np.full(n,np.sign(target)*np.sqrt(p/n));v=b*np.sqrt(d)
    a=brentq(lambda a:float(np.sum(b*mills(a+b*target))-target),-100,20)
    cases.append((f'boundary_n{n}_p{p}_mode{target}',a*np.sqrt(d),d,v,np.ones(n)))
nativepath=ROOT/'runs/research_2026_09_26/astra_design/raisin_signed_refit/fixed-c-profile/global-profile/native-profiles.npz'
with np.load(nativepath) as z:C=z['covariance_B'].copy()
nd,nv=k.decompose_covariance(C)
for snr,pattern in itertools.product([-10.,-3.,0.,3.,10.],['allpositive','alternating','first9negative']):
    s={'allpositive':np.ones(74),'alternating':np.where(np.arange(74)%2,1.,-1.),'first9negative':np.r_[-np.ones(9),np.ones(65)]}[pattern]
    cases.append((f'native_snr{snr}_{pattern}',snr*np.sqrt(np.diag(C)),nd,nv,s))

results=[]
for name,mu,d,v,s in cases:
    exact,mode,quaderr=reference(mu,d,v,s);precision=float(np.sum(v*v/d));record=dict(name=name,factor_precision=precision,mode=mode,adaptive_log_probability=float(exact),adaptive_relative_error_estimate=float(quaderr))
    try:
        values=[k.log_orthant(mu,d,v,s,order) for order in [64,128,256]]
        record.update(status='accepted',logprob=values,max_abs_error=float(np.max(abs(np.array(values)-exact))))
        assert abs(mode)<=4+1e-10 and precision<=.5
        assert record['max_abs_error']<2e-10,(name,record)
    except ValueError as exc:
        record.update(status='refused',reason=str(exc))
        assert abs(mode)>4-1e-10 or precision>.5,(name,record)
    results.append(record)

conditional=[]
for snr in [0.,2.,5.]:
    mu=snr*np.sqrt(np.diag(C));obs=np.arange(74)<65;y=mu[obs]+np.sqrt(np.diag(C)[obs])*np.linspace(-1,1,65)
    # Entirely synthetic values; generic censor density does not require observed signs.
    V=1/(1+np.sum(nv[obs]**2/nd[obs]));M=V*np.sum(nv[obs]*(y-mu[obs])/nd[obs])
    logp,mode,qe=reference(mu[~obs]+nv[~obs]*M,nd[~obs],nv[~obs]*np.sqrt(V),-np.ones(9))
    exact=multivariate_normal.logpdf(y,mean=mu[obs],cov=C[np.ix_(obs,obs)])+logp
    values=[k.log_censored(y,mu,nd,nv,obs,-np.ones(9),o) for o in [64,128,256]]
    err=float(np.max(abs(np.array(values)-exact)));assert err<2e-10
    conditional.append(dict(SNR=snr,adaptive_logdensity=float(exact),values=values,max_abs_error=err))

stress=[]
for a in [-1e6,-1e12,-1e150]:
    mu=np.full(2,a);d=np.ones(2);v=np.full(2,.1);score4=float(-4+np.sum(v*mills(mu+4*v)))
    with warnings.catch_warnings(record=True) as ww:
        warnings.simplefilter('always')
        try: val=k.log_orthant(mu,d,v,np.ones(2));result=dict(status='accepted',value=float(val) if np.isfinite(val) else str(val))
        except ValueError as exc:result=dict(status='refused',reason=str(exc))
    stress.append(dict(mean=a,stable_score_at4=score4,mathematical_mode_outside=True,**result,warnings=[str(w.message) for w in ww]))
analytic=[]
for loading in [.3,3.,1000.]:
    exact=float(log_ndtr(-.2/np.sqrt(1+loading**2)));value=k.log_orthant([-.2],[1.],[loading],[1.]);assert value==exact
    analytic.append(dict(loading=loading,exact=exact,value=value))
result=dict(kernel_sha256=sha(SOURCE),prespecified_cases=len(results),accepted=sum(r['status']=='accepted' for r in results),refused=sum(r['status']=='refused' for r in results),
    max_accepted_adaptive_error=max(r.get('max_abs_error',0) for r in results),max_adaptive_reported_relative_error=max(r['adaptive_relative_error_estimate'] for r in results),cases=results,conditional_checks=conditional,
    strong_scalar_exact=analytic,extreme_tail_guard_stresses=stress,protocol_sha256=sha(OUT/'protocol.json'))
save('result.json',result)
print(json.dumps({key:result[key] for key in ['prespecified_cases','accepted','refused','max_accepted_adaptive_error','max_adaptive_reported_relative_error','conditional_checks','extreme_tail_guard_stresses']},indent=2))
