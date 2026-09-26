"""Synthetic independent mathematical checks; never use observed flux/means."""
from pathlib import Path
import hashlib, importlib.util, itertools, json
import numpy as np
from scipy.integrate import quad
from scipy.special import log_ndtr
from scipy.stats import norm, multivariate_normal

ROOT=Path(__file__).resolve().parents[4];OUT=Path(__file__).resolve().parent
SOURCE=ROOT/'scripts/research_2026_09_26/onefactor_sign_likelihood.py'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
save=lambda n,o:(OUT/n).write_text(json.dumps(o,indent=2,allow_nan=False)+'\n')
spec=importlib.util.spec_from_file_location('reviewed',SOURCE);k=importlib.util.module_from_spec(spec);spec.loader.exec_module(k)
protocol=dict(scope='Synthetic math and saved covariance only. No observed fluxes or fitted means used.',source_sha256=sha(SOURCE),
    checks=['all/none observed limits','conditional likelihood normalization by quadrature','latent loading sign gauge and permutation',
      'native-C synthetic SNR grid [-10,-3,0,3,10] with positive, alternating, first9negative signs;64/128/256 nodes',
      'generic strong-factor one-dimensional exact CDF stress','small/sparse covariance decomposition edge cases'],
    moderate_numeric_tolerance=1e-10,diagnostic_only='Generic strong-factor and extreme-tail stresses are not the native supported model domain.')
save('protocol.json',protocol)

mean=np.array([.3,.8]);d=np.array([.7,1.1]);v=np.array([.25,-.35]);C=np.diag(d)+np.outer(v,v);y=np.array([.4,1.])
lp=k.log_gaussian(y,mean,d,v);g=multivariate_normal.logpdf(y,mean=mean,cov=C)
allobs=k.log_censored(y,mean,d,v,[True,True],[],128)
none=k.log_censored([],mean,d,v,[False,False],[1,-1],128)
orth=k.log_orthant(mean,d,v,[1,-1],128)
onegauss=lambda t:norm.pdf(t,loc=mean[0],scale=np.sqrt(C[0,0]))
condmean=lambda t:mean[1]+C[1,0]/C[0,0]*(t-mean[0])
condvar=C[1,1]-C[1,0]**2/C[0,0]
exact_censor=lambda t:onegauss(t)*norm.cdf(-condmean(t)/np.sqrt(condvar))
pattern_probability=quad(exact_censor,0,np.inf,epsabs=1e-12,epsrel=1e-12)[0]
integral_censor=quad(lambda t:np.exp(k.log_censored([t],mean,d,v,[True,False],[-1],128)),0,np.inf,epsabs=1e-12,epsrel=1e-12)[0]
integral_pattern_conditioned=integral_censor/np.exp(orth)
integral_retained_only=quad(lambda t:np.exp(k.log_positive_conditioned([t],mean[:1],d[:1],v[:1],128)),0,np.inf,epsabs=1e-12,epsrel=1e-12)[0]
gauge=max(abs(k.log_orthant(mean,d,v,[1,-1],128)-k.log_orthant(mean,d,-v,[1,-1],128)),abs(k.log_gaussian(y,mean,d,v)-k.log_gaussian(y,mean,d,-v)))
perm=max(abs(k.log_orthant(mean,d,v,[1,-1],128)-k.log_orthant(mean[::-1],d[::-1],v[::-1],[-1,1],128)),abs(k.log_gaussian(y,mean,d,v)-k.log_gaussian(y[::-1],mean[::-1],d[::-1],v[::-1])))
core=dict(gaussian_error=lp-g,all_observed_error=allobs-lp,none_observed_error=none-orth,integrated_censored_minus_exact=integral_censor-pattern_probability,
    sign_pattern_probability=pattern_probability,orthant_probability=np.exp(orth),pattern_conditioned_integral=integral_pattern_conditioned,retained_only_positive_integral=integral_retained_only,gauge_max_error=gauge,permutation_max_error=perm)
assert max(abs(lp-g),abs(allobs-lp),abs(none-orth),abs(integral_censor-pattern_probability),abs(integral_pattern_conditioned-1),abs(integral_retained_only-1),gauge,perm)<1e-10
nativepath=ROOT/'runs/research_2026_09_26/astra_design/raisin_signed_refit/fixed-c-profile/global-profile/native-profiles.npz'
with np.load(nativepath) as a:NC=a['covariance_B'].copy()
nd,nv=k.decompose_covariance(NC);sig=np.sqrt(np.diag(NC));rows=[]
signs={'positive':np.ones(74),'alternating':np.where(np.arange(74)%2,1.,-1.),'first9negative':np.r_[-np.ones(9),np.ones(65)]}
for snr, (name,s) in itertools.product([-10.,-3.,0.,3.,10.],signs.items()):
    logs=[k.log_orthant(snr*sig,nd,nv,s,order) for order in [64,128,256]]
    rows.append(dict(synthetic_SNR=snr,sign_pattern=name,logprob=logs,spread=max(logs)-min(logs)))
stress=[]
for loading in [.3,3.,1000.]:
    m=-.2;var=1.;exact=float(log_ndtr(m/np.sqrt(var+loading**2)))
    values=[k.log_orthant([m],[var],[loading],[1],order) for order in [64,128,256]]
    stress.append(dict(mean=m,diagonal_variance=var,loading=loading,exact_logprob=exact,logprob=values,errors=[x-exact for x in values]))
edges=[]
for name,cov in [('scalar',np.array([[2.]])),('correlated_two',np.array([[2.,.5],[.5,3.]])),('sparse_three',np.array([[2.,.5,0.],[.5,3.,0.],[0.,0.,1.]]))]:
    try:
        dd,vv=k.decompose_covariance(cov);edges.append(dict(name=name,status='returned',error=float(np.max(abs(cov-np.diag(dd)-np.outer(vv,vv))))))
    except Exception as exc:edges.append(dict(name=name,status='explicit failure',exception=type(exc).__name__,message=str(exc)))
result=dict(scope=protocol['scope'],core=core,native_covariance_only_synthetic_grid=rows,native_grid_max_logprob_spread=max(r['spread'] for r in rows),generic_strong_factor_stress=stress,decomposition_edge_cases=edges,
    source_sha256=sha(SOURCE),native_covariance_file_sha256=sha(nativepath),protocol_sha256=sha(OUT/'protocol.json'))
assert sha(SOURCE)==protocol['source_sha256']
save('result.json',result);print(json.dumps(result,indent=2))
