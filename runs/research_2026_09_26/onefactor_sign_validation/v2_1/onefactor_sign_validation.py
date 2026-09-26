"""Frozen synthetic validation of a one-factor sign/censor likelihood."""
from pathlib import Path
import hashlib
import itertools
import json
import sys
import numpy as np
from scipy.special import log_ndtr
from scipy.stats import multivariate_normal, norm
import onefactor_sign_likelihood as kernel

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'runs/research_2026_09_26/onefactor_sign_validation'
NATIVE = ROOT / 'runs/research_2026_09_26/astra_design/raisin_signed_refit/fixed-c-profile/global-profile/native-profiles.npz'
SHA = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()


def dump(p, o):
    p.write_text(json.dumps(o, indent=2)+'\n')


def prepare():
    OUT.mkdir(exist_ok=False)
    dump(OUT/'protocol.json', {
        'scope': 'Synthetic mathematical validation only. Native covariance structure already inspected; no observed fluxes fitted or scored.',
        'seed': 120776, 'draws': 100000,
        'inputs_sha256': {str(p.relative_to(ROOT)): SHA(p) for p in [NATIVE, Path(__file__), Path(kernel.__file__)]},
        'checks': ['independent-product limit', 'sum of all 16 four-dimensional orthants',
                   'full Gaussian density against SciPy', 'one missing row versus exact conditional normal CDF',
                   '64/128/256 Gauss-Hermite convergence', 'paired synthetic orthant frequency with six-Monte-Carlo-SE gate',
                   'native C diagonal-plus-one-factor closure, synthetic mean only'],
        'numeric_tolerance': 1e-10, 'native_covariance_relative_tolerance': 1e-12,
        'limitations': 'Valid only for stated Gaussian fixed covariance and known sign observation process. Does not establish real extraction noise, event selection, photometric calibration or cosmology.',
    })


def run():
    p=json.loads((OUT/'protocol.json').read_text())
    for f,h in p['inputs_sha256'].items(): assert SHA(ROOT/f)==h
    mean=np.array([-.5,0.,.5,1.]);d=np.array([1.,.7,1.3,.8]);v=np.array([.3,.2,-.25,.4]);s=np.array([1.,-1.,1.,1.])
    independent=kernel.log_orthant(mean,d,0*v,s)
    independent_error=abs(independent-log_ndtr(s*mean/np.sqrt(d)).sum())
    orthants=np.array([kernel.log_orthant(mean,d,v,np.array(signs),128) for signs in itertools.product([-1.,1.],repeat=4)])
    partition_error=abs(np.exp(orthants).sum()-1)
    y=np.array([.2,-.1,.7,1.2]);C=np.diag(d)+np.outer(v,v)
    gaussian_error=abs(kernel.log_gaussian(y,mean,d,v)-multivariate_normal.logpdf(y,mean=mean,cov=C))
    observed=np.array([1,1,1,0],bool)
    co=C[-1,:3]; cm=mean[-1]+co@np.linalg.solve(C[:3,:3],y[:3]-mean[:3]);cv=C[-1,-1]-co@np.linalg.solve(C[:3,:3],co)
    exact=multivariate_normal.logpdf(y[:3],mean=mean[:3],cov=C[:3,:3])+log_ndtr(-cm/np.sqrt(cv))
    censor_error=abs(kernel.log_censored(y[:3],mean,d,v,observed,[-1.],128)-exact)
    resolution=[kernel.log_orthant(mean,d,v,s,n) for n in [64,128,256]]
    rng=np.random.default_rng(p['seed']);u=rng.normal(size=(p['draws'],1));e=rng.normal(size=(p['draws'],4))
    z=mean+u*v+e*np.sqrt(d);success=np.all(z*s>0,axis=1)
    pred=float(np.exp(resolution[-1]));freq=float(success.mean());se=float(np.sqrt(pred*(1-pred)/len(success)))
    native=np.load(NATIVE);NC=native['covariance_B'];nd,nv=kernel.decompose_covariance(NC)
    native_error=float(np.linalg.norm(NC-np.diag(nd)-np.outer(nv,nv))/np.linalg.norm(NC))
    # Synthetic positive mean; no observed flux or fitted mean enters this check.
    nm=2.5*np.sqrt(np.diag(NC));native_resolution=[kernel.log_orthant(nm,nd,nv,np.ones(len(nd)),n) for n in [64,128,256]]
    nu=rng.normal(size=(p['draws'],1));ne=rng.normal(size=(p['draws'],len(nd)))
    nz=nm+nu*nv+ne*np.sqrt(nd);nf=float(np.all(nz>0,axis=1).mean());npred=float(np.exp(native_resolution[-1]));nse=float(np.sqrt(npred*(1-npred)/len(nz)))
    checks=[independent_error,partition_error,gaussian_error,censor_error,np.ptp(resolution),np.ptp(native_resolution)]
    passed=max(checks)<p['numeric_tolerance'] and abs(freq-pred)<6*se and abs(nf-npred)<6*nse
    result={'passed':bool(passed),'independent_logprob_error':float(independent_error),'orthant_partition_error':float(partition_error),
        'Gaussian_logpdf_error':float(gaussian_error),'single_censored_logdensity_error':float(censor_error),
        'orthant_logprob_64_128_256':resolution,'synthetic_probability':pred,'synthetic_frequency':freq,'MC_SE':se,
        'native_covariance_relative_closure':native_error,'native_synthetic_logprob_64_128_256':native_resolution,
        'native_synthetic_probability':npred,'native_synthetic_frequency':nf,'native_MC_SE':nse,
        'native_max_pair_correlation':float(np.max(np.abs(NC-np.diag(np.diag(NC)))/np.sqrt(np.outer(np.diag(NC),np.diag(NC))))),
        'limitations':p['limitations'],'protocol_sha256':SHA(OUT/'protocol.json')}
    dump(OUT/'result.json',result)
    for path in [Path(__file__),Path(kernel.__file__)]: (OUT/path.name).write_bytes(path.read_bytes())
    dump(OUT/'manifest.json',{'files_sha256':{str(f.relative_to(OUT)):SHA(f) for f in OUT.iterdir() if f.is_file() and f.name!='manifest.json'}})
    print(json.dumps(result,indent=2));assert passed


if __name__=='__main__': {'prepare':prepare,'run':run}[sys.argv[1]]()
