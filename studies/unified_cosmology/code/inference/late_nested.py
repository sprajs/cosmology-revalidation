"""Independent nested-sampling check of late-time CPL tails under identical priors.

This checks the declared late-time likelihood; it supplies no CMB information,
no empirical host correction and no proof that every possible mode was found.
"""
from __future__ import annotations
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import argparse,hashlib,json,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import logsumexp
import dynesty
from dynesty import utils as dyfunc
from late_geometry import Geometry,ROOT,BAO,SN

WORK=ROOT/'.work/unified-cosmology/inference/late-nested'
RESULTS=ROOT/'studies/unified_cosmology/results/inference'
SEEDS=[927085,927086]

def sha(path):
    with Path(path).open('rb')as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

def design():
    record=dict(scope='Baseline Dovekie+DESI DR2, flat matter+CPL atlate redshift, no radiation and no w0+wa cut; same uniform priors asGeometry. Exploratory independent check following difficult ensemble tail mixing.',
      names=['Omega_m','H0_rdrag','w0','wa'],bounds=[[.01,.99],[5000.,15000.],[-3.,1.],[-5.,3.]],seeds=SEEDS,
      live_points=1000,method='dynesty static nested,multi-ellipsoid bounds,rwalk30',dlogz=.03,maximum_calls=1500000,
      uncertainty='64 prior-volume shrinkage jitter realizations perrun, plus twoindependent live-point populations. Jitter omits missed-mode uncertainty; seedagreement is necessary, not proof of complete modecoverage.',
      diagnostic_gate='Each run weightedESS>=2000, completed before callcap; two-run maximum marginal weighted-CDF distance<=.05 and logZ difference<=3combined reportederrors. These are numerical diagnostic gates, not scientific assumptions or universal convergence guarantees.',
      evidence='logZ is for an unnormalized distance likelihood with analytically marginalized commonSNintercept. Do not compare evidence across changedSNcovariance/data/evolution models without restoring normalizationconstants.',
      frozen_utc=datetime.now(timezone.utc).isoformat())
    p=RESULTS/'late-nested-design.json';RESULTS.mkdir(parents=True,exist_ok=True)
    if p.exists():
        old=json.loads(p.read_text());record['frozen_utc']=old['frozen_utc'];assert old==record
    else:p.write_text(json.dumps(record,indent=2)+'\n')
    return record


def qj(theta,z):
    om,w,wa=theta[:,0],theta[:,2],theta[:,3];wz=w+wa*z/(1+z)
    m=om*(1+z)**3;de=(1-om)*np.exp(3*(1+w+wa)*np.log1p(z)-3*wa*z/(1+z));fraction=de/(m+de)
    return .5+1.5*fraction*wz,1+4.5*fraction*wz*(1+wz)+1.5*fraction*wa/(1+z)


def quantities(theta):
    d={name:theta[:,i]for i,name in enumerate(['Omega_m','H0_rdrag','w0','wa'])}
    d['w0_plus_wa']=theta[:,2]+theta[:,3]
    for z,label in [(0.,'0'),(.5,'05'),(1.,'1')]:d['q'+label],d['j'+label]=qj(theta,z)
    return d


def summary(samples,weights):
    weights=np.asarray(weights);weights=weights/weights.sum();result={}
    for name,x in quantities(samples).items():
        mean=float(weights@x);result[name]=dict(mean=mean,sd=float(np.sqrt(weights@((x-mean)**2))),
          quantiles_025_16_50_84_975=list(map(float,dyfunc.quantile(x,[.025,.16,.5,.84,.975],weights=weights))))
    result['probability_q0_nonnegative']=float(weights@(quantities(samples)['q0']>=0))
    result['probability_w0_plus_wa_positive']=float(weights@(samples[:,2]+samples[:,3]>0))
    return result


class Target:
    def __init__(self):
        self.g=Geometry('cpl','none',True,sample='dovekie');self.bounds=np.asarray(self.g.bounds)
    def prior(self,u):return self.bounds[:,0]+u*(self.bounds[:,1]-self.bounds[:,0])
    def loglike(self,x):return float(self.g.logp(x)[0])


def check_likelihood(target,seed):
    rng=np.random.default_rng(seed);points=[np.array([.31,9900,-.8,-.7])]+[target.prior(rng.uniform(.001,.999,4))for _ in range(3)];checks=[]
    cov=np.loadtxt(BAO/'desi_gaussian_bao_ALL_GCcomb_cov.txt')
    for theta in points:
        om,scale,w,wa=theta
        def E(z):return np.sqrt(om*(1+z)**3+(1-om)*(1+z)**(3*(1+w+wa))*np.exp(-3*wa*z/(1+z)))
        pred=[]
        for z,kind in zip(target.g.bz,target.g.btype):
            dm=299792.458/scale*quad(lambda x:1/E(x),0,z,epsabs=1e-11,epsrel=1e-11)[0];dh=299792.458/scale/E(z)
            pred.append({'DM_over_rs':dm,'DH_over_rs':dh,'DV_over_rs':np.cbrt(dm*dm*dh*z)}[kind])
        r=np.array(pred)-target.g.bmean;directbao=float(r@np.linalg.solve(cov,r));sn,bao=target.g.components(theta);directsn=target.g.direct_chi2(theta)
        assert abs(sn[0]-directsn)<1e-4 and abs(bao[0]-directbao)<1e-7
        checks.append(dict(theta=theta.tolist(),SN_chi2_error=float(sn[0]-directsn),BAO_chi2_error=float(bao[0]-directbao)))
    return checks


def run(seed):
    d=design();WORK.mkdir(parents=True,exist_ok=True);out=WORK/f'seed-{seed}';out.mkdir(exist_ok=True)
    recordfile=RESULTS/f'late-nested-{seed}.json'
    if recordfile.exists():raise RuntimeError('Completed record exists; refuse silent overwrite')
    target=Target();assert target.g.bounds==d['bounds'] and target.g.names==d['names']
    checks=check_likelihood(target,seed);rng=np.random.default_rng(seed);started=time.monotonic()
    sampler=dynesty.NestedSampler(target.loglike,target.prior,4,nlive=d['live_points'],bound='multi',sample='rwalk',walks=30,rstate=rng)
    sampler.run_nested(dlogz=d['dlogz'],maxcall=d['maximum_calls'],print_progress=False,
      checkpoint_file=str(out/'checkpoint.pkl'),checkpoint_every=60)
    r=sampler.results;weight=np.exp(r.logwt-logsumexp(r.logwt));ess=float(1/(weight@weight))
    file=out/'samples.npz';np.savez_compressed(file,samples=r.samples,loglikelihood=r.logl,logwt=r.logwt,weights=weight,logz=r.logz,logzerr=r.logzerr,logvol=r.logvol,ncall=r.ncall)
    post=summary(r.samples,weight);jitter=[]
    for _ in range(64):
        jr=dyfunc.jitter_run(r,rstate=rng);jw=np.exp(jr.logwt-logsumexp(jr.logwt));js=summary(jr.samples,jw)
        jitter.append([js['probability_q0_nonnegative'],js['probability_w0_plus_wa_positive'],js['w0']['quantiles_025_16_50_84_975'][2],js['wa']['quantiles_025_16_50_84_975'][2]])
    jitter=np.asarray(jitter);np.savez_compressed(out/'shrinkage-jitter.npz',values=jitter)
    complete=int(np.sum(r.ncall))<d['maximum_calls'];passed=complete and ess>=2000
    record=dict(status='completed_run' if complete else 'call_cap_incomplete',single_run_diagnostic_passed=passed,seed=seed,live_points=d['live_points'],iterations=int(r.niter),likelihood_calls=int(np.sum(r.ncall)),weighted_ESS=ess,seconds=time.monotonic()-started,
       logZ_unnormalized=float(r.logz[-1]),logZ_error=float(r.logzerr[-1]),posterior=post,
       prior_edge_fractions={name:float(weight@((r.samples[:,i]<lo+.01*(hi-lo))|(r.samples[:,i]>hi-.01*(hi-lo))))for i,(name,(lo,hi))in enumerate(zip(d['names'],d['bounds']))},
       shrinkage_jitter_labels=['P(q0>=0)','P(w0+wa>0)','median_w0','median_wa'],shrinkage_jitter_sd=np.std(jitter,axis=0,ddof=1).tolist(),
       independent_likelihood_checks=checks,environment=dict(dynesty=dynesty.__version__,numpy=np.__version__),
       code_sha256=sha(__file__),dependencies_sha256={str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__).with_name('late_geometry.py'),SN,BAO/'desi_gaussian_bao_ALL_GCcomb_mean.txt',BAO/'desi_gaussian_bao_ALL_GCcomb_cov.txt',RESULTS/'late-nested-design.json']},
       outputs_sha256={str(p.relative_to(ROOT)):sha(p)for p in [file,out/'shrinkage-jitter.npz']},scope=d['scope'],evidence_scope=d['evidence'])
    recordfile.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n');print(json.dumps(record),flush=True)


def compare():
    design();records=[json.loads((RESULTS/f'late-nested-{seed}.json').read_text())for seed in SEEDS];files=[WORK/f'seed-{seed}'/'samples.npz'for seed in SEEDS];a,b=[np.load(p)for p in files];distances={}
    for name in quantities(a['samples']):
        xa,xb=quantities(a['samples'])[name],quantities(b['samples'])[name];ia,ib=np.argsort(xa),np.argsort(xb);grid=np.unique(np.r_[xa,xb])
        ca=np.r_[0,np.cumsum(a['weights'][ia])];cb=np.r_[0,np.cumsum(b['weights'][ib])]
        da=ca[np.searchsorted(xa[ia],grid,side='right')];db=cb[np.searchsorted(xb[ib],grid,side='right')];distances[name]=float(np.max(abs(da-db)))
    error=np.hypot(*[r['logZ_error']for r in records]);zdiff=abs(records[0]['logZ_unnormalized']-records[1]['logZ_unnormalized'])
    passed=all(r['single_run_diagnostic_passed']for r in records) and max(distances.values())<=.05 and zdiff<=3*error
    out=dict(status='passed_independent_seed_diagnostics' if passed else 'diagnostic_failed',seeds=SEEDS,maximum_marginal_CDF_differences=distances,absolute_logZ_difference=zdiff,combined_reported_logZ_error=float(error),
       source_records_sha256={str((RESULTS/f'late-nested-{seed}.json').relative_to(ROOT)):sha(RESULTS/f'late-nested-{seed}.json')for seed in SEEDS},code_sha256=sha(__file__),
       scope='Independentseedagreement supports this bounded posterior calculation; cannot prove absent unidentified or verysmall modes. NoCMB/early-timecut/empiricalhostprior added.')
    (RESULTS/'late-nested-comparison.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--seed',type=int,choices=SEEDS);p.add_argument('--compare',action='store_true');a=p.parse_args()
    if a.compare:compare()
    elif a.seed is not None:run(a.seed)
    else:p.error('Specify --seed or --compare')
