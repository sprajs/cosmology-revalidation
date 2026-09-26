"""Selected-sample conditional magnitude prediction; not parent-population inference.

The training/heldout design and interpretation limits are frozen in
docs/phase2/conditional-prediction-plan.md. No corrected distance enters.
"""
import argparse,datetime,hashlib,json,os,time
os.environ.setdefault('JAX_PLATFORMS','cpu')
os.environ.setdefault('XLA_FLAGS','--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1')
from pathlib import Path
import numpy as np,pandas as pd
from scipy.special import ndtr
from scipy.stats import t as student_t
import jax,jax.numpy as jnp
from jax.scipy.special import logsumexp
import numpyro,numpyro.distributions as dist
from numpyro.infer import MCMC,NUTS,init_to_median
from numpyro.diagnostics import summary
from core import distance_kinematic

ROOT=Path(__file__).resolve().parents[3]
FEATURES=['x','c','host','host_colour','red_hinge','x_evolution','c_evolution','x_squared']
PRIORS={'x':(-.15,.2),'c':(3.,2.),'host':(0.,.3),'host_colour':(0.,2.),'red_hinge':(0.,2.),'x_evolution':(0.,.5),'c_evolution':(0.,3.),'x_squared':(0.,.1)}
MODELS={'none':[],'stretch':['x'],'colour':['c'],'tripp':['x','c'],'host':['x','c','host'],'host_colour':['x','c','host','host_colour'],'broken_colour':['x','c','host','red_hinge'],'evolution':['x','c','host','x_evolution','c_evolution'],'flexible':FEATURES}
KNOTS=jnp.array([.01,.10,.20,.35,.50,.70,.90,1.20])

def features(data,h):
    x=data['y'][:,1];c=data['y'][:,2];hc=h-.5;g=data['z']/(1+data['z'])-.2;one=jnp.ones_like(x);zero=jnp.zeros_like(x)
    f=jnp.stack([x,c,hc*one,hc*c,jnp.maximum(c,0),x*g,c*g,x*x],axis=-1)
    dx=jnp.stack([one,zero,zero,zero,zero,g,zero,2*x],axis=-1)
    dc=jnp.stack([zero,one,zero,hc*one,(c>0).astype(float),zero,g,zero],axis=-1)
    return f,dx,dc

def components(params,data,name):
    indices=jnp.asarray([FEATURES.index(x) for x in MODELS[name]],dtype=int)
    baseline=data['distance_reference']+params['M']+jnp.interp(data['z'],KNOTS,jnp.concatenate([jnp.zeros(1),params['offsets']]))
    coef=params.get('coefficients',jnp.empty(0));means=[];variances=[]
    for h in [0,1]:
        f,dx,dc=features(data,h);mu=baseline+f[:,indices]@coef
        # Keep the dot products separate from stack. With the pinned JAX 0.6.2,
        # fusing both indexed products into stack gives a different compiled v.
        dmean_dx=dx[:,indices]@coef
        dmean_dc=dc[:,indices]@coef
        v=jnp.stack([jnp.ones(len(mu)),-dmean_dx,-dmean_dc],axis=-1)
        variance=jnp.einsum('ni,nij,nj->n',v,data['cov'],v)+params['scatter'][h]**2
        means.append(mu);variances.append(variance)
    return jnp.stack(means,axis=-1),jnp.stack(variances,axis=-1)

def predictive(params,data,name,noise):
    means,variances=components(params,data,name)
    law=dist.Normal(means,jnp.sqrt(variances)) if noise=='gaussian' else dist.StudentT(4.,means,jnp.sqrt(variances/2.))
    lps=law.log_prob(data['y'][:,0,None])
    ph=jnp.clip(data['host_prob'],1e-10,1-1e-10)
    lp=jnp.logaddexp(jnp.log1p(-ph)+lps[:,0],jnp.log(ph)+lps[:,1])
    mean=(1-ph)*means[:,0]+ph*means[:,1]
    var=(1-ph)*variances[:,0]+ph*variances[:,1]+ph*(1-ph)*(means[:,1]-means[:,0])**2
    return lp,mean,var

def model(data,name,noise):
    p={'M':numpyro.sample('M',dist.Normal(-19.3,1.)),
       'offsets':numpyro.sample('offsets',dist.Normal(0.,.5).expand([7]).to_event(1)),
       'scatter':numpyro.sample('scatter',dist.HalfNormal(.7).expand([2]).to_event(1))}
    if MODELS[name]:
        prior=jnp.asarray([PRIORS[x] for x in MODELS[name]])
        p['coefficients']=numpyro.sample('coefficients',dist.Normal(prior[:,0],prior[:,1]).to_event(1))
    numpyro.factor('conditional_selected_magnitude',predictive(p,data,name,noise)[0].sum())

def serial(x):
    if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(np.ndarray,jax.Array)):return np.asarray(x).tolist()
    if isinstance(x,np.generic):return x.item()
    return x

def convergence_diagnostics(samples,extra,max_tree_depth=10):
    """Reject chains before any held-out prediction is evaluated."""
    stats=summary(samples,group_by_chain=True)
    rhats=np.concatenate([np.asarray(v['r_hat']).ravel() for v in stats.values()])
    neffs=np.concatenate([np.asarray(v['n_eff']).ravel() for v in stats.values()])
    acceptance=np.asarray(extra['accept_prob'])
    steps=np.asarray(extra['num_steps'])
    divergences=int(np.asarray(extra['diverging']).sum())
    reasons=[]
    if not np.all(np.isfinite(rhats)) or np.max(rhats)>1.05:reasons.append('nonfinite or R-hat > 1.05')
    if not np.all(np.isfinite(neffs)) or np.min(neffs)<100:reasons.append('nonfinite or bulk effective sample size < 100')
    if divergences:reasons.append('divergent transitions')
    if not np.all(np.isfinite(acceptance)) or float(acceptance.mean())<.6:reasons.append('nonfinite or mean acceptance < 0.6')
    if np.any(steps>=2**max_tree_depth-1):reasons.append('maximum tree depth reached')
    for name,values in samples.items():
        if not np.all(np.isfinite(values)) or np.any(np.std(values,axis=1)==0):
            reasons.append(f'{name} has a nonfinite or immobile chain')
    report={'valid_for_scoring':not reasons,'reasons':reasons,'max_r_hat':float(np.max(rhats)),
            'min_n_eff':float(np.min(neffs)),'divergences':divergences,
            'max_steps':int(np.max(steps)),'mean_acceptance':float(np.mean(acceptance)),
            'parameters':serial(stats)}
    return report

def main():
    p=argparse.ArgumentParser();p.add_argument('--data',required=True);p.add_argument('--name',required=True);p.add_argument('--model',choices=MODELS,required=True);p.add_argument('--noise',choices=['gaussian','student4'],default='gaussian');p.add_argument('--warmup',type=int,default=600);p.add_argument('--draws',type=int,default=1000);p.add_argument('--chains',type=int,default=2);p.add_argument('--seed',type=int,default=2026092160);p.add_argument('--reference',choices=['coasting','eds','desitter'],default='coasting');p.add_argument('--cov-scale',type=float,default=1.);p.add_argument('--omit-training-cid',default=None);p.add_argument('--progress',action='store_true');p.add_argument('--output-root',default='phase2/hierarchy/conditional');a=p.parse_args()
    folder=ROOT/a.data;arrays=dict(np.load(folder/'data.npz',allow_pickle=False));rows=pd.read_csv(folder/'rows.csv',dtype={'CID':str})
    # Selection is frozen independently of any model's residual or score.
    cohortpath=ROOT/'phase2/hierarchy/conditional-cohort.json';cohort=json.loads(cohortpath.read_text())
    keep=(arrays['survey']==10)&(arrays['pIa']>.999)
    actual={str(cid):int(fold) for cid,fold in zip(rows.loc[keep,'CID'],arrays['fold'][keep])}
    assert actual==cohort['cid_to_fold'],'Every arm must preserve exactly the frozen CID/fold cohort'
    arrays['cov']=arrays['cov']*a.cov_scale
    train=keep&(arrays['fold']!=0);test=keep&(arrays['fold']==0)
    if a.omit_training_cid:
        hit=rows.CID==a.omit_training_cid;assert hit.sum()==1 and bool(train[hit][0]);train &= ~hit.to_numpy()
    assert train.any() and test.any();out=ROOT/a.output_root/a.name
    if out.exists():raise RuntimeError('Refusing to replace an existing conditional experiment')
    out.mkdir(parents=True);configuration={'arguments':vars(a),'n_train':int(train.sum()),'n_test':int(test.sum()),'features':MODELS[a.model],'plan':'docs/phase2/conditional-prediction-plan.md','qualification':'Empirical conditional apparent-magnitude prediction in selected high-purity DES sample, from pre-BBC fits; no population/selection-free cosmology/dust identification.'};(out/'configuration.json').write_text(json.dumps(configuration,indent=2)+'\n')
    files=[folder/'data.npz',folder/'rows.csv',cohortpath,ROOT/configuration['plan'],Path(__file__),Path(__file__).with_name('core.py'),ROOT/'phase2/hierarchy/uv.lock'];sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest();frozen={str(f.relative_to(ROOT)):sha(f) for f in files}
    for f in [Path(__file__),Path(__file__).with_name('core.py')]:
        (ROOT/'phase2/hierarchy/provenance'/(sha(f)+'.py')).write_bytes(f.read_bytes())
    def subset(mask):
        d={k:jnp.asarray(v[mask]) for k,v in arrays.items()};qref={'coasting':0.,'eds':.5,'desitter':-1.}[a.reference];d['distance_reference']=distance_kinematic(d['z'],d['zhel'],qref,0.);return d
    dtrain=subset(train);dtest=subset(test);began=time.time()
    mcmc=MCMC(NUTS(model,target_accept_prob=.9,max_tree_depth=10,dense_mass=True,init_strategy=init_to_median(num_samples=20)),num_warmup=a.warmup,num_samples=a.draws,num_chains=a.chains,chain_method='sequential',progress_bar=a.progress)
    mcmc.run(jax.random.PRNGKey(a.seed),dtrain,a.model,a.noise,extra_fields=('diverging','num_steps','accept_prob'))
    samples=mcmc.get_samples(group_by_chain=True);extra=mcmc.get_extra_fields(group_by_chain=True)
    np.savez_compressed(out/'chains.npz',**{k:np.asarray(v) for k,v in samples.items()})
    gate=convergence_diagnostics(samples,extra)
    (out/'convergence.json').write_text(json.dumps(gate,indent=2)+'\n')
    if not gate['valid_for_scoring']:
        (out/'INVALID_FOR_INFERENCE.json').write_text(json.dumps({'status':'failed_convergence','reasons':gate['reasons'],'instruction':'No held-out scores were calculated; diagnose sampler and rerun under a new name.'},indent=2)+'\n')
        (out/'manifest.json').write_text(json.dumps({'time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'inputs_sha256':frozen,'outputs_sha256':{str(f.relative_to(ROOT)):sha(f) for f in out.iterdir() if f.name!='manifest.json'}},indent=2)+'\n')
        raise RuntimeError(f"Convergence gate failed before held-out scoring: {gate['reasons']}")
    flat=mcmc.get_samples()
    fn=jax.jit(jax.vmap(lambda param:predictive(param,dtest,a.model,a.noise)))
    ll,means,variances=fn(flat);score=np.asarray(logsumexp(ll,axis=0)-jnp.log(ll.shape[0]));prediction=np.asarray(means.mean(axis=0));predvar=np.asarray(variances.mean(axis=0)+means.var(axis=0));obs=np.asarray(dtest['y'][:,0])
    table=rows.loc[test,['CID','IDSURVEY','field','depth','fold']].copy();table['observed_mB']=obs;table['predictive_mean']=prediction;table['predictive_sd']=np.sqrt(predvar);table['log_predictive_density']=score;table['standardized_residual']=(obs-prediction)/np.sqrt(predvar);table.to_csv(out/'heldout.csv',index=False)
    cm,cv=jax.jit(jax.vmap(lambda param:components(param,dtest,a.model)))(flat);cm=np.asarray(cm);cv=np.asarray(cv);ph=np.clip(np.asarray(dtest['host_prob']),1e-10,1-1e-10);weights=np.stack([1-ph,ph],axis=-1)[None]
    def cdf(values):
        t=(values[None,:,None]-cm)/np.sqrt(cv if a.noise=='gaussian' else cv/2.)
        probs=ndtr(t) if a.noise=='gaussian' else student_t.cdf(t,4.)
        return np.sum(weights*probs,axis=-1).mean(axis=0)
    def quantile(prob):
        lo=prediction-100*np.sqrt(predvar);hi=prediction+100*np.sqrt(predvar)
        for _ in range(50):
            mid=(lo+hi)/2;left=cdf(mid)<prob;lo=np.where(left,mid,lo);hi=np.where(left,hi,mid)
        return (lo+hi)/2
    table['predictive_cdf_at_observation']=cdf(obs);table['predictive_q025']=quantile(.025);table['predictive_q975']=quantile(.975)
    # Batch estimates expose sampling uncertainty in the integrated log score.
    lla=np.asarray(ll).reshape(a.chains,a.draws,len(obs));batch_scores=[]
    for chain in lla:
        for batch in np.array_split(chain,10):batch_scores.append(np.asarray(logsumexp(jnp.asarray(batch),axis=0)-np.log(len(batch))))
    batch_scores=np.asarray(batch_scores);table['log_score_batch_mcse']=batch_scores.std(axis=0,ddof=1)/np.sqrt(len(batch_scores));table.to_csv(out/'heldout.csv',index=False);np.savez_compressed(out/'predictive_batch_scores.npz',batch_scores=batch_scores)
    # Moment coverage is separately labelled; exact mixture quantiles are above.
    report={'parameters':serial(summary(samples,group_by_chain=True)),'divergences':int(np.asarray(extra['diverging']).sum()),'max_steps':int(np.asarray(extra['num_steps']).max()),'mean_acceptance':float(np.asarray(extra['accept_prob']).mean()),'heldout_log_score_sum':float(score.sum()),'heldout_log_score_mean':float(score.mean()),'heldout_rmse':float(np.sqrt(np.mean((obs-prediction)**2))),'fraction_within_two_predictive_sd':float(np.mean(np.abs(table.standardized_residual)<2)),'runtime_seconds':time.time()-began}
    report.update(heldout_95_interval_coverage=float(np.mean((obs>=table.predictive_q025)&(obs<=table.predictive_q975))),heldout_log_score_batch_mcse=float(batch_scores.sum(axis=1).std(ddof=1)/np.sqrt(len(batch_scores))),sampling_uncertainty_qualification='20 within-chain batches; approximate MCSE, not sampling uncertainty across supernovae or fields.')
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n');(out/'manifest.json').write_text(json.dumps({'time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'inputs_sha256':frozen,'outputs_sha256':{str(f.relative_to(ROOT)):sha(f) for f in out.iterdir() if f.name!='manifest.json'}},indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='parameters'},indent=2))
if __name__=='__main__':main()
