"""Selected-sample conditional magnitude prediction; not parent-population inference.

The training/heldout design and interpretation limits are frozen in
docs/phase2/conditional-prediction-plan.md. No corrected distance enters.
"""
import argparse,datetime,hashlib,json,os,time
os.environ.setdefault('JAX_PLATFORMS','cpu')
os.environ.setdefault('XLA_FLAGS','--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1')
from pathlib import Path
import numpy as np,pandas as pd
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

def predictive(params,data,name,noise):
    indices=jnp.asarray([FEATURES.index(x) for x in MODELS[name]],dtype=int)
    baseline=data['distance_reference']+params['M']+jnp.interp(data['z'],KNOTS,jnp.concatenate([jnp.zeros(1),params['offsets']]))
    coef=params.get('coefficients',jnp.empty(0));lps=[];means=[];variances=[]
    for h in [0,1]:
        f,dx,dc=features(data,h);mu=baseline+f[:,indices]@coef
        v=jnp.stack([jnp.ones(len(mu)),-dx[:,indices]@coef,-dc[:,indices]@coef],axis=-1)
        variance=jnp.einsum('ni,nij,nj->n',v,data['cov'],v)+params['scatter'][h]**2
        law=dist.Normal(mu,jnp.sqrt(variance)) if noise=='gaussian' else dist.StudentT(4.,mu,jnp.sqrt(variance/2.))
        lps.append(law.log_prob(data['y'][:,0]));means.append(mu);variances.append(variance)
    ph=jnp.clip(data['host_prob'],1e-10,1-1e-10)
    lp=jnp.logaddexp(jnp.log1p(-ph)+lps[0],jnp.log(ph)+lps[1])
    mean=(1-ph)*means[0]+ph*means[1]
    var=(1-ph)*variances[0]+ph*variances[1]+ph*(1-ph)*(means[1]-means[0])**2
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

def main():
    p=argparse.ArgumentParser();p.add_argument('--data',required=True);p.add_argument('--name',required=True);p.add_argument('--model',choices=MODELS,required=True);p.add_argument('--noise',choices=['gaussian','student4'],default='gaussian');p.add_argument('--warmup',type=int,default=600);p.add_argument('--draws',type=int,default=1000);p.add_argument('--chains',type=int,default=2);p.add_argument('--seed',type=int,default=2026092160);a=p.parse_args()
    folder=ROOT/a.data;arrays=dict(np.load(folder/'data.npz',allow_pickle=False));rows=pd.read_csv(folder/'rows.csv',dtype={'CID':str})
    # Selection is frozen independently of any model's residual or score.
    keep=(arrays['survey']==10)&(arrays['pIa']>.999)
    if 'common_comparison_member' in arrays:keep &= arrays['common_comparison_member'].astype(bool)
    train=keep&(arrays['fold']!=0);test=keep&(arrays['fold']==0)
    assert train.any() and test.any();out=ROOT/'phase2/hierarchy/conditional'/a.name
    if out.exists():raise RuntimeError('Refusing to replace an existing conditional experiment')
    out.mkdir(parents=True);configuration={'arguments':vars(a),'n_train':int(train.sum()),'n_test':int(test.sum()),'features':MODELS[a.model],'plan':'docs/phase2/conditional-prediction-plan.md','qualification':'Empirical conditional apparent-magnitude prediction in selected high-purity DES sample, from pre-BBC fits; no population/selection-free cosmology/dust identification.'};(out/'configuration.json').write_text(json.dumps(configuration,indent=2)+'\n')
    files=[folder/'data.npz',folder/'rows.csv',ROOT/configuration['plan'],Path(__file__),Path(__file__).with_name('core.py'),ROOT/'phase2/hierarchy/uv.lock'];sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest();frozen={str(f.relative_to(ROOT)):sha(f) for f in files}
    for f in [Path(__file__),Path(__file__).with_name('core.py')]:
        (ROOT/'phase2/hierarchy/provenance'/(sha(f)+'.py')).write_bytes(f.read_bytes())
    def subset(mask):
        d={k:jnp.asarray(v[mask]) for k,v in arrays.items()};d['distance_reference']=distance_kinematic(d['z'],d['zhel'],0.,0.);return d
    dtrain=subset(train);dtest=subset(test);began=time.time()
    mcmc=MCMC(NUTS(model,target_accept_prob=.9,max_tree_depth=10,dense_mass=True,init_strategy=init_to_median(num_samples=20)),num_warmup=a.warmup,num_samples=a.draws,num_chains=a.chains,chain_method='sequential',progress_bar=False)
    mcmc.run(jax.random.PRNGKey(a.seed),dtrain,a.model,a.noise,extra_fields=('diverging','num_steps','accept_prob'))
    samples=mcmc.get_samples(group_by_chain=True);flat=mcmc.get_samples();extra=mcmc.get_extra_fields(group_by_chain=True)
    np.savez_compressed(out/'chains.npz',**{k:np.asarray(v) for k,v in samples.items()})
    fn=jax.jit(jax.vmap(lambda param:predictive(param,dtest,a.model,a.noise)))
    ll,means,variances=fn(flat);score=np.asarray(logsumexp(ll,axis=0)-jnp.log(ll.shape[0]));prediction=np.asarray(means.mean(axis=0));predvar=np.asarray(variances.mean(axis=0)+means.var(axis=0));obs=np.asarray(dtest['y'][:,0])
    table=rows.loc[test,['CID','IDSURVEY','field','depth','fold']].copy();table['observed_mB']=obs;table['predictive_mean']=prediction;table['predictive_sd']=np.sqrt(predvar);table['log_predictive_density']=score;table['standardized_residual']=(obs-prediction)/np.sqrt(predvar);table.to_csv(out/'heldout.csv',index=False)
    # Moment coverage is labelled; Student mixtures do not have Normal quantiles.
    report={'parameters':serial(summary(samples,group_by_chain=True)),'divergences':int(np.asarray(extra['diverging']).sum()),'max_steps':int(np.asarray(extra['num_steps']).max()),'mean_acceptance':float(np.asarray(extra['accept_prob']).mean()),'heldout_log_score_sum':float(score.sum()),'heldout_log_score_mean':float(score.mean()),'heldout_rmse':float(np.sqrt(np.mean((obs-prediction)**2))),'fraction_within_two_predictive_sd':float(np.mean(np.abs(table.standardized_residual)<2)),'runtime_seconds':time.time()-began}
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n');(out/'manifest.json').write_text(json.dumps({'time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'inputs_sha256':frozen,'outputs_sha256':{str(f.relative_to(ROOT)):sha(f) for f in out.iterdir() if f.name!='manifest.json'}},indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='parameters'},indent=2))
if __name__=='__main__':main()
