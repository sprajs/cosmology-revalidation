"""Fit saved independent hierarchical configurations with NUTS diagnostics."""
import argparse,datetime,hashlib,json,os,time
os.environ.setdefault('JAX_PLATFORMS','cpu')
os.environ.setdefault('XLA_FLAGS','--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1')
from pathlib import Path
import numpy as np
import jax,jax.numpy as jnp
from numpyro.infer import MCMC,NUTS,init_to_median
from numpyro.diagnostics import summary
from inference import model,MODELS

ROOT=Path(__file__).resolve().parents[3]
def main():
    p=argparse.ArgumentParser();p.add_argument('--data',required=True);p.add_argument('--name',required=True);p.add_argument('--model',choices=list(MODELS),default='D0');p.add_argument('--distance',choices=['fixed','kinematic','qbins','offsetbins'],default='fixed');p.add_argument('--selection',choices=['complete','synthetic_probit'],default='complete');p.add_argument('--warmup',type=int,default=800);p.add_argument('--draws',type=int,default=1200);p.add_argument('--chains',type=int,default=2);p.add_argument('--seed',type=int,default=2026092103);a=p.parse_args()
    out=ROOT/'phase2/hierarchy/fits'/a.name;out.mkdir(parents=True,exist_ok=True);(out/'configuration.json').write_text(json.dumps(vars(a),indent=2)+'\n')
    files=[ROOT/a.data,ROOT/'phase2/hierarchy/uv.lock',ROOT/'docs/phase2/plan.md',*Path(__file__).parent.glob('*.py')]
    frozen={str(x.relative_to(ROOT)):hashlib.sha256(x.read_bytes()).hexdigest() for x in files};snap=ROOT/'phase2/hierarchy/provenance';snap.mkdir(exist_ok=True)
    for f in Path(__file__).parent.glob('*.py'):(snap/(frozen[str(f.relative_to(ROOT))]+'.py')).write_bytes(f.read_bytes())
    data={k:jnp.asarray(v) for k,v in np.load(ROOT/a.data,allow_pickle=False).items()};data['selection']=a.selection
    if a.distance=='qbins':
        z=np.asarray(data['z']);zh=np.asarray(data['zhel'])
        if np.any(~np.isfinite(z)) or np.any((z<=0)|(z>1.3)) or np.any(~np.isfinite(zh)) or np.any(zh<=-1):
            raise ValueError('qbins requires finite 0 < z <= 1.3 and zhel > -1; no extrapolation is defined')
    started=time.time();sampler=MCMC(NUTS(model,target_accept_prob=.9,max_tree_depth=10,dense_mass=True,init_strategy=init_to_median(num_samples=20)),num_warmup=a.warmup,num_samples=a.draws,num_chains=a.chains,chain_method='sequential',progress_bar=False)
    sampler.run(jax.random.PRNGKey(a.seed),data,a.model,a.distance,extra_fields=('diverging','num_steps','accept_prob','potential_energy'))
    samples={k:np.asarray(v) for k,v in sampler.get_samples(group_by_chain=True).items()};extra={k:np.asarray(v) for k,v in sampler.get_extra_fields(group_by_chain=True).items()}
    np.savez_compressed(out/'chains.npz',**samples);np.savez_compressed(out/'sampler.npz',**extra)
    stats=summary(samples,group_by_chain=True)
    def serial(v):
        if isinstance(v,dict):return {k:serial(x) for k,x in v.items()}
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        return v
    report={'parameters':serial(stats),'divergences':int(extra['diverging'].sum()),'maximum_num_steps':int(extra['num_steps'].max()),'mean_acceptance':float(extra['accept_prob'].mean()),'runtime_seconds':time.time()-started,'qualification':'Independent effective-SALT hierarchical model; synthetic selection explicitly labelled; no automatic claim of physical dust identification.'}
    if 'q0' in samples:report['P_q0_negative_sample_fraction']=float((samples['q0']<0).mean())
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
    (out/'manifest.json').write_text(json.dumps({'time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'purpose':'Independent hierarchical inference','configuration':vars(a),'inputs_sha256':frozen,'code_capture':'At process start; exact snapshots under phase2/hierarchy/provenance','outputs_sha256':{str(f.relative_to(ROOT)):sha(f) for f in out.iterdir() if f.name!='manifest.json'}},indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
