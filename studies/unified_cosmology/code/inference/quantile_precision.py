"""Conditional Monte Carlo precision of untrimmed weighted marginal quantiles.

Resample contiguous blocks separately within each independently run chain.
This is a numerical diagnostic, not additional astrophysical uncertainty, a
new convergence gate, or assurance that an unvisited posterior mode is absent.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
from scipy.special import logsumexp
from scipy.stats import norm
from scipy.integrate import quad

ROOT=Path(__file__).resolve().parents[4]
PROBABILITIES=np.array([.025,.16,.5,.84,.975])
SEEDS={10:928271,20:928272,500:928273}
DESIGN={'probabilities':PROBABILITIES.tolist(),'default_blocks_per_chain':[10,20],
        'replicates':500,'seeds':SEEDS,'chain_count':4,
        'estimator':'Existing weighted midpoint-CDF interpolation, with explicit duplicated rows.',
        'resampling':'Independently draw the same number of nonoverlapping contiguous blocks within each original chain, with replacement. Preserve row order inside each block.',
        'weights':'Original untrimmed logweights repeated with their rows; normalize anew in every bootstrap replicate.',
        'interpretation':'Conditional numerical precision only. Dependence longer than a block, only four observed chains, sparse tails and unvisited modes are not resolved. No acceptance gate is introduced.'}


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def weighted_quantiles(values,weights):
    order=np.argsort(values)
    x=np.asarray(values)[order];w=weights[order]
    return np.interp(PROBABILITIES,np.cumsum(w)-w/2,x)


def block_quantile_precision(values,logweights,groups,block_counts=(10,20),replicates=500,seeds=None,return_draws=False):
    """Input rows must be in chronological selection order within each chain.

    Equal-sized complete blocks are required, avoiding a random sample-size
    change or discarding tail rows. The default 500 rows/chain permits both
    prescribed partitions. Arbitrary input order is not inferable from arrays.
    """
    lw=np.asarray(logweights,dtype=float);groups=np.asarray(groups)
    if lw.ndim!=1 or groups.shape!=lw.shape or not np.isfinite(lw).all():
        raise ValueError('Finite one-dimensional logweights and matching groups required.')
    arrays={name:np.asarray(x,dtype=float)for name,x in values.items()}
    if not arrays or any(x.shape!=lw.shape or not np.isfinite(x).all()for x in arrays.values()):
        raise ValueError('Every finite marginal array must match the logweights.')
    chain_ids=np.unique(groups)
    if len(chain_ids)!=4:raise ValueError('Four independent parent chains are required.')
    if replicates<2:raise ValueError('At least two bootstrap replicates are required.')
    weights=np.exp(lw-logsumexp(lw))
    baseline={name:weighted_quantiles(x,weights)for name,x in arrays.items()}
    independent_chains={}
    for group in chain_ids:
        keep=groups==group
        chain_weights=np.exp(lw[keep]-logsumexp(lw[keep]))
        independent_chains[str(group)]={'points':int(keep.sum()),
            'global_weight_fraction':float(weights[keep].sum()),
            'raw_weight_ESS':float(1/(chain_weights@chain_weights)),
            'quantiles':{name:weighted_quantiles(x[keep],chain_weights).tolist()for name,x in arrays.items()}}
    out={'status':'conditional_quantile_precision_diagnostic','probabilities':PROBABILITIES.tolist(),
         'points':len(lw),'chains':len(chain_ids),'replicates':replicates,
         'raw_weight_ESS':float(1/(weights@weights)),
         'independent_chain_quantiles':independent_chains,
         'baseline_quantiles':{name:x.tolist()for name,x in baseline.items()},'partitions':{},
         'interpretation':DESIGN['interpretation']}
    all_draws={}
    for count in block_counts:
        if count<2:raise ValueError('At least two blocks per chain are required.')
        blocks=[]
        for group in chain_ids:
            indices=np.flatnonzero(groups==group)
            if len(indices)%count or len(indices)<count:
                raise ValueError('Each chain length must divide exactly into the requested complete blocks; no rows are silently dropped.')
            blocks.append(indices.reshape(count,-1))
        seed=(SEEDS if seeds is None else seeds)[count]
        rng=np.random.default_rng(seed)
        draws={name:np.empty((replicates,len(PROBABILITIES)))for name in arrays}
        ess=np.empty(replicates)
        for b in range(replicates):
            # Keeping explicit copies matches the original midpoint estimator;
            # combining equal-valued rows into one weight would change it.
            sampled=np.concatenate([block[rng.integers(count,size=count)].ravel()for block in blocks])
            w=np.exp(lw[sampled]-logsumexp(lw[sampled]))
            ess[b]=1/(w@w)
            for name,x in arrays.items():draws[name][b]=weighted_quantiles(x[sampled],w)
        out['partitions'][str(count)]={'seed':int(seed),'blocks_per_chain':count,
            'rows_per_block_per_chain':[int(block.shape[1])for block in blocks],
            'bootstrap_raw_weight_ESS_range':[float(ess.min()),float(ess.max())],
            'marginals':{name:{'bootstrap_SD':np.std(x,axis=0,ddof=1).tolist(),
                'bootstrap_quantiles_025_975':np.quantile(x,[.025,.975],axis=0).T.tolist(),
                'bootstrap_mean_minus_original':(x.mean(axis=0)-baseline[name]).tolist()}for name,x in draws.items()}}
        all_draws[str(count)]=draws
    return (out,all_draws)if return_draws else out


def validate():
    """Analytic independent-importance target and stationary correlated mocks."""
    work=ROOT/'.work/unified-cosmology/inference/quantile-precision-validation'
    work.mkdir(parents=True,exist_ok=True)
    groups=np.repeat(np.arange(4),500)
    n=len(groups);mu=.4;sigma=1.25
    target=norm.ppf(PROBABILITIES)
    theoretical=[]
    for p,z in zip(PROBABILITIES,target):
        density=lambda x:np.exp(2*norm.logpdf(x)-norm.logpdf(x,mu,sigma))
        variance=((1-p)**2*quad(density,-np.inf,z,epsabs=1e-10)[0]+p**2*quad(density,z,np.inf,epsabs=1e-10)[0])/(norm.pdf(z)**2*n)
        theoretical.append(np.sqrt(variance))
    rng=np.random.default_rng(928604)
    cohorts=32;rows=[];checks={}
    for rho in [0.,.9]:
        estimates=[];predicted={str(k):[]for k in [10,20,500]}
        for cohort in range(cohorts):
            innovation=rng.normal(size=(4,500))
            x=np.empty_like(innovation);x[:,0]=innovation[:,0]
            for i in range(1,500):x[:,i]=rho*x[:,i-1]+np.sqrt(1-rho*rho)*innovation[:,i]
            x=(mu+sigma*x).ravel()
            lw=norm.logpdf(x)-norm.logpdf(x,mu,sigma)
            report,draws=block_quantile_precision({'x':x},lw,groups,block_counts=(10,20,500),return_draws=True)
            estimates.append(report['baseline_quantiles']['x'])
            for count in [10,20,500]:predicted[str(count)].append(report['partitions'][str(count)]['marginals']['x']['bootstrap_SD'])
            rows.append({'rho':rho,'cohort':cohort,'estimate':report['baseline_quantiles']['x'],'bootstrap_SD':{k:report['partitions'][k]['marginals']['x']['bootstrap_SD']for k in predicted}})
            if rho==0 and cohort==0:
                translated=block_quantile_precision({'x':x},lw+137.,groups)
                shift=max(np.max(abs(np.array(translated['partitions'][str(k)]['marginals']['x']['bootstrap_SD'])-np.array(report['partitions'][str(k)]['marginals']['x']['bootstrap_SD'])))for k in [10,20])
                assert shift<1e-12
                checks['additive_logweight_invariance_max_SD_difference']=float(shift)
                constant=block_quantile_precision({'constant':np.full(n,7.)},lw,groups)
                assert all(np.max(abs(np.array(v['marginals']['constant']['bootstrap_SD'])))==0 for v in constant['partitions'].values())
                checks['constant_quantile_exact_zero_precision_variance']=True
        estimates=np.array(estimates)
        empirical=estimates.std(axis=0,ddof=1)
        checks[str(rho)]={'independent_cohorts':cohorts,'target_quantiles':target.tolist(),
            'mean_estimate':estimates.mean(axis=0).tolist(),'empirical_SD_across_cohorts':empirical.tolist(),
            'mean_bootstrap_SD':{k:np.mean(v,axis=0).tolist()for k,v in predicted.items()},
            'mean_bootstrap_SD_divided_by_empirical':{k:(np.mean(v,axis=0)/empirical).tolist()for k,v in predicted.items()}}
        if rho==.9:
            ratio=np.mean(predicted['10'],axis=0)/np.mean(predicted['500'],axis=0)
            assert np.all(ratio>1.5),ratio
            checks['correlated_10block_vs_single_row_SD_ratio']=ratio.tolist()
        print(json.dumps({'validated_rho':rho,'cohorts':cohorts}),flush=True)
    checks['iid_asymptotic_importance_SD']=theoretical
    checks['iid_mean_bootstrap_SD_divided_by_asymptotic']={k:(np.array(checks['0.0']['mean_bootstrap_SD'][k])/theoretical).tolist()for k in ['10','20','500']}
    rp=work/'mock-cohorts.json';rp.write_text(json.dumps(rows,indent=2)+'\n')
    out={'status':'passed_implementation_and_sensitivity_checks','scientific_cosmology_result':False,
         'design':DESIGN,'mock_design':{'target':'N(0,1)','proposal':'N(0.4,1.25^2)','rows_per_chain':500,'chains':4,'stationary_AR1_rho':[0,.9],'independent_cohorts_per_rho':cohorts,'seed':928604,'replicates_per_partition':500,'single_row_partition':500,'purpose':'Known-target numerical calibration and correlated-mock block-size sensitivity; no gate is chosen from these outcomes.'},
         'checks':checks,'source_sha256':{str(Path(__file__).relative_to(ROOT)):digest(__file__)},'output_sha256':{str(rp.relative_to(ROOT)):digest(rp)}}
    return out


def actual(folder,summary_path):
    # The consumer rechecks native records, parent convergence, all weight gates,
    # dependency bytes, and target/cohort identities before returning.
    from measurement_summary import summarize_run
    folder=folder.resolve();summary_path=summary_path.resolve()
    qualified=summarize_run(folder,summary_path)
    parent_summary=json.loads(summary_path.read_text())
    selection_path=ROOT/parent_summary['selection_path']
    assert selection_path.parent.parent.resolve()==folder
    selection=json.loads(selection_path.read_text())
    groups=np.asarray(selection['groups'])
    for group in np.unique(groups):
        positions=[selection['locations'][i]['expanded_index']for i in np.flatnonzero(groups==group)]
        if not np.all(np.diff(positions)>0):raise ValueError('Rows are not in chronological selection order within each chain.')
    records=[json.loads((selection_path.parent/f'{i:05d}.json').read_text())for i in range(len(groups))]
    names=qualified['weighted_covariance']['parameter_order']
    values={name:np.array([dict(r['point'],**r['derived'])[name]for r in records])for name in names}
    report=block_quantile_precision(values,[r['log_weight']for r in records],groups)
    for name in names:
        assert np.allclose(report['baseline_quantiles'][name],qualified['posterior'][name]['quantiles_025_16_50_84_975'],rtol=0,atol=1e-12)
    report.update(target_identity=qualified['target_identity'],settings=qualified['settings'],
                  parent_measurement_gates_passed=True,design=DESIGN,
                  input_sha256=qualified['input_sha256'],
                  source_sha256={str(p.relative_to(ROOT)):digest(p)for p in [Path(__file__),Path(__file__).with_name('measurement_summary.py')]})
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--validate',action='store_true')
    parser.add_argument('--chain-folder',type=Path)
    parser.add_argument('--correction-summary',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    if args.validate:
        if args.chain_folder or args.correction_summary:parser.error('Validation uses synthetic inputs only.')
        report=validate()
    else:
        if not args.chain_folder or not args.correction_summary:parser.error('A qualified chain folder and correction summary are required.')
        report=actual(args.chain_folder,args.correction_summary)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':report['status'],'output':str(args.output)}),flush=True)


if __name__=='__main__':main()
