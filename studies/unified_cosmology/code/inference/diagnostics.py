"""Independent rank-normalised chain diagnostics and transparent posterior summaries."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import arviz as az

ROOT = Path(__file__).resolve().parents[4]
PARAMETERS = ['H0','ombh2','omch2','logA','ns','tau','w','wa','epsilon',
              'omegam','rdrag','q0','q05','q1','j0','j05','j1','sn_chi2']


def summarize(samples):
    x = np.asarray(samples).ravel()
    return {'mean':float(x.mean()),'sd':float(x.std(ddof=1)),
            'quantiles_025_16_50_84_975':np.quantile(x,[.025,.16,.5,.84,.975]).tolist(),
            'range':[float(x.min()),float(x.max())]}


def metrics(values):
    # Input chains are independent MPI chains, with Metropolis frequency weights
    # expanded to their original integer multiplicities. No iid resampling.
    return {'rhat':float(az.rhat(values,method='rank')),
            'bulk_ess':float(az.ess(values,method='bulk')),
            'tail_ess':float(az.ess(values,method='tail')),
            'mean_mcse':float(az.mcse(values,method='mean'))}


def check_mpi(folder,discard=.3):
    folder = folder.resolve()
    paths = sorted(folder.glob('chain.[0-9]*.txt'))
    assert len(paths)>=4,'At least four independent chains are required.'
    names = paths[0].open().readline().lstrip('#').split()
    chains = []
    for path in paths:
        assert names==path.open().readline().lstrip('#').split()
        raw = np.loadtxt(path,ndmin=2)
        weights = raw[:,0]
        assert np.allclose(weights,np.round(weights)) and np.all(weights>0)
        expanded = np.repeat(raw,weights.astype(int),axis=0)
        chains.append(expanded[int(discard*len(expanded)):])
    length = min(map(len,chains))
    assert length>=100
    equal = np.array([x[-length:] for x in chains])
    diagnostic = {}
    posterior = {}
    # Summaries use exactly the same final common chain segments diagnosed below.
    allrows = equal.reshape(-1,len(names))
    sampled = []
    # Updated YAML identifies true sampled parameters; nuisance variables must
    # pass too, not just the cosmological subset selected for the manuscript.
    from cobaya.yaml import yaml_load_file
    config = yaml_load_file(str(folder/'chain.updated.yaml'))
    for key,value in config['params'].items():
        if isinstance(value,dict) and 'prior' in value:
            sampled.append(key)
    for key in sorted(set(sampled+PARAMETERS)&set(names)):
        i = names.index(key)
        if np.ptp(equal[:,:,i])==0:
            continue
        diagnostic[key] = metrics(equal[:,:,i])
        posterior[key] = summarize(allrows[:,i])
    failed = {k:v for k,v in diagnostic.items()
              if not (v['rhat']<=1.01 and min(v['bulk_ess'],v['tail_ess'])>=400)}
    failures = {}
    for path in folder.glob('camb-failures-*.jsonl'):
        failures[path.name] = len(path.read_text().splitlines())
    qprob = {}
    for key in ['q0','q05','q1','j0']:
        if key in names:
            qprob[key] = {'fraction_below_zero':float(np.mean(allrows[:,names.index(key)]<0)),
                          'warning':'Finite Monte Carlo fraction, not a certainty statement when no draws cross zero.'}
    return {'status':'passed' if not failed else 'not_converged','discard_fraction_after_sampler_burnin':discard,
            'independent_chains':len(chains),'retained_weighted_draws':list(map(len,chains)),
            'equal_chain_length_for_diagnostics':length,'diagnostics':diagnostic,'failed_gates':failed,
            'posterior':posterior,'sign_fractions':qprob,'camb_failure_records':failures,
            'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
            'diagnostic_source':'https://arxiv.org/abs/1903.08008',
            'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('folder',type=Path)
    p.add_argument('--discard',type=float,default=.3)
    p.add_argument('--output',type=Path)
    a = p.parse_args()
    result = check_mpi(a.folder,a.discard)
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
