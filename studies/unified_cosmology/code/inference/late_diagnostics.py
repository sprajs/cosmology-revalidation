"""Diagnostics from independent emcee ensembles, not independent-walker fiction."""
import argparse
import json
import hashlib
from pathlib import Path
import numpy as np
import emcee
from diagnostics import metrics, summarize, ROOT


def qj(theta, names, z):
    om = theta[...,names.index('Omega_m')]
    w0 = theta[...,names.index('w0')] if 'w0' in names else -np.ones_like(om)
    wa = theta[...,names.index('wa')] if 'wa' in names else np.zeros_like(om)
    w = w0+wa*z/(1+z)
    m = om*(1+z)**3
    de = (1-om)*np.exp(3*(1+w0+wa)*np.log1p(z)-3*wa*z/(1+z))
    fraction = de/(m+de)
    q = .5+1.5*fraction*w
    j = 1+4.5*fraction*w*(1+w)+1.5*fraction*wa/(1+z)
    return q,j


def inspect(folder,discard=.3):
    folder = folder.resolve()
    fit = json.loads((folder/'fit.json').read_text())
    names = fit['parameter_names']
    files = sorted(folder.glob('ensemble-*.h5'))
    assert len(files)>=4,'Require four independently seeded ensembles.'
    chains = []
    for path in files:
        backend = emcee.backends.HDFBackend(str(path),read_only=True)
        chains.append(backend.get_chain(discard=int(backend.iteration*discard)))
    steps = min(len(x) for x in chains)
    # shape independent ensemble, iteration, walker, parameter
    equal = np.array([x[-steps:] for x in chains])
    derived = []
    derived_names = []
    for z,label in [(0.,'0'),(.5,'05'),(1.,'1')]:
        q,j = qj(equal,names,z)
        derived.extend([q,j]);derived_names.extend(['q'+label,'j'+label])
    values = np.concatenate([equal,np.stack(derived,axis=-1)],axis=-1)
    all_names = names+derived_names
    diagnostic = {}
    posterior = {}
    for index,name in enumerate(all_names):
        posterior[name] = summarize(values[:,:,:,index])
        if np.ptp(values[:,:,:,index])==0:
            diagnostic[name] = {'fixed_by_model':True}
            continue
        # For each fixed walker index, the four series come from four independent
        # ensembles. Never treat the48 mutually interacting walkers as48
        # independent chains, nor flatten them before autocorrelation estimation.
        per_walker = [metrics(values[:,:,walker,index]) for walker in range(values.shape[2])]
        diagnostic[name] = {'maximum_rhat_across_fixed_walker_indices':max(x['rhat'] for x in per_walker),
                            'minimum_bulk_ess_across_indices':min(x['bulk_ess'] for x in per_walker),
                            'minimum_tail_ess_across_indices':min(x['tail_ess'] for x in per_walker),
                            'largest_fixed_walker_mean_mcse':max(x['mean_mcse'] for x in per_walker)}
    failed = {k:v for k,v in diagnostic.items() if not v.get('fixed_by_model') and (v['maximum_rhat_across_fixed_walker_indices']>1.01 or min(v['minimum_bulk_ess_across_indices'],v['minimum_tail_ess_across_indices'])<400)}
    flat = equal.reshape(-1,len(names))
    redshifts = np.linspace(0,2.33,61)
    history = []
    for z in redshifts:
        q,j = qj(flat,names,z)
        history.append({'z':float(z),'q_quantiles_025_16_50_84_975':np.quantile(q,[.025,.16,.5,.84,.975]).tolist(),
                        'j_quantiles_025_16_50_84_975':np.quantile(j,[.025,.16,.5,.84,.975]).tolist()})
    result = {'status':'passed' if not failed else 'not_converged','fit':fit,
              'independent_ensembles':len(chains),'walkers_per_ensemble':values.shape[2],
              'retained_iterations_per_ensemble':steps,'discard_fraction':discard,
              'diagnostics':diagnostic,'failed_gates':failed,'posterior':posterior,'history':history,
              'prior_edge_fractions':{name:float(np.mean((flat[:,i]<bounds[0]+.01*(bounds[1]-bounds[0]))|(flat[:,i]>bounds[1]-.01*(bounds[1]-bounds[0])))) for i,(name,bounds) in enumerate(zip(names,fit['bounds']))},
              'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'inputs':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in files}}
    if 'wa' in names:
        result['fraction_w0_plus_wa_positive'] = float(np.mean(flat[:,names.index('w0')]+flat[:,names.index('wa')]>0))
    result['fraction_q0_nonnegative'] = float(np.mean(values[:,:,:,all_names.index('q0')]>=0))
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument('folder',type=Path)
    p.add_argument('--output',type=Path)
    p.add_argument('--discard',type=float,default=.3)
    a = p.parse_args()
    result = inspect(a.folder,a.discard)
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['history','inputs']},indent=2))


if __name__=='__main__':
    main()
