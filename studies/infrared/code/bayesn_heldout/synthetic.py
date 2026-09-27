"""Prepare all twelve synthetic datasets; run one declared case at a time.

Run with the archived isolated BayeSN Python. No observed data fit exists here.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import datetime
import json
import multiprocessing
import os
from pathlib import Path
import sys
import time

from model import ARCHIVE, GRAY, HERE, ROOT, RUN, Kernel, check_parent, configure_cpu, dump, filter_config, numpyro_model, sha, source_hashes


def prepare(work, archive, validation):
    configure_cpu()
    import numpy as np
    work.mkdir(parents=True, exist_ok=True)
    assert not (work/'prepared.json').exists()
    design, gate = check_parent(archive)
    checked = json.loads(validation.read_text())
    assert checked['status']=='passed' and checked['source_sha256']==source_hashes()
    config = filter_config(archive,work)
    directory = archive/RUN/'astra_design/bayesn_signed_pilot'
    cohort = {r['CID']: r for r in csv.DictReader((directory/'frozen-cohort.csv').open())}
    mask = list(csv.DictReader((directory/'frozen-mask-metadata.csv').open()))
    preparation = {'scope':'Synthetic only; no measured flux parsed.', 'source_sha256':source_hashes(),
                   'design_sha256':sha(HERE/'design.json'), 'parent_gate_sha256':sha(ROOT/'studies/infrared/results/bayesn-nir-forward-gate.json'),
                   'archive':str(archive), 'filter_config_sha256':sha(config),
                   'validation_path':str(validation.resolve()),'validation_sha256':sha(validation),
                   'cases':[], 'forward_closure':[]}
    dump(work/'preparation-plan.json', {k:v for k,v in preparation.items() if k!='cases'})
    for object_index,cid in enumerate(design['objects']):
        c = cohort[cid]
        metadata = {'CID':cid}|{k:float(c[k]) for k in ('zHEL','zHD','MWEBV','optical_trigger','mu_LCDM','sigma_external')}
        rows = []
        for r in mask:
            if r['CID'] != cid or r['primary_keep'] != 'True':
                continue
            path = archive/r['source_path']
            # Parse only metadata and the quoted error. Never convert FLUXCAL.
            with path.open() as stream:
                for line_number,line in enumerate(stream,1):
                    if line_number == int(r['source_line']):
                        tokens = line.split()
                        assert tokens[0]=='OBS:' and float(tokens[1])==float(r['MJD']) and tokens[2]==r['band']
                        error = float(tokens[5]); assert error>0
                        break
                else:
                    raise ValueError('Missing source row')
            rows.append(dict(role=r['role'],band=r['band'],trigger_rest_time=float(r['trigger_rest_time']),
                             MJD=float(r['MJD']),error=error,source_path=r['source_path'],source_line=int(r['source_line'])))
        rows.sort(key=lambda r:(r['role']!='optical',r['source_path'],r['source_line']))
        nopt = sum(r['role']=='optical' for r in rows)
        assert nopt == (15 if object_index==0 else 14)
        assert len(rows)-nopt == (6 if object_index==0 else 5)
        kernel = Kernel(archive,config,metadata,rows,design['integration_bins'])
        errors = np.array([r['error'] for r in rows])
        for truth_index,truth in enumerate(design['truths']):
            parameters = np.zeros(47)
            parameters[:4] = [metadata['mu_LCDM']+truth['gray_offset'],truth['AV'],truth['RV'],truth['theta']]
            parameters[-1] = design['truth_tau']
            if truth['epsilon_seed'] is not None:
                parameters[4:46] = np.random.default_rng(truth['epsilon_seed']).normal(size=(2,42))[object_index]
            means = np.asarray(kernel.forward(parameters))
            assert np.isfinite(means).all()
            # The new one-object layout also closes against the official API.
            middle = np.asarray(kernel.model.L_Sigma)@parameters[4:46]
            eps = np.zeros((1,9,6)); eps[0,1:-1,:] = middle.reshape(7,6,order='F')
            reference,_,_ = kernel.model.simulate_light_curve(
                np.array([r['trigger_rest_time'] for r in rows]),1,
                ['RAISIN_DES_'+r['band'] for r in rows],yerr=0.,err_type='flux',
                z=metadata['zHEL'],mu=parameters[0],ebv_mw=metadata['MWEBV'],RV=parameters[2],
                tmax=parameters[46],del_M=0.,AV=parameters[1],theta=parameters[3],eps=eps,mag=False)
            gap = float(np.max(abs(means-np.asarray(reference)[:,0])/errors))
            assert gap < 1e-8
            preparation['forward_closure'].append(dict(CID=cid,truth_index=truth_index,max_quoted_error=gap))
            for noise_seed in design['noise_seeds']:
                index = len(preparation['cases'])
                case_dir = work/f'case-{index:02d}'
                case_dir.mkdir(exist_ok=True)
                for role_index,(name,slc) in enumerate([('optical',slice(0,nopt)),('synthetic-nir',slice(nopt,None))]):
                    rng = np.random.default_rng(np.random.SeedSequence([noise_seed,object_index,truth_index,role_index]))
                    flux = means[slc]+errors[slc]*rng.normal(size=len(means[slc]))
                    dataset = dict(scope='synthetic',case_index=index,metadata=metadata,rows=rows[slc],
                                   flux=flux.tolist(),errors=errors[slc].tolist(),
                                   source_error_semantics='Original quoted errors; generated signed Gaussian measurements, no sign cut.')
                    dump(case_dir/f'{name}.json',dataset)
                dump(case_dir/'truth.json',dict(case_index=index,object_index=object_index,truth_index=truth_index,
                     noise_seed=noise_seed,parameters=parameters.tolist(),noiseless_optical=means[:nopt].tolist(),noiseless_nir=means[nopt:].tolist()))
                preparation['cases'].append(dict(case_index=index,CID=cid,truth_index=truth_index,noise_seed=noise_seed,
                    optical_sha256=sha(case_dir/'optical.json'),nir_sha256=sha(case_dir/'synthetic-nir.json'),truth_sha256=sha(case_dir/'truth.json')))
        print(json.dumps({'prepared_object':cid,'cases':len(preparation['cases'])}),flush=True)
    assert len(preparation['cases'])==12
    check_parent(archive)
    assert preparation['source_sha256']==source_hashes()
    dump(work/'prepared.json',preparation)


def worker_spec(work, case_index, arm, cpu, preparation):
    """The optical worker identity deliberately contains no NIR payload hash."""
    case = preparation['cases'][case_index]
    return dict(work=str(work),case_index=case_index,arm=arm,cpu=cpu,
                optical_sha256=case['optical_sha256'],archive=preparation['archive'],
                source_sha256=preparation['source_sha256'],design_sha256=preparation['design_sha256'],
                filter_config_sha256=preparation['filter_config_sha256'])


def fit_arm(spec):
    configure_cpu(spec['cpu'])
    import numpy as np
    import jax
    import numpyro
    from numpyro.infer import MCMC,NUTS,init_to_value
    work = Path(spec['work'])
    case_dir = work/f"case-{spec['case_index']:02d}"
    output = case_dir/spec['arm']
    output.mkdir(exist_ok=True)
    assert not (output/'worker-start.json').exists(), 'No automatic retry of a partially executed worker.'
    assert source_hashes()==spec['source_sha256']
    assert sha(HERE/'design.json')==spec['design_sha256']
    design = json.loads((HERE/'design.json').read_text())
    data_path = case_dir/'optical.json'
    assert sha(data_path)==spec['optical_sha256']
    data = json.loads(data_path.read_text())
    assert data['scope']=='synthetic' and all(r['role']=='optical' for r in data['rows'])
    config = work/'release-filter-config.yaml'
    assert sha(config)==spec['filter_config_sha256']
    def isolate(event,arguments):
        if event != 'open':
            return
        path = str(arguments[0]).lower()
        if 'synthetic-nir.json' in path or 'nir-sealed' in path or '/photometry/' in path or path.endswith('truth.json'):
            raise RuntimeError('Optical worker attempted to open held-out outcomes or source photometry')
    sys.addaudithook(isolate)
    dump(output/'worker-start.json',spec|{'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
    start = time.monotonic()
    kernel = Kernel(Path(spec['archive']),config,data['metadata'],data['rows'],design['integration_bins'])
    likelihood = numpyro_model(kernel,data,spec['arm'])
    chain_records = []
    for chain_index in range(4):
        arm_index = 0 if spec['arm']=='LCDM' else 1
        seed = design['sampler']['chain_seeds'][chain_index]+100*spec['case_index']+10*arm_index
        initial = dict(D=data['metadata']['mu_LCDM'],AV=.2,RV=3.1,theta=0.,epsilon_white=np.zeros(42),tau=design['sampler']['initial_tau'][chain_index])
        sampler = MCMC(NUTS(likelihood,init_strategy=init_to_value(values=initial),target_accept_prob=.9,max_tree_depth=10,dense_mass=False),
                       num_warmup=1000,num_samples=1000,num_chains=1,progress_bar=False)
        tic = time.monotonic()
        sampler.run(jax.random.PRNGKey(seed),extra_fields=('energy','potential_energy','num_steps','diverging','accept_prob'))
        samples = {k:np.asarray(v) for k,v in sampler.get_samples().items()}
        fields = {k:np.asarray(v) for k,v in sampler.get_extra_fields().items()}
        path = output/f'chain-{chain_index}.npz'
        assert not path.exists()
        np.savez_compressed(path,**samples,**{f'stat_{k}':v for k,v in fields.items()})
        record = dict(chain=chain_index,seed=seed,warmup=1000,draws=1000,seconds=time.monotonic()-tic,
                      archive=path.name,sha256=sha(path),initial_tau=initial['tau'],
                      divergences=int(fields['diverging'].sum()),max_num_steps=int(fields['num_steps'].max()))
        dump(output/f'chain-{chain_index}.json',record)
        chain_records.append(record)
        print(json.dumps({'case':spec['case_index'],'arm':spec['arm'],**record}),flush=True)
        if time.monotonic()-start > design['resources']['max_wall_seconds_per_case']:
            raise TimeoutError('Preserve full finished chains; declared case wall cap exceeded')
    assert source_hashes()==spec['source_sha256'] and sha(data_path)==spec['optical_sha256']
    result = dict(scope='Synthetic optical posterior only; NIR file access blocked.',worker_spec=spec,chains=chain_records,
                  seconds=time.monotonic()-start,versions=dict(numpy=np.__version__,jax=jax.__version__,numpyro=numpyro.__version__))
    dump(output/'posterior-manifest.json',result)
    return result


def diagnose(work,case_index,arm):
    import numpy as np
    import arviz as az
    from scipy.stats import truncnorm
    folder = work/f'case-{case_index:02d}'/arm
    manifest = json.loads((folder/'posterior-manifest.json').read_text())
    assert manifest['worker_spec']['source_sha256']==source_hashes()
    files = []
    for r in manifest['chains']:
        p = folder/r['archive']; assert sha(p)==r['sha256']; files.append(np.load(p))
    params = ('D','AV','RV','theta','epsilon_white','tau')
    posterior = {key:np.stack([x[key] for x in files]) for key in params}
    energy = np.stack([x['stat_energy'] for x in files])
    summary = az.summary(az.from_dict(posterior=posterior),kind='all',round_to='none')
    g = json.loads((HERE/'design.json').read_text())['diagnostics']
    rhat = float(summary.r_hat.max()); bulk = float(summary.ess_bulk.min()); tail = float(summary.ess_tail.min())
    mcse = float((summary.mcse_mean/summary.sd).max())
    bfmi = np.mean(np.diff(energy,axis=1)**2,axis=1)/np.var(energy,axis=1,ddof=1)
    divergences = int(sum(x['stat_diverging'].sum() for x in files))
    tree = [float(np.mean(x['stat_num_steps']>=2**10-1)) for x in files]
    finite = bool(np.isfinite(summary[['r_hat','ess_bulk','ess_tail','mcse_mean','sd']].to_numpy()).all())
    passed = finite and rhat<g['Rhat_max'] and bulk>=g['bulk_ESS_min'] and tail>=g['tail_ESS_min'] and mcse<g['mean_MCSE_over_SD_max'] and divergences==0 and float(bfmi.min())>g['EBFMI_min'] and max(tree)<g['max_tree_fraction_max']
    tau = posterior['tau']; timing = float(np.mean((tau<-9)|(tau>19)))
    broad = None
    if arm=='broad':
        d = posterior['D']; a,b=(20-d)/GRAY,(50-d)/GRAY
        broad = dict(lower=float(np.mean(truncnorm.cdf(21,a,b,loc=d,scale=GRAY))),
                     upper=float(np.mean(truncnorm.sf(49,a,b,loc=d,scale=GRAY))))
    result = dict(sampler_passed=bool(passed),case_index=case_index,arm=arm,
                  posterior_manifest_sha256=sha(folder/'posterior-manifest.json'),Rhat_max=rhat,bulk_ESS_min=bulk,tail_ESS_min=tail,
                  mean_MCSE_over_SD_max=mcse,EBFMI_by_chain=bfmi.tolist(),divergences=divergences,max_tree_fraction_by_chain=tree,
                  timing_boundary_probability=timing,timing_limited=timing>.01,broad_mu_boundary_probabilities=broad,
                  scalar_parameters=summary.reset_index(names='parameter').to_dict(orient='records'),
                  interpretation='Monte Carlo diagnostics for one synthetic realization; no coverage claim or observed fit.')
    dump(folder/'diagnostics.json',result)
    return result


def run_case(work,case_index,workers):
    assert workers in (1,2)
    plan = json.loads((work/'prepared.json').read_text())
    assert plan['source_sha256']==source_hashes()
    assert sha(plan['validation_path'])==plan['validation_sha256']
    checked = json.loads(Path(plan['validation_path']).read_text())
    assert checked['status']=='passed' and checked['source_sha256']==source_hashes()
    assert 0<=case_index<12
    check_parent(Path(plan['archive']))
    # The real-NIR parent hashes are checked only here, not in fit workers.
    # No synthetic held-out value or hash is a worker dependency.
    cpus = sorted(os.sched_getaffinity(0))
    specs = [worker_spec(work,case_index,arm,cpus[i%workers],plan) for i,arm in enumerate(('LCDM','broad'))]
    dump(work/f'case-{case_index:02d}'/'fit-plan.json',{'source_sha256':source_hashes(),'worker_specs':specs})
    with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as pool:
        list(pool.map(fit_arm,specs))
    check_parent(Path(plan['archive']))
    for arm in ('LCDM','broad'):
        print(json.dumps({k:v for k,v in diagnose(work,case_index,arm).items() if k!='scalar_parameters'}),flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('prepare','run-case','diagnose'))
    parser.add_argument('--work',type=Path,default=ROOT/'.work/infrared/heldout-synthetic')
    parser.add_argument('--archive',type=Path,default=ARCHIVE)
    parser.add_argument('--case',type=int,default=0)
    parser.add_argument('--workers',type=int,choices=(1,2),default=2)
    parser.add_argument('--arm',choices=('LCDM','broad'),default='LCDM')
    parser.add_argument('--validation',type=Path,default=ROOT/'.work/infrared/heldout-synthetic-validation-final.json')
    args = parser.parse_args()
    args.work = args.work.resolve(); args.archive = args.archive.resolve()
    if args.action=='prepare':
        prepare(args.work,args.archive,args.validation)
    elif args.action=='run-case':
        run_case(args.work,args.case,args.workers)
    else:
        print(json.dumps(diagnose(args.work,args.case,args.arm)),flush=True)


if __name__=='__main__':
    main()
