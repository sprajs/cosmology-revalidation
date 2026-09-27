"""Separate dense-metric retry of the unchanged synthetic optical target."""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor
import datetime
import json
import multiprocessing
import os
from pathlib import Path
import shutil
import sys
import time

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'bayesn_heldout'))
import model
import synthetic
DESIGN=HERE/'bayesn-synthetic-retune-design.json'
ROOT=model.ROOT


def dependencies():
    return {str(path.relative_to(ROOT)):model.sha(path) for path in
            [Path(__file__).resolve(),DESIGN,HERE/'bayesn_synthetic_retune_validate.py']}


def settings():
    design=json.loads(DESIGN.read_text())['sampler']
    return ({key:design[key] for key in ('target_accept_prob','dense_mass','max_tree_depth')},
            {key:design[key] for key in ('num_warmup','num_samples')}|
            {'num_chains':design['num_chains_per_invocation'],'progress_bar':False})


def target(kernel,data,arm):
    """The scientific density is the original callable, without modifications."""
    return model.numpyro_model(kernel,data,arm)


def make_sampler(density, initial):
    from numpyro.infer import MCMC, NUTS, init_to_value
    nuts_settings, run_settings = settings()
    return MCMC(NUTS(density, init_strategy=init_to_value(values=initial), **nuts_settings),
                **run_settings)


def prepare(original,work,validation):
    assert not work.exists(),'Use a fresh attempt directory.'
    checked=json.loads(validation.read_text())
    assert checked['status']=='passed_same_density_and_sampler_settings'
    assert checked['new_source_sha256']==dependencies()
    prior=json.loads((original/'prepared.json').read_text())
    assert prior['source_sha256']==model.source_hashes()
    design=json.loads(DESIGN.read_text())
    assert model.sha(model.HERE/'design.json')==design['original_design_sha256']
    gradient=ROOT/'studies/infrared/results/bayesn-synthetic-gradients.json'
    assert model.sha(gradient)==design['gradient_prerequisite_sha256']
    assert json.loads(gradient.read_text())['native_potential_gradients_bitwise_equal']
    model.check_parent(Path(prior['archive']))
    assert model.sha(original/'case-00/optical.json')==prior['cases'][0]['optical_sha256']
    assert model.sha(original/'release-filter-config.yaml')==prior['filter_config_sha256']
    (work/'case-00').mkdir(parents=True)
    for name in ('prepared.json','release-filter-config.yaml','case-00/optical.json'):
        shutil.copyfile(original/name,work/name)
        assert model.sha(original/name)==model.sha(work/name)
    plan=dict(original_work=str(original),source_sha256=dependencies(),
              original_prepared_sha256=model.sha(original/'prepared.json'),
              validation_path=str(validation),validation_sha256=model.sha(validation),
              original_target_source_sha256=prior['source_sha256'],
              original_optical_sha256=prior['cases'][0]['optical_sha256'],
              original_gradient_prerequisite_sha256=model.sha(gradient),settings=settings(),
              data_scope='Only original optical payload and generation metadata copied. No truth or NIR payload copied.',
              prepared_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    model.dump(work/'retuning-plan.json',plan)
    return plan


def fit_arm(spec):
    model.configure_cpu(spec['cpu'])
    import numpy as np
    import jax
    import numpyro
    assert dependencies()==spec['retuning_source_sha256']
    assert model.source_hashes()==spec['source_sha256']
    work=Path(spec['work']);folder=work/'case-00'/spec['arm'];folder.mkdir(exist_ok=False)
    optical=work/'case-00/optical.json'
    assert model.sha(optical)==spec['optical_sha256']
    data=json.loads(optical.read_text());assert data['scope']=='synthetic'
    assert all(row['role']=='optical' for row in data['rows'])
    config=work/'release-filter-config.yaml';assert model.sha(config)==spec['filter_config_sha256']
    def guard(event,arguments):
        if event=='open':
            name=str(arguments[0]).lower()
            if 'synthetic-nir.json' in name or 'nir-sealed' in name or '/photometry/' in name or name.endswith('truth.json'):
                raise RuntimeError('Optical-only retry cannot read held-out/truth/source-photometry files.')
    sys.addaudithook(guard)
    model.dump(folder/'worker-start.json',spec|{'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
    start=time.monotonic();records=[]
    design=json.loads(DESIGN.read_text())
    original_design=json.loads((model.HERE/'design.json').read_text())
    try:
        kernel=model.Kernel(Path(spec['archive']),config,data['metadata'],data['rows'],original_design['integration_bins'])
        density=target(kernel,data,spec['arm']);nuts_settings,run_settings=settings()
        for index in range(design['sampler']['chains_per_arm']):
            assert dependencies()==spec['retuning_source_sha256'] and model.source_hashes()==spec['source_sha256']
            seed=design['sampler']['chain_seeds'][index]+design['sampler']['arm_seed_offset'][spec['arm']]
            initial=dict(D=data['metadata']['mu_LCDM'],AV=.2,RV=3.1,theta=0.,epsilon_white=np.zeros(42),
                         tau=original_design['sampler']['initial_tau'][index])
            sampler=make_sampler(density,initial)
            tic=time.monotonic()
            sampler.run(jax.random.PRNGKey(seed),extra_fields=('energy','potential_energy','num_steps','diverging','accept_prob'))
            samples={key:np.asarray(value) for key,value in sampler.get_samples().items()}
            fields={key:np.asarray(value) for key,value in sampler.get_extra_fields().items()}
            adapted=sampler.last_state.adapt_state
            matrices={};metrics=[]
            for number,(names,value) in enumerate(adapted.inverse_mass_matrix.items()):
                matrix=np.asarray(value);key=f'adapt_inverse_mass_{number}';matrices[key]=matrix
                assert matrix.shape==(47,47),'The declared dense metric must span all47unconstrained dimensions.'
                metrics.append(dict(latent_site_order=list(names),array=key,shape=list(matrix.shape)))
            path=folder/f'chain-{index}.npz';assert not path.exists()
            np.savez_compressed(path,**samples,**{f'stat_{key}':value for key,value in fields.items()},**matrices,
                                adapt_step_size=np.asarray(adapted.step_size))
            record=dict(chain=index,seed=seed,warmup=run_settings['num_warmup'],draws=run_settings['num_samples'],
                        seconds=time.monotonic()-tic,archive=path.name,sha256=model.sha(path),initial_tau=initial['tau'],
                        divergences=int(fields['diverging'].sum()),max_num_steps=int(fields['num_steps'].max()),
                        final_adapted_step_size=float(adapted.step_size),dense_metric=metrics)
            model.dump(folder/f'chain-{index}.json',record);records.append(record)
            print(json.dumps({'arm':spec['arm'],**record}),flush=True)
            assert dependencies()==spec['retuning_source_sha256'] and model.source_hashes()==spec['source_sha256']
            if time.monotonic()-start>design['resources']['max_seconds_per_arm']:
                raise TimeoutError('Declared finite retry cap exceeded; preserve all completed chains.')
        assert model.sha(optical)==spec['optical_sha256']
        result=dict(scope='Separate sampler retry; unchanged synthetic optical density, no NIR access.',worker_spec=spec,
                    chains=records,seconds=time.monotonic()-start,
                    versions=dict(numpy=np.__version__,jax=jax.__version__,numpyro=numpyro.__version__))
        model.dump(folder/'posterior-manifest.json',result)
        return result
    except BaseException as error:
        model.dump(folder/'worker-failure.json',dict(error=repr(error),completed_chains=records,
                   seconds=time.monotonic()-start,source_sha256=dependencies()))
        raise


def run(work):
    plan=json.loads((work/'retuning-plan.json').read_text())
    assert plan['source_sha256']==dependencies()
    assert model.sha(plan['validation_path'])==plan['validation_sha256']
    prior=json.loads((work/'prepared.json').read_text())
    assert model.sha(work/'prepared.json')==plan['original_prepared_sha256']
    assert prior['source_sha256']==model.source_hashes()
    model.check_parent(Path(prior['archive']))
    cpus=sorted(os.sched_getaffinity(0));assert len(cpus)>=2
    specs=[]
    for index,arm in enumerate(('LCDM','broad')):
        spec=synthetic.worker_spec(work,0,arm,cpus[index],prior)
        spec['retuning_source_sha256']=dependencies();specs.append(spec)
    model.dump(work/'retuning-execution.json',dict(source_sha256=dependencies(),worker_specs=specs,settings=settings()))
    with ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context('spawn')) as pool:
        list(pool.map(fit_arm,specs))
    model.check_parent(Path(prior['archive']))
    assert dependencies()==plan['source_sha256']
    for arm in ('LCDM','broad'):
        result=synthetic.diagnose(work,0,arm)
        print(json.dumps({key:value for key,value in result.items() if key!='scalar_parameters'}),flush=True)


def prepare_score(work):
    plan=json.loads((work/'retuning-plan.json').read_text())
    assert plan['source_sha256']==dependencies()
    original=Path(plan['original_work']);prior=json.loads((work/'prepared.json').read_text())
    frozen={}
    for arm in ('LCDM','broad'):
        folder=work/'case-00'/arm;diagnostic=json.loads((folder/'diagnostics.json').read_text())
        assert diagnostic['sampler_passed']
        assert diagnostic['posterior_manifest_sha256']==model.sha(folder/'posterior-manifest.json')
        manifest=json.loads((folder/'posterior-manifest.json').read_text())
        assert manifest['worker_spec']['retuning_source_sha256']==dependencies()
        for record in manifest['chains']:
            path=folder/record['archive'];assert model.sha(path)==record['sha256'];frozen[str(path)]=record['sha256']
        frozen[str(folder/'diagnostics.json')]=model.sha(folder/'diagnostics.json')
    model.dump(work/'before-synthetic-nir-copy.json',dict(optical_sha256=frozen,source_sha256=dependencies()))
    source=original/'case-00/synthetic-nir.json';destination=work/'case-00/synthetic-nir.json'
    assert not destination.exists() and model.sha(source)==prior['cases'][0]['nir_sha256']
    shutil.copyfile(source,destination);assert model.sha(destination)==prior['cases'][0]['nir_sha256']


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('prepare','run','prepare-score'))
    parser.add_argument('--work',type=Path,required=True)
    parser.add_argument('--original',type=Path)
    parser.add_argument('--validation',type=Path)
    args=parser.parse_args();work=args.work.resolve()
    if args.action=='prepare':
        assert args.original and args.validation
        prepare(args.original.resolve(),work,args.validation.resolve())
    elif args.action=='run':run(work)
    else:prepare_score(work)
