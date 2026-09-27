"""Verify the separate retry's unchanged density and actual sampler settings.

This uses only a declared linear synthetic kernel. No SED, sampling, source
photometry, held-out vector, or truth values are evaluated.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import bayesn_synthetic_retune as retry


def validate(original):
    retry.model.configure_cpu()
    import jax
    import jax.numpy as jnp
    import numpy as np
    from numpyro.infer.util import log_density
    design = json.loads(retry.DESIGN.read_text())
    original_design = json.loads((retry.model.HERE/'design.json').read_text())
    prepared = json.loads((original/'prepared.json').read_text())
    assert prepared['source_sha256'] == retry.model.source_hashes()
    assert retry.model.sha(retry.model.HERE/'design.json') == design['original_design_sha256']
    assert retry.synthetic.numpyro_model is retry.model.numpyro_model
    assert retry.synthetic.Kernel is retry.model.Kernel
    assert len(design['sampler']['chain_seeds']) == design['sampler']['chains_per_arm'] == 4
    old_seeds = {seed+100*case+offset for seed in original_design['sampler']['chain_seeds']
                 for case in range(12) for offset in (0,10)}
    new_seeds = {seed+offset for seed in design['sampler']['chain_seeds'] for offset in (0,10)}
    assert len(new_seeds) == 8 and new_seeds.isdisjoint(old_seeds)
    optical=original/'case-00/optical.json'
    assert retry.model.sha(optical) == prepared['cases'][0]['optical_sha256']
    assert retry.model.sha(original/'release-filter-config.yaml') == prepared['filter_config_sha256']
    assert design['resources'] == dict(CPU_workers=2,max_seconds_per_arm=7200,
                                      external_process_group_limit_seconds=7200,termination_grace_seconds=30)
    rng=np.random.default_rng(273499)
    matrix=rng.normal(size=(3,47))*.02
    class Linear:
        def forward(self, parameters):
            return jnp.asarray(matrix)@parameters
    data=dict(metadata=dict(mu_LCDM=35.,sigma_external=.1),errors=[.4,.6,.8],flux=[1.,-.2,.3])
    cases=[];actual_settings=[]
    for arm in ('LCDM','broad'):
        density=retry.target(Linear(),data,arm)
        reference=retry.model.numpyro_model(Linear(),data,arm)
        distances=([34.6,34.9,35.,35.1,35.4,35.7] if arm=='LCDM'
                   else [19.8,20.05,34.8,35.,49.95,50.2])
        for index,distance in enumerate(distances):
            point=dict(D=distance,AV=.05+.12*index,RV=1.3+.7*index,theta=-.5+.2*index,
                       epsilon_white=jnp.asarray(rng.normal(size=42)),tau=-8+5*index)
            a,ta=log_density(density,(),{},point);b,tb=log_density(reference,(),{},point)
            assert float(a)==float(b)
            assert set(ta)==set(tb)
            for key in ta:
                assert np.array_equal(ta[key]['value'],tb[key]['value'])
            cases.append(dict(arm=arm,D=distance,absolute_log_density_difference=float(abs(a-b))))
        initial=dict(D=35.,AV=.2,RV=3.1,theta=0.,epsilon_white=np.zeros(42),tau=-5.)
        sampler=retry.make_sampler(density,initial)
        kernel=sampler.sampler
        assert kernel._target_accept_prob == .99 and kernel._dense_mass is True
        assert kernel._max_tree_depth == 10
        assert sampler.num_warmup==1000 and sampler.num_samples==1000 and sampler.num_chains==1
        # Initialize the synthetic linear kernel only: verifies native dense
        # metric layout and all latent sites without advancing a chain.
        state=kernel.init(jax.random.PRNGKey(273499),num_warmup=1000,model_args=(),model_kwargs={})
        metric=state.adapt_state.inverse_mass_matrix
        assert len(metric)==1
        names,value=next(iter(metric.items()))
        assert np.asarray(value).shape==(47,47)
        assert set(names)=={'D','AV','RV','theta','epsilon_white','tau'}
        assert np.isfinite(value).all() and float(state.adapt_state.step_size)>0
        actual_settings.append(dict(arm=arm,target_accept_prob=kernel._target_accept_prob,
                                    dense_mass=kernel._dense_mass,max_tree_depth=kernel._max_tree_depth,
                                    warmup=sampler.num_warmup,draws=sampler.num_samples,chains=sampler.num_chains,
                                    metric_shape=list(value.shape),latent_site_order=list(names)))
    return dict(status='passed_same_density_and_sampler_settings',new_source_sha256=retry.dependencies(),
                original_source_sha256=retry.model.source_hashes(),case_checks=cases,actual_sampler_settings=actual_settings,
                density_max_absolute_difference=max(x['absolute_log_density_difference'] for x in cases),
                source_identity='Original Kernel and numpyro_model are imported unchanged; actual MCMC factory inspected.',
                input_sha256={str(p.relative_to(retry.ROOT)):retry.model.sha(p) for p in
                              (original/'prepared.json',optical,original/'release-filter-config.yaml')},
                original_seeds=sorted(old_seeds),new_seeds=sorted(new_seeds),
                resources=design['resources'],SED_calls=0,sampling_steps=0,
                interpretation='Engineering identity/settings validation, not a test of posterior convergence or physical model adequacy.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--original',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=validate(args.original.resolve())
    retry.model.dump(args.output.resolve(),result)
    print(json.dumps(result,indent=2,allow_nan=False))
