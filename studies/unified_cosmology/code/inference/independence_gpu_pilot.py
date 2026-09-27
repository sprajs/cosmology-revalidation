"""Bounded initialized-GPU efficiency pilot; never a cosmological measurement."""
import argparse
import os
import json
from pathlib import Path
import time

import numpy as np
from scipy.special import logsumexp
from independence_proposal import ROOT, FrozenMixture, digest, relative
from independence_run import runtime
from independence_runtime import native_initialize, numerical_runtime

HERE = Path(__file__).resolve().parent


def run_pilot(proposal_folder, output, seed=273280, acceptance_seed=273281, calls=400):
    os.sched_setaffinity(0, [max(os.sched_getaffinity(0))])
    assert calls == 400, 'This bounded pilot design allows exactly400 requested points.'
    proposal_folder, output = Path(proposal_folder).resolve(), Path(output).resolve()
    assert not output.exists(); output.mkdir(parents=True)
    from independence_gpu_run import load_gpu_proposal,dependencies,VALIDATION
    driver_dependencies=dependencies()
    validation=json.loads(VALIDATION.read_text())
    assert validation['status']=='passed_configuration_only_no_likelihood_calls'
    assert validation['driver_dependency_sha256']==driver_dependencies
    info,identity,scientific,parent,mixture,sidecar=load_gpu_proposal(proposal_folder)
    assert scientific['evolution']=='linear', 'Only the allocated linear pilot is authorized.'
    training=json.loads((proposal_folder/'proposal.json').read_text())
    settings=parent['arguments']
    sources = [Path(__file__), *[HERE/name for name in ['independence_proposal.py', 'independence_runtime.py', 'independence_run.py','independence_gpu_run.py','independence-gpu-design.json']]]
    design = {'scope': 'GPU efficiency and support diagnostic only, no posterior measurement.',
              'GPU':True,'cpu_affinity':list(os.sched_getaffinity(0)),
              'GPU_driver_dependency_sha256':driver_dependencies,
              'GPU_configuration_validation_sha256':digest(VALIDATION),
              'GPU_proposal_sidecar_sha256':digest(proposal_folder/'gpu-proposal.json'),
              'model': settings['model'], 'evolution': settings['evolution'], 'requested_target_calls': calls,
              'one_additional_native_initialization': True, 'proposal_seed': seed, 'acceptance_seed': acceptance_seed,
              'proposal_record_sha256': digest(proposal_folder/'proposal.json'),
              'proposal_sha256': training['proposal_sha256'], 'scientific_target_identity': parent['target_identity']['identity'],
              'runtime': runtime(), 'initialization': 'First finite independently requested proposal; earlier failures retained and not counted as MH transitions.',
              'acceptance': 'min(1,exp(logpi(y)-logpi(x)+logq(x)-logq(y))); no learning/clipping/redraw',
              'scope_of_decision': 'Acceptance, holding concentration, candidate importance weights, short-chain ESS and target-call cost only; no production launch.',
              'source_sha256': {relative(p): digest(p) for p in sources}}
    (output/'design.json').write_text(json.dumps(design, indent=2)+'\n')
    rng = np.random.default_rng(seed)
    requested = [mixture.draw(rng) for _ in range(calls)]
    points = np.array([r[0] for r in requested]); components = [r[1] for r in requested]
    logqs = mixture.logpdf(points); logus = np.log(np.random.default_rng(acceptance_seed).random(calls))
    np.savez(output/'requests.npz', points=points, components=components, logq=logqs, log_uniform=logus, names=mixture.names)
    from measurement_summary import verify_current_target
    inputs = verify_current_target(parent)
    assert identity == parent['target_identity']
    initialization = native_initialize(info)
    (output/'native-initialization.json').write_text(json.dumps(initialization, indent=2)+'\n')
    from cobaya.model import get_model
    started = time.perf_counter(); current = None; logposts = []; occupied = []; records = []
    with get_model(info) as model:
        setup_seconds = time.perf_counter()-started
        assert list(model.parameterization.sampled_params()) == mixture.names
        theory = model.theory['spectral_surrogate']
        for index, point in enumerate(points):
            t = time.perf_counter(); exact_before = theory.exact_calls
            value = model.logposterior(dict(zip(mixture.names, map(float, point))))
            seconds = time.perf_counter()-t; logp = float(value.logpost); finite = bool(np.isfinite(logp))
            logposts.append(logp)
            if current is None:
                logalpha = None; accepted = finite
            else:
                logalpha = min(0., logp-logposts[current]+logqs[current]-logqs[index])
                accepted = bool(finite and logus[index] < logalpha)
            if accepted: current = index
            occupied.append(current)
            row = {'index': index, 'point': dict(zip(mixture.names, map(float, point))),
                   'component': components[index], 'logq': float(logqs[index]), 'target_logpost': logp if finite else None,
                   'finite': finite, 'prior_rejected': bool(value.logprior == -np.inf),
                   'log_uniform': float(logus[index]), 'log_acceptance': logalpha if logalpha is None or np.isfinite(logalpha) else None,
                   'accepted': accepted, 'occupied_candidate': current, 'seconds': seconds,
                   'successful_native_fallback_calls': theory.exact_calls-exact_before,
                   'component_logpriors': [float(v) if np.isfinite(v) else None for v in value.logpriors],
                   'component_loglikes': {k: float(v) if np.isfinite(v) else None for k, v in zip(model.likelihood, value.loglikes)}}
            records.append(row)
            with (output/'records.jsonl').open('a') as stream: stream.write(json.dumps(row, allow_nan=False)+'\n')
            if (index+1) % 25 == 0:
                print(json.dumps({'calls': index+1, 'finite': sum(r['finite'] for r in records),
                                  'accepted_including_initialization': sum(r['accepted'] for r in records),
                                  'successful_native_fallbacks': theory.exact_calls, 'elapsed': time.perf_counter()-started}), flush=True)
        final_runtime = numerical_runtime()
    assert dependencies()==driver_dependencies
    for path, expected in inputs.items(): assert digest(ROOT/path) == expected
    for path, expected in design['source_sha256'].items(): assert digest(ROOT/path) == expected
    finite = np.isfinite(logposts)
    lw = np.asarray(logposts)-logqs
    weights = np.exp(lw-logsumexp(lw)) if finite.any() else np.zeros(calls)
    active = np.array([i for i in occupied if i is not None]); _, counts = np.unique(active, return_counts=True)
    first = next((i for i,r in enumerate(records) if r['accepted']), None)
    ess = {}
    if len(active) >= 20:
        import arviz as az
        ess = {name: {'identity': float(az.ess(points[active, k], method='identity')),
                      'bulk': float(az.ess(points[active, k], method='bulk'))}
               for k, name in enumerate(mixture.names)}
    result = {'status': 'completed_efficiency_pilot_no_measurement', 'target_identity': parent['target_identity']['identity'],
              'model': settings['model'], 'evolution': settings['evolution'], 'target_calls': calls,
              'finite_candidates': int(finite.sum()), 'prior_rejections': sum(r['prior_rejected'] for r in records),
              'other_invalid_target_candidates': sum(not r['finite'] and not r['prior_rejected'] for r in records),
              'native_initialization_calls': 1, 'native_counter_semantics': 'Frozen exact_calls counter counts successful native fallbacks; all invalid target attempts remain separately recorded.', 'successful_native_fallback_calls': sum(r['successful_native_fallback_calls'] for r in records),
              'accepted_including_initialization': sum(r['accepted'] for r in records),
              'initialization_attempts': first+1 if first is not None else calls,
              'acceptance_excluding_initialization': (sum(r['accepted'] for r in records)-1)/(len(active)-1) if len(active)>1 else None,
              'unique_occupied_states': len(counts), 'longest_hold': int(counts.max()) if len(counts) else None,
              'largest_occupancy_fraction': float(counts.max()/len(active)) if len(active) else None,
              'occupancy_effective_states': float(len(active)**2/(counts@counts)) if len(active) else None,
              'iid_proposal_raw_importance_ESS': float(1/(weights@weights)) if finite.any() else None,
              'largest_candidate_normalized_weight': float(weights.max()) if finite.any() else None,
              'coordinate_short_single_chain_ESS_not_qualification': ess,
              'setup_seconds': setup_seconds, 'evaluation_seconds': sum(r['seconds'] for r in records),
              'evaluation_time_quantiles': np.quantile([r['seconds'] for r in records], [0,.5,.9,.99,1]).tolist(),
              'runtime_after_model': final_runtime,
              'source_sha256': design['source_sha256'], 'scientific_inputs_sha256': inputs,
              'artifact_sha256': {relative(p): digest(p) for p in output.iterdir() if p.is_file()},
              'qualification': 'Frozen-proposal efficiency only. Finite single-chain diagnostics cannot establish posterior convergence or undiscovered mode coverage.'}
    (output/'result.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['status', 'target_calls', 'finite_candidates', 'successful_native_fallback_calls', 'acceptance_excluding_initialization', 'iid_proposal_raw_importance_ESS', 'evaluation_seconds']}, indent=2))
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('proposal_folder', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--seed', type=int, default=273280); p.add_argument('--acceptance-seed', type=int, default=273281)
    a = p.parse_args(); run_pilot(a.proposal_folder, a.output, a.seed, a.acceptance_seed)
