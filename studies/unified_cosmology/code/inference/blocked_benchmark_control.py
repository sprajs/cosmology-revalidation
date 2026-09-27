"""Known-Gaussian blocking control with measured whole-block reference costs.

Separate from the capped real-target pilot: it remains executable if that pilot
needs too many native fallbacks. The cost model is not measured cosmology ESS.
"""
import argparse
import copy
import inspect
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import solve_triangular

from blocked_benchmark import (ROOT, HERE, DESIGN, COSMOLOGY, FAST, configuration,
    identify, sample_path, digest, proposal, kernel, metrics)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--surrogate', type=Path, required=True)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--pilot-results', type=Path, required=True,
                   help='Replay the preserved known-Gaussian draws, not real-target draws.')
    args = p.parse_args()
    args.work = args.work.resolve(); args.surrogate = args.surrogate.resolve()
    assert args.work.is_relative_to(ROOT/'.work') and not args.work.exists()
    args.work.mkdir(parents=True)
    import camb
    from cobaya.model import get_model
    from cobaya.samplers.mcmc.proposal import BlockedProposer
    original = camb.get_results
    def forbidden(*a, **k):
        raise RuntimeError('No native CMB spectrum calls allowed in reference/control stage.')
    camb.get_results = forbidden
    design = json.loads(DESIGN.read_text())
    report = {'status': 'in_progress', 'scope': 'Known Gaussian control; modeled ESS/CPU, not measured cosmology efficiency.',
              'design': design, 'whole_block_costs': {}, 'gaussian': [], 'native_spectrum_calls': 0}
    try:
        info = configuration('cpl', 'linear', 'dovekie', 'official_planck', args.surrogate)
        report['target_identity'] = identify(info, sample_path('dovekie'), args.surrogate)
        info['debug'] = 40
        with get_model(info) as model:
            names = list(model.parameterization.sampled_params())
            start, covariance, inputs = proposal(names, info)
            report['proposal'] = inputs
            theory = model.theory['spectral_surrogate']
            assert np.isfinite(model.logposterior(start, return_derived=False).logpost)
            for label, coordinates in [('slow', COSMOLOGY), ('fast', FAST)]:
                records = []
                for delta in np.linspace(-.05, .05, 12):
                    point = start.copy()
                    for i, name in enumerate(names):
                        if name in coordinates:
                            point[i] += delta*np.sqrt(covariance[i, i])
                    before = theory.surrogate_calls
                    begin = time.process_time()
                    density = model.logposterior(point, return_derived=False)
                    elapsed = time.process_time()-begin
                    assert np.isfinite(density.logpost)
                    added = theory.surrogate_calls-before
                    if label == 'fast':
                        assert added == 0
                    records.append({'CPU_seconds': elapsed, 'theory_calls': added})
                    model.logposterior(start, return_derived=False)
                report['whole_block_costs'][label] = records
            slowcost = float(np.median([r['CPU_seconds'] for r in report['whole_block_costs']['slow']]))
            fastcost = float(np.median([r['CPU_seconds'] for r in report['whole_block_costs']['fast']]))
        chol = np.linalg.cholesky(covariance)
        original_report = json.loads(args.pilot_results.read_text())
        assert original_report['status'] == 'completed_performance_only'
        report['pilot_results_sha256'] = digest(args.pilot_results)
        assert original_report['design'] == design
        for seed in design['gaussian_seeds']:
            for factor in design['oversampling_factors']:
                matches = [r for r in original_report['gaussian'] if r['seed'] == seed
                           and len(r['blocking']) == (1 if factor == 0 else 2)
                           and (factor == 0 or r['blocking'][1][0] == factor)]
                assert len(matches) == 1
                result = copy.deepcopy(matches[0])
                rawpath = ROOT/result['draws_path']
                assert digest(rawpath) == result['draws_sha256']
                with np.load(rawpath, allow_pickle=False) as data:
                    rows = data['rows']; assert data['names'].tolist() == names
                assert len(rows) == design['gaussian_attempts']
                modeled = result['slow_proposals']*slowcost+result['fast_proposals']*fastcost
                result['modeled_seconds'] = modeled
                result['cost_model'] = {'whole_slow_block_CPU_seconds': slowcost,
                                       'whole_fast_block_CPU_seconds': fastcost}
                result['metrics'] = metrics(rows, names, result['CPU_seconds'], modeled)
                white = solve_triangular(chol, rows.T, lower=True, check_finite=False).T
                result['whitened_coordinate_order'] = names
                result['whitened_metrics'] = metrics(white, names, result['CPU_seconds'], modeled)
                result['maximum_mean_error_in_known_marginal_SD'] = float(np.max(abs(rows[len(rows)//4:].mean(axis=0))/np.sqrt(np.diag(covariance))))
                report['gaussian'].append(result)
        paths = [Path(__file__), HERE/'blocked_benchmark.py', DESIGN, Path(inspect.getfile(BlockedProposer))]
        report['source_sha256'] = {str(path): digest(path) for path in paths}
        report['status'] = 'completed_known_Gaussian_control'
    finally:
        camb.get_results = original
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': report['status'], 'slow_seconds': slowcost, 'fast_seconds': fastcost,
                      'runs': len(report['gaussian'])}))


if __name__ == '__main__':
    main()
