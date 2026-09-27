"""Bounded native refinement with profiled nuisance parameters and explicit gates.

Original candidate records and active target implementations are read-only.
The correction changes the numerical search objective, never the likelihood.
"""
import os
for _key in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ.setdefault(_key, '1')
os.environ.setdefault('CLIPY_NOJAX', '1')
from pathlib import Path
import datetime
import hashlib
import json
import time
import numpy as np
from scipy.optimize import minimize
from cobaya.model import get_model
from modern_fast import configuration, identify
from modern_profile import evaluate, sha, SURROGATE
from late_geometry import sample_path

ROOT = Path(__file__).resolve().parents[4]
WORK = ROOT / '.work/unified-cosmology/inference/modern-profile-refine'
PRIOR = ROOT / 'studies/unified_cosmology/results/inference/modern-profile.json'
RESULT = ROOT / 'studies/unified_cosmology/results/inference/modern-profile-refinement.json'
CAP = 80
NUIS = ['A_planck', 'P_act', 'Tcal', 'Ecal', 'A_fg']
COSMO = ['H0', 'ombh2', 'omch2', 'logA', 'ns', 'tau']
H = .08
GTOL = .10
NOISE_TOL = .05


class BudgetReached(RuntimeError):
    pass


class Ledger:
    def __init__(self):
        self.native_calls = 0
        self.events = []

    def charge(self, label, parameters):
        if self.native_calls >= CAP:
            raise BudgetReached('80 native spectrum calculations exhausted')
        self.native_calls += 1
        event = {'number': self.native_calls, 'model': label, 'parameters': parameters}
        self.events.append(event)
        print(json.dumps({'native_spectrum': event}), flush=True)


def solve(label, prior, ledger):
    started = time.monotonic()
    prior_fit = prior['proposal'][label]
    baseline = prior['native'][label]['best_evaluated_candidate']
    c_names = COSMO + (['w', 'wa'] if label == 'cpl' else [])
    names = c_names + NUIS
    covfile = ROOT / prior_fit['author_covariance_file']
    columns = covfile.open().readline().lstrip('#').split()
    raw = np.loadtxt(covfile)
    ix = [columns.index(n) for n in names]
    covariance = raw[np.ix_(ix, ix)]
    centre = np.array([baseline['point'][n] for n in names])
    full_chol = np.linalg.cholesky(covariance)
    d = len(c_names)
    c_chol = np.linalg.cholesky(covariance[:d, :d])
    regression = np.linalg.solve(covariance[:d, :d], covariance[:d, d:]).T
    n_chol = np.linalg.cholesky(covariance[d:, d:] - regression @ covariance[:d, d:])
    gaussian = prior_fit['gaussian_calibration_priors']
    proposal_info = configuration(model=label, evolution='none', sample='dovekie', calibration='official_planck', surrogate=SURROGATE)
    native_info = configuration(model=label, evolution='none', sample='dovekie', calibration='official_planck')
    proposal_info['debug'] = native_info['debug'] = False
    identities = {'proposal': identify(proposal_info, sample_path('dovekie'), SURROGATE), 'native': identify(native_info, sample_path('dovekie'))}
    records = open(WORK / f'{label}-evaluations.jsonl', 'w')
    calls = {'proposal': 0, 'native': 0}

    with get_model(proposal_info) as proposal, get_model(native_info) as native:
        original_calculate = native.theory['camb'].calculate

        def counted_calculate(state, want_derived=True, **parameters):
            ledger.charge(label, {k: float(v) for k, v in parameters.items()})
            return original_calculate(state, want_derived=want_derived, **parameters)

        native.theory['camb'].calculate = counted_calculate

        def score(which, values, phase):
            model = proposal if which == 'proposal' else native
            calls[which] += 1
            before = proposal.theory['spectral_surrogate'].exact_calls
            result = evaluate(model, dict(zip(names, map(float, values))), gaussian)
            if proposal.theory['spectral_surrogate'].exact_calls != before:
                ledger.charge(label + '_unexpected_surrogate_fallback', {})
                raise BudgetReached('Surrogate exact fallback: stop rather than consume unplanned native points')
            records.write(json.dumps({'phase': phase, 'theory': which, 'result': result, 'input': values.tolist()}) + '\n')
            records.flush()
            return result

        cache = {}

        def full_score(x):
            key = np.asarray(x, dtype='<f8').tobytes()
            if key not in cache:
                r = score('proposal', centre + full_chol @ x, 'finish_surrogate')
                cache[key] = r
            return cache[key]['penalized_objective'] if cache[key] else 1e20

        initial = np.zeros(len(names))
        complete = minimize(full_score, initial, method='L-BFGS-B', bounds=[(-3., 3.)] * len(names),
                            options={'maxfun': 8000, 'maxiter': 300, 'eps': .002, 'ftol': 1e-11, 'gtol': .02, 'maxls': 25})
        surrogate_best = min((r for r in cache.values() if r), key=lambda r: r['penalized_objective'])
        anchor = np.array([surrogate_best['point'][n] for n in names])
        print(json.dumps({'model': label, 'surrogate_finished': float(surrogate_best['penalized_objective']), 'success': bool(complete.success), 'message': str(complete.message), 'calls': complete.nfev}), flush=True)
        # Coordinates are centred on the completed surrogate candidate. Their
        # cosmological scale is the marginal author covariance, not a new prior.
        profiles = {'proposal': {}, 'native': {}}

        def profile(which, x, phase):
            x = np.asarray(x)
            key = x.astype('<f8').tobytes()
            if key in profiles[which]:
                return profiles[which][key]
            cosmology = anchor[:d] + c_chol @ x
            n_base = anchor[d:] + regression @ (cosmology - anchor[:d])
            local = {}

            def function(y):
                k = np.asarray(y, dtype='<f8').tobytes()
                if k not in local:
                    values = np.r_[cosmology, n_base + n_chol @ y]
                    local[k] = score(which, values, phase)
                return local[k]['penalized_objective'] if local[k] else 1e20

            def gradient(y):
                step = .02
                return np.array([(function(y + np.eye(5)[j] * step) - function(y - np.eye(5)[j] * step)) / (2 * step) for j in range(5)])

            optimum = minimize(function, np.zeros(5), jac=gradient, method='BFGS', options={'maxiter': 40, 'gtol': .01})
            # Requiring the measured gradient avoids trusting an optimizer's
            # success string when the component likelihood has finite precision.
            y = optimum.x
            grad = gradient(y)
            best = local[y.astype('<f8').tobytes()] if y.astype('<f8').tobytes() in local else None
            if best is None:
                function(y)
                best = local[y.astype('<f8').tobytes()]
            if best is None:
                raise RuntimeError('No finite fixed-spectrum nuisance optimum')
            out = {'coordinates': x.tolist(), 'record': best, 'nuisance_coordinates': y.tolist(), 'nuisance_gradient': grad.tolist(),
                   'nuisance_stationary': bool(np.max(abs(grad)) <= .05), 'optimizer_success': bool(optimum.success), 'optimizer_message': str(optimum.message), 'optimizer_evaluations': int(optimum.nfev)}
            profiles[which][key] = out
            if which == 'native':
                print(json.dumps({'model': label, 'native_profile': len(profiles[which]), 'phase': phase, 'S': best['penalized_objective'], 'nuisance_gradient_max': float(np.max(abs(grad)))}), flush=True)
            return out

        def value(which, x, phase):
            return profile(which, x, phase)['record']['penalized_objective']

        zero = np.zeros(d)
        native0 = value('native', zero, 'initial_correction_centre')
        proposal0 = value('proposal', zero, 'initial_correction_centre')
        delta0 = native0 - proposal0
        dg = np.array([(value('native', H * np.eye(d)[j], 'initial_forward_correction') - value('proposal', H * np.eye(d)[j], 'initial_forward_correction') - delta0) / H for j in range(d)])
        iterations = []
        current = zero

        for iteration, radius in enumerate([.35, .25]):
            # Linear exact-minus-surrogate correction in a bounded local box.
            # Optimize nuisances jointly here; nesting a complete nuisance
            # optimizer inside every outer finite difference is unnecessary.
            start_nuisance = np.array(profile('proposal', current, 'trust_search_start')['nuisance_coordinates'])
            joint_cache = {}

            def corrected(joint):
                x, y = joint[:d], joint[d:]
                key = np.asarray(joint, dtype='<f8').tobytes()
                if key not in joint_cache:
                    cosmology = anchor[:d] + c_chol @ x
                    nuisances = anchor[d:] + regression @ (cosmology - anchor[:d]) + n_chol @ y
                    joint_cache[key] = score('proposal', np.r_[cosmology, nuisances], f'corrected_search_{iteration}')
                record = joint_cache[key]
                return record['penalized_objective'] + float(dg @ (x - current)) if record else 1e20

            fit = minimize(corrected, np.r_[current, start_nuisance], method='L-BFGS-B',
                           bounds=[(v - radius, v + radius) for v in current] + [(-4., 4.)] * 5,
                           options={'maxfun': 3000, 'maxiter': 150, 'eps': .002, 'ftol': 1e-10, 'gtol': .02, 'maxls': 25})
            next_x = np.asarray(fit.x[:d])
            point = profile('native', next_x, f'trust_candidate_{iteration}')
            native_gradient = np.empty(d)
            delta_gradient = np.empty(d)
            curvature = np.empty(d)
            for j in range(d):
                step = H * np.eye(d)[j]
                np_ = value('native', next_x + step, f'trust_gradient_{iteration}')
                nm_ = value('native', next_x - step, f'trust_gradient_{iteration}')
                sp_ = value('proposal', next_x + step, f'trust_gradient_{iteration}')
                sm_ = value('proposal', next_x - step, f'trust_gradient_{iteration}')
                native_gradient[j] = (np_ - nm_) / (2 * H)
                delta_gradient[j] = ((np_ - sp_) - (nm_ - sm_)) / (2 * H)
                curvature[j] = (np_ + nm_ - 2 * point['record']['penalized_objective']) / H**2
            trial = {'iteration': iteration, 'radius': radius, 'centre': current.tolist(), 'candidate': point, 'correction_gradient': dg.tolist(), 'native_gradient': native_gradient.tolist(), 'native_diagonal_curvature': curvature.tolist(), 'trust_boundary': bool(np.max(abs(next_x - current)) >= radius - 1e-4), 'search_success': bool(fit.success), 'search_message': str(fit.message), 'search_evaluations': int(fit.nfev)}
            iterations.append(trial)
            print(json.dumps({'model': label, 'iteration_complete': iteration, 'S': point['record']['penalized_objective'], 'native_gradient_max': float(np.max(abs(native_gradient))), 'native_calls_total': ledger.native_calls}), flush=True)
            current = next_x
            dg = delta_gradient

        final = iterations[-1]
        grad = np.asarray(final['native_gradient'])
        worst = int(np.argmax(abs(grad)))
        half = H * np.eye(d)[worst] / 2
        half_gradient = (value('native', current + half, 'final_step_halving') - value('native', current - half, 'final_step_halving')) / H
        discrepancy = abs(half_gradient - grad[worst])
        all_native = list(profiles['native'].values())
        nuisance_gate = all(r['nuisance_stationary'] for r in all_native)
        gates = {'all_fixed_spectrum_nuisance_gradients_pass': nuisance_gate,
                 'final_native_gradient_max_below_0p1': bool(np.max(abs(grad)) <= GTOL),
                 'worst_direction_step_halving_difference_below_0p05': bool(discrepancy <= NOISE_TOL),
                 'final_candidate_not_on_trust_boundary': not final['trust_boundary']}
        best_native = min(all_native, key=lambda r: r['record']['penalized_objective'])
        output = {'identities': identities, 'cosmology_coordinates': c_names, 'nuisance_coordinates': NUIS,
                  'author_covariance': str(covfile.relative_to(ROOT)), 'author_covariance_sha256': sha(covfile),
                  'surrogate_completion': {'success': bool(complete.success), 'message': str(complete.message), 'nfev': int(complete.nfev), 'candidate': surrogate_best},
                  'anchor': dict(zip(names, anchor.tolist())), 'cosmology_cholesky': c_chol.tolist(), 'nuisance_conditional_cholesky': n_chol.tolist(), 'nuisance_regression': regression.tolist(),
                  'initial_native': profile('native', zero, 'cached_initial'), 'iterations': iterations,
                  'step_halving': {'direction_index': worst, 'direction': c_names[worst], 'step': H, 'half_step_gradient': float(half_gradient), 'full_step_gradient': float(grad[worst]), 'difference': float(discrepancy)},
                  'gates': gates, 'local_stationarity_gate_passed': all(gates.values()),
                  'stationarity_scope': 'Finite differences of the nuisance-profiled objective in the declared local scaled coordinates; half-step validation in the largest-gradient direction only. No global optimality certification.',
                  'best_native_tested_point': best_native, 'native_profiles': all_native, 'likelihood_calls': calls, 'seconds': time.monotonic() - started}
    records.close()
    (WORK / f'{label}-complete.json').write_text(json.dumps(output, indent=2) + '\n')
    return output


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    prior = json.loads(PRIOR.read_text())
    frozen = dict(prior['source_sha256'])
    frozen[str(PRIOR.relative_to(ROOT))] = sha(PRIOR)
    for p, h in frozen.items():
        assert sha(ROOT / p) == h, p
    design = {'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'native_spectrum_cap_total': CAP, 'workers': 1,
              'target': 'Unchanged modern_fast.configuration, official_planck, Dovekie, no SN luminosity evolution',
              'objective': 'Same S as modern-profile.json; exact spectra with numerically profiled five calibration/foreground nuisances.',
              'surrogate_finish_call_cap_per_model': 8000, 'initial_Taylor_difference': 'Forward derivative of profiled native-minus-surrogate objective at step0.08 in marginal-author-Cholesky cosmology coordinates',
              'trust_radii': [.35, .25], 'native_gradient_step': H, 'gradient_max_tolerance': GTOL,
              'step_halving_difference_tolerance': NOISE_TOL, 'nuisance_gradient_max_tolerance': .05,
              'planned_native_points': {'lcdm': 35, 'cpl': 45},
              'interpretation': 'Local finite-difference stationarity gates only; no global-optimum, posterior, evidence or significance claim.'}
    (WORK / 'execution-design.json').write_text(json.dumps(design, indent=2) + '\n')
    ledger = Ledger()
    fits = {}
    failure = None
    try:
        for label in ['lcdm', 'cpl']:
            fits[label] = solve(label, prior, ledger)
    except Exception as error:
        failure = {'type': type(error).__name__, 'message': str(error)}
        print(json.dumps({'failure': failure}), flush=True)
    for p, h in frozen.items():
        assert sha(ROOT / p) == h, p
    comparison = None
    if len(fits) == 2:
        l = fits['lcdm']['iterations'][-1]['candidate']['record']
        c = fits['cpl']['iterations'][-1]['candidate']['record']
        comparison = {'final_local_candidates_LCDM_minus_CPL_S': l['penalized_objective'] - c['penalized_objective'], 'both_local_stationarity_gates_passed': all(f['local_stationarity_gate_passed'] for f in fits.values()), 'meaning': 'Native difference of final trust candidates; not a global profile improvement bound or significance.'}
    out = {'schema': 'modern-native-refinement-v1', 'design': design, 'status': 'completed' if failure is None else 'stopped_with_preserved_failure', 'failure': failure, 'models': fits,
           'native_spectrum_calculations': ledger.native_calls, 'native_events': ledger.events, 'comparison': comparison,
           'source_sha256': {**frozen, str(Path(__file__).relative_to(ROOT)): sha(__file__)}, 'surrogate_sha256': sha(SURROGATE),
           'work_records_sha256': {str(p.relative_to(ROOT)): sha(p) for p in sorted(WORK.glob('*')) if p.is_file()}}
    RESULT.write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'status': out['status'], 'comparison': comparison, 'native_spectrum_calculations': ledger.native_calls}), flush=True)


if __name__ == '__main__':
    main()
