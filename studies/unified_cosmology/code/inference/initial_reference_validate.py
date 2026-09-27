"""Independent starting-point support checks; never evaluate a likelihood."""
import copy
import hashlib
import json
from pathlib import Path

import numpy as np

from modern_sample import initial_reference
from modern_fast import configuration

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def constraints(info, point):
    # Separately written scalar/support predicate, without clipping.
    failed = []
    for name, value in point.items():
        prior = info['params'][name]['prior']
        if not np.isfinite(value):
            failed.append(name)
        elif 'min' in prior and not value > prior['min']:
            failed.append(name)
        elif 'max' in prior and not value < prior['max']:
            failed.append(name)
    if point.get('w', -1.) + point.get('wa', 0.) > 0.:
        failed.append('CAMB_w_plus_wa')
    return failed


def main():
    import camb
    def forbidden(*args, **kwargs):
        raise AssertionError('Initialization validation must not evaluate CAMB.')
    camb.get_results = forbidden
    camb.get_background = forbidden
    author_path = ROOT/'studies/unified_cosmology/results/external_probes/author-configuration.json'
    author = json.loads(author_path.read_text())
    sources = [Path(__file__), HERE/'modern_sample.py', HERE/'modern_fast.py',
               HERE/'modern_run.py', author_path]
    rows = []
    for model in ['lcdm', 'cpl']:
        info = configuration(model=model)
        unchanged = copy.deepcopy(info)
        proposal = author['models'][model]['proposal_only']
        covariance_path = ROOT/proposal['path']
        assert digest(covariance_path) == proposal['sha256']
        sources.append(covariance_path)
        names = covariance_path.open().readline().lstrip('#').split()
        covariance = np.loadtxt(covariance_path)
        mean = proposal['transformed_mean']
        for seed_base in ([272652, 272672] if model == 'lcdm' else [272652]):
            for rank in range(4):
                seed = seed_base + rank
                result = initial_reference(info, names, mean, covariance, seed)
                rng = np.random.default_rng(seed)
                first = None
                rejected = []
                for attempt in range(10000):
                    vector = rng.multivariate_normal(np.zeros(len(names)), covariance)
                    point = {key: float(mean[key] + vector[i]) for i, key in enumerate(names)
                             if isinstance(info['params'].get(key), dict)
                             and 'prior' in info['params'][key]}
                    if first is None:
                        first = dict(point)
                    failed = constraints(info, point)
                    if not failed:
                        break
                    rejected.append(failed)
                else:
                    raise AssertionError('Independent reconstruction exhausted support attempts.')
                assert result['accepted_attempt'] == attempt + 1
                assert result['point'] == point
                assert result['rejected_support_constraints'] == rejected
                assert constraints(info, result['point']) == []
                if not constraints(info, first):
                    assert result['accepted_attempt'] == 1 and result['point'] == first
                rows.append({'model': model, 'base_seed': seed_base, 'rank': rank,
                             'first_draw': first, 'first_draw_constraints': constraints(info, first),
                             'accepted_attempt': result['accepted_attempt'],
                             'accepted_point': result['point'],
                             'all_draw_components_preserved_without_clipping': True})
        assert info == unchanged, 'Initialization mutated scientific configuration.'
    historical = [r for r in rows if r['model'] == 'lcdm' and r['base_seed'] == 272652 and r['rank'] == 2][0]
    assert historical['first_draw_constraints'] == ['A_fg']
    assert historical['accepted_attempt'] > 1
    # A zero-variance vector permanently outside the box must hit the explicit
    # finite retry bound; an endpoint itself is also not silently nudged inward.
    failure_cases = []
    original_rng = np.random.default_rng
    class Counter:
        def __init__(self):
            self.calls = 0
            self.generator = original_rng(8)
        def multivariate_normal(self, *args, **kwargs):
            self.calls += 1
            return self.generator.multivariate_normal(*args, **kwargs)
    for mean in [1., 2., float('nan')]:
        counter = Counter()
        np.random.default_rng = lambda seed: counter
        try:
            initial_reference({'params': {'x': {'prior': {'min': 0., 'max': 1.}}}},
                              ['x'], {'x': mean}, np.zeros((1, 1)), 8)
        except RuntimeError as error:
            assert str(error) == 'Could not draw an initial point inside the declared support.'
            assert counter.calls == 10000
            failure_cases.append({'case': 'nonfinite' if not np.isfinite(mean) else 'endpoint' if mean == 1. else 'outside_box',
                                  'attempts': counter.calls})
        else:
            raise AssertionError('Impossible starting support did not terminate with the declared error.')
        finally:
            np.random.default_rng = original_rng
    counter = Counter()
    np.random.default_rng = lambda seed: counter
    try:
        initial_reference({'params': {k: {'prior': {'min': -2., 'max': 2.}} for k in ['w', 'wa']}},
                          ['w', 'wa'], {'w': -.4, 'wa': .5}, np.zeros((2, 2)), 8)
    except RuntimeError:
        assert counter.calls == 10000
        failure_cases.append({'case': 'CAMB_w_plus_wa', 'attempts': counter.calls})
    else:
        raise AssertionError('Invalid coupled CAMB support was not rejected.')
    finally:
        np.random.default_rng = original_rng
    report = {'status': 'passed', 'native_theory_calls': 0,
              'original_seed_cases': 8, 'replacement_lcdm_seed_cases': 4,
              'rows': rows, 'unchanged_scientific_configurations': True,
              'retry_cap_failure_cases': failure_cases,
              'source_sha256': {str(p.relative_to(ROOT)): digest(p) for p in sources},
              'scope': 'Initialization-only numerical support checks. No target, prior, likelihood or active campaign was modified.'}
    path = ROOT/'studies/unified_cosmology/results/inference/initial-reference-validation.json'
    path.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ['rows', 'source_sha256']}, indent=2))


if __name__ == '__main__':
    main()
