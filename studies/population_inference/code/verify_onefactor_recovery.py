"""Independent arithmetic and direct latent-integral review of synthetic recovery.

Imports no execution or production-kernel code. This is a post-outcome review,
not a newly registered scientific estimator experiment.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.special import log_ndtr

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'runs/research_2026_09_26/onefactor_recovery'
OUT = BASE / 'root-independent-review'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def write(p, value):
    p.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def main():
    OUT.mkdir(exist_ok=False)
    manifest = json.loads((BASE / 'manifest.json').read_text())
    for p, digest in manifest['files_sha256'].items():
        assert sha(ROOT / p) == digest, p
    z = np.load(BASE / 'draws.npz')
    y, h, d, v = [z[k] for k in ['flux', 'mean', 'diagonal_variance', 'loading']]
    C = np.diag(d) + np.outer(v, v)
    records = [json.loads(s) for s in (BASE / 'draw-results.jsonl').read_text().splitlines()]
    summary = json.loads((BASE / 'summary.json').read_text())
    assert len(records) == len(y) == 128
    assert np.array_equal(z['positive'], y > 0)
    assert np.array_equal(z['archived_positive'], z['archived_flux'] > 0)
    assert np.allclose(y, h + z['latent'][:, None]*v + z['independent_noise']*np.sqrt(d), rtol=1e-14, atol=1e-13)
    gaussian_gap = []
    summary_gap = []
    likelihood_gap = []
    estimator_gap = []
    # Fixed positions plus every reported zero boundary. Explicitly retrospective.
    selected = sorted({0, 31, 63, 95, 127} | {
        r['draw'] for r in records if r['estimators']['full_pattern_conditioned']['amplitude'] == 0})
    write(OUT / 'review-design.json', {
        'status': 'Post-outcome independent verification, no estimator changes',
        'inputs_manifest_sha256': sha(BASE / 'manifest.json'),
        'executable_sha256': sha(Path(__file__)),
        'method': 'Direct integration of Gaussian latent density times observed-row normal densities and missing-row CDFs; no posterior-latent reduction or Hermite nodes. Direct matrix Gaussian solutions. All reported candidate likelihoods, all summaries, and independent 65-point bounded search for selected draws.',
        'selected_draws': selected,
        'selection': 'five fixed index positions plus all reported zero boundaries; not an unbiased scientific subsample',
        'gates': {'amplitude_absolute': 1e-5, 'loglike_absolute': 1e-9, 'arithmetic_absolute': 1e-10},
    })
    for i, r in enumerate(records):
        assert r['draw'] == i and np.array_equal(r['positive_mask'], y[i] > 0)
        p = y[i] > 0
        for name, mask in [('signed_gaussian', np.ones(74, bool)),
                           ('naive_positive_gaussian', p),
                           ('archive_mask_gaussian', z['archived_positive'])]:
            cc = C[np.ix_(mask, mask)]
            hh, yy = h[mask], y[i, mask]
            a = np.clip(hh @ np.linalg.solve(cc, yy) / (hh @ np.linalg.solve(cc, hh)), 0, 4)
            gaussian_gap.append(abs(a-r['estimators'][name]['amplitude']))

        def integrate(log_density):
            atzero = log_density(0.0)
            val, err = quad(lambda u: np.exp(log_density(u)-atzero), -12, 12,
                            epsabs=2e-12, epsrel=2e-12, limit=100, points=[0])
            assert val > 0 and err/val < 1e-9
            return atzero+np.log(val)

        def direct(a, kind):
            mu = a*h
            def joint(u):
                observed = -.5*np.sum(np.log(2*np.pi*d[p]) + (y[i,p]-mu[p]-v[p]*u)**2/d[p])
                missing = np.sum(log_ndtr((-mu[~p]-v[~p]*u)/np.sqrt(d[~p])))
                return -.5*u*u-.5*np.log(2*np.pi)+observed+missing
            answer = integrate(joint)
            if kind == 'full_pattern_conditioned':
                signs = np.where(p, 1., -1.)
                answer -= integrate(lambda u: -.5*u*u-.5*np.log(2*np.pi)+
                                    np.sum(log_ndtr(signs*(mu+v*u)/np.sqrt(d))))
            return float(answer)

        for name in ['joint_censored', 'full_pattern_conditioned']:
            e = r['estimators'][name]
            for candidate in e['candidate_checks']:
                ll = direct(candidate['amplitude'], name)
                likelihood_gap.append(abs(ll-candidate['hermite']['256']))
            if i in selected:
                grid = np.linspace(0., 4., 65)
                vals = np.array([direct(a, name) for a in grid])
                candidates = [(vals[0], 0.), (vals[-1], 4.)]
                for j in range(1, 64):
                    if vals[j] >= vals[j-1] and vals[j] >= vals[j+1]:
                        fit = minimize_scalar(lambda a: -direct(a, name),
                                              bounds=(grid[j-1], grid[j+1]), method='bounded',
                                              options={'xatol': 2e-10})
                        assert fit.success
                        candidates.append((-fit.fun, fit.x))
                _, best = max(candidates)
                estimator_gap.append({'draw': i, 'estimator': name, 'independent_amplitude': float(best),
                                      'absolute_gap': abs(best-e['amplitude'])})

    arrays = {}
    for name, reported in summary['estimators'].items():
        amplitudes = np.array([r['estimators'][name]['amplitude'] for r in records])
        arrays[name] = amplitudes
        summary_gap += [abs(amplitudes.mean()-1-reported['amplitude_bias']),
                        abs(amplitudes.std(ddof=1)/np.sqrt(128)-reported['amplitude_bias_mcse'])]
        assert np.count_nonzero(amplitudes == 0) == reported['zero_boundary']
        assert np.count_nonzero(amplitudes == 4) == reported['upper_boundary']
        assert np.count_nonzero(amplitudes > 0) == reported['log_distance_finite_count']
        finite_dist = -2.5*np.log10(amplitudes[amplitudes > 0])
        summary_gap += [abs(finite_dist.mean()-reported['log_distance_mean_finite_only']),
                        abs(finite_dist.std(ddof=1)/np.sqrt(len(finite_dist))-reported['log_distance_mcse_finite_only'])]
    for name, reported in summary['paired_amplitude_changes'].items():
        a = arrays[name.removesuffix('_minus_signed')]-arrays['signed_gaussian']
        summary_gap += [abs(a.mean()-reported['mean']), abs(a.std(ddof=1)/np.sqrt(len(a))-reported['mcse'])]
    assert max(gaussian_gap + summary_gap) < 1e-10
    assert max(likelihood_gap) < 1e-9
    assert max(x['absolute_gap'] for x in estimator_gap) < 1e-5
    result = {
        'pass': True, 'verified_manifest_entries': len(manifest['files_sha256']),
        'direct_matrix_amplitudes': len(gaussian_gap), 'max_gaussian_amplitude_gap': max(gaussian_gap),
        'direct_latent_candidate_likelihoods': len(likelihood_gap), 'max_candidate_loglike_gap': max(likelihood_gap),
        'independent_searches': estimator_gap, 'max_summary_gap': max(summary_gap),
        'zero_conditional_draws': [r['draw'] for r in records if r['estimators']['full_pattern_conditioned']['amplitude'] == 0],
        'scope': 'Synthetic experiment implementation checked; does not validate observed noise or physical corrections.',
    }
    write(OUT / 'result.json', result)
    (OUT / 'executed-source.py').write_bytes(Path(__file__).read_bytes())
    write(OUT / 'manifest.json', {'files_sha256': {p.name: sha(p) for p in OUT.iterdir() if p.is_file()}})
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
