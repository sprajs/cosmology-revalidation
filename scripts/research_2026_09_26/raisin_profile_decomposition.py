"""Post-outcome conditional geometry of the first predeclared RAISIN profile."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'runs/research_2026_09_26/astra_design/raisin_signed_refit/fixed-c-profile/global-profile'
OUT = ROOT / 'runs/research_2026_09_26/raisin_profile_decomposition'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def dump(p, obj):
    p.write_text(json.dumps(obj, indent=2) + '\n')


def prepare():
    OUT.mkdir(exist_ok=False)
    dump(OUT / 'protocol.json', {
        'status': 'Descriptive follow-up; first A/B shift and best-fit coordinates already seen.',
        'question': 'Separate amplitude-only response from shape/extinction flexibility and decompose common-C Gaussian metric by added rows.',
        'method': 'At saved A/B minima, profile positive amplitude on A and B with common B covariance. Decompose full B quadratic into A marginal plus N|A conditional. No new fits, model calls, domains, uncertainty calibration or outcome-selected exclusion.',
        'limits': 'Empirical AV, fixed header peak, state-frozen covariance, one object. Geometry is not a likelihood-ratio comparison across row sets, physical correction, dust inference or cosmology.',
        'inputs_sha256': {str(p.relative_to(ROOT)): sha(p) for p in [
            SOURCE / 'result.json', SOURCE / 'native-profiles.npz', Path(__file__)
        ]},
    })


def run():
    protocol = json.loads((OUT / 'protocol.json').read_text())
    for p, h in protocol['inputs_sha256'].items():
        assert sha(ROOT / p) == h
    dat = np.load(SOURCE / 'native-profiles.npz')
    summary = json.loads((SOURCE / 'result.json').read_text())
    ix = dat['A_indices']; ni = np.setdiff1d(np.arange(len(dat['data_flux'])), ix)
    C = dat['covariance_B']; y = dat['data_flux']; dref = float(dat['reference_parameters'][0])
    assert np.all(y[ix] > 0) and np.all(y[ni] < 0)
    A = C[np.ix_(ix, ix)]; NA = C[np.ix_(ni, ix)]
    transform = np.linalg.solve(A, NA.T).T
    conditional = C[np.ix_(ni, ni)] - transform @ NA.T
    L = np.linalg.cholesky(C); LA = np.linalg.cholesky(A); LN = np.linalg.cholesky(conditional)
    metrics = {}; errors = []
    for origin in ['A', 'B']:
        point = summary['metrics']['Banchor_' + origin]['finest']
        key = np.array([point['shape'], point['AV']])
        where = np.flatnonzero(np.all(dat['coordinates'] == key, axis=1))
        assert len(where) == 1
        h = dat['model_means'][where[0]]
        for fit_arm, ids, chol in [('A', ix, LA), ('B', np.arange(len(y)), L)]:
            hw = np.linalg.solve(chol, h[ids]); yw = np.linalg.solve(chol, y[ids])
            a = float(hw @ yw / (hw @ hw)); assert a > 0
            residual = y - a * h
            ar = np.linalg.solve(LA, residual[ix])
            nr = np.linalg.solve(LN, residual[ni] - transform @ residual[ix])
            br = np.linalg.solve(L, residual)
            qa, qn, qb = float(ar @ ar), float(nr @ nr), float(br @ br)
            errors.append(abs(qb - qa - qn))
            D = dref - 2.5 * np.log10(a)
            metrics[origin + '_shapeAV_' + fit_arm + '_amplitude'] = {
                'shape': point['shape'], 'AV': point['AV'], 'amplitude': a,
                'DLMAG': float(D), 'Q_positive_marginal': qa,
                'Q_negative_conditional': qn, 'Q_full': qb,
            }
            if fit_arm == origin:
                assert abs(D - point['DLMAG']) < 1e-10
                assert abs((qa if origin == 'A' else qb) - point['Q']) < 1e-9
    aa = metrics['A_shapeAV_A_amplitude']; ab = metrics['A_shapeAV_B_amplitude']
    bb = metrics['B_shapeAV_B_amplitude']
    result = {
        'scope': protocol['limits'], 'n_positive': len(ix), 'n_negative': len(ni),
        'points': metrics,
        'amplitude_only_B_minus_A_at_A_shapeAV': ab['DLMAG'] - aa['DLMAG'],
        'additional_shapeAV_response_with_B_data': bb['DLMAG'] - ab['DLMAG'],
        'total_B_minus_A': bb['DLMAG'] - aa['DLMAG'],
        'Q_B_reduction_A_best_to_B_best': aa['Q_full'] - bb['Q_full'],
        'positive_Q_cost_at_B_best': bb['Q_positive_marginal'] - aa['Q_positive_marginal'],
        'negative_conditional_Q_gain_at_B_best': aa['Q_negative_conditional'] - bb['Q_negative_conditional'],
        'max_block_quadratic_closure_error': max(errors),
        'protocol_sha256': sha(OUT / 'protocol.json'),
    }
    rows = []
    for j in ni:
        rows.append({'band': str(dat['band'][j]), 'MJD': float(dat['MJD'][j]),
                     'flux': float(y[j]), 'quoted_error': float(dat['data_fluxerr'][j]),
                     'quoted_SNR': float(y[j] / dat['data_fluxerr'][j])})
    with (OUT / 'negative-epochs.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    dump(OUT / 'result.json', result)
    (OUT / 'executed_source.py').write_bytes(Path(__file__).read_bytes())
    dump(OUT / 'manifest.json', {'files_sha256': {
        str(p.relative_to(OUT)): sha(p) for p in OUT.rglob('*')
        if p.is_file() and p.name != 'manifest.json'}})
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    import sys
    {'prepare': prepare, 'run': run}[sys.argv[1]]()
