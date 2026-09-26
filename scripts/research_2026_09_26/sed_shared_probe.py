"""Fixed constructed mean probes in an existing shared predictive distribution."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import solve_triangular

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'runs/research_2026_09_26/astra_design'
NONLINEAR = ROOT / 'runs/research_2026_09_26/sed_nonlinear_validation/resolved'
OUT = ROOT / 'runs/research_2026_09_26/sed_shared_probe'
PROTOCOL = ROOT / 'docs/research-2026-09-26/sed-shared-protocol.md'
INPUTS = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def record(path):
    INPUTS[str(Path(path).relative_to(ROOT))] = sha(path)
    return path


def arrays(path, keys):
    with np.load(record(path), allow_pickle=False) as data:
        return {key: data[key] for key in keys}


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def inverse_functions(a, covariance):
    small = np.linalg.inv(covariance) + a.T @ a
    q, r = np.linalg.qr(a, mode='reduced')
    qr_cov = np.eye(r.shape[0]) + r @ covariance @ r.T
    def woodbury(v):
        return v - a @ np.linalg.solve(small, a.T @ v)
    def qr(v):
        projected = q.T @ v
        return v - q @ projected + q @ np.linalg.solve(qr_cov, projected)
    return woodbury, qr


def design():
    OUT.mkdir(exist_ok=False)
    record(__file__); record(PROTOCOL)
    shared = arrays(BASE / 'expanded12/validation-projected-modes.npz',
                    ['CID', 'row_offsets', 'calibration_design', 'observer_design'])
    post = arrays(BASE / 'expanded12/posterior-models.npz',
                  ['systematics_only_discovery_posterior_covariance'])
    coefficients = arrays(BASE / 'validation1020/frozen-discovery-coefficients.npz', ['basis_mean'])
    a = shared['calibration_design']
    covariance = post['systematics_only_discovery_posterior_covariance']
    ids = shared['CID'].astype(str)
    assert len(ids) == len(set(ids)) == 1020 and a.shape[1] == 12
    output = {'observer': [], 'sed': []}
    projection_lost = {'observer': 0., 'sed': 0., 'difference': 0.}
    inputs_verified = 0
    for i, cid in enumerate(ids):
        object_path = BASE / f'validation1020/objectives/objective_{cid}.npz'
        cache_path = BASE / f'validation1020/analysis/objects/{cid}.npz'
        info = json.loads(record(NONLINEAR / f'objects/{cid}.json').read_text())
        assert info['all_gates_pass'] and info['CID'] == cid
        for path in (object_path, cache_path):
            assert sha(path) == info['input_sha256'][str(path.relative_to(ROOT))]
            inputs_verified += 1
        pred_path = NONLINEAR / f'objects/{cid}-predictions.npz'
        assert sha(pred_path) == info['prediction_sha256']
        predictions = arrays(pred_path, ['native_noiseless__' + x for x in ('nominal', 'observer', 'sed')])
        nominal = arrays(object_path, ['MJD', 'band', 'model_flux', 'frozen_flux_covariance'])
        cache = arrays(cache_path, ['MJD', 'band', 'jacobian_flux', 'exact_covariance', 'native_flux_model'])
        order = sorted(range(len(nominal['MJD'])), key=lambda j: (nominal['MJD'][j], nominal['band'][j]))
        assert np.array_equal(nominal['MJD'][order], cache['MJD'])
        assert np.array_equal(nominal['band'][order], cache['band'])
        assert np.array_equal(nominal['frozen_flux_covariance'][np.ix_(order, order)], cache['exact_covariance'])
        base = predictions['native_noiseless__nominal'][order]
        assert np.allclose(base, cache['native_flux_model'], rtol=1e-9, atol=1e-9)
        chol = np.linalg.cholesky(cache['exact_covariance'])
        jac = cache['jacobian_flux'].copy()
        jac[:, 0] = -.4 * np.log(10) * nominal['model_flux'][order]
        u, s, _ = np.linalg.svd(solve_triangular(chol, jac, lower=True), full_matrices=True)
        assert (s > s[0] * 1e-10).sum() == 4
        q = u[:, 4:]
        lo, hi = shared['row_offsets'][i:i+2]
        assert hi - lo == q.shape[1]
        white = {}
        for mode in ('observer', 'sed'):
            white[mode] = solve_triangular(chol, predictions['native_noiseless__'+mode][order] - base, lower=True)
            vector = q.T @ white[mode]
            output[mode].append(vector)
            projection_lost[mode] += float(np.sum((u[:, :4].T @ white[mode])**2))
        projection_lost['difference'] += float(np.sum((u[:, :4].T @ (white['observer'] - white['sed']))**2))
    observer = np.concatenate(output['observer']); sed = np.concatenate(output['sed'])
    old = shared['observer_design'] @ coefficients['basis_mean']
    apply, qr = inverse_functions(a, covariance)
    vectors = np.column_stack([observer, sed, observer-sed, old])
    assert np.max(abs(apply(vectors) - qr(vectors))) < 1e-10
    info_plain = vectors.T @ vectors
    info_shared = vectors.T @ apply(vectors)
    assert np.linalg.eigvalsh(info_plain-info_shared).min() > -1e-9
    np.savez_compressed(OUT / 'design.npz', observer=observer, sed=sed, original_observer=old,
                        calibration_design=a, posterior_covariance=covariance, CID=ids,
                        row_offsets=shared['row_offsets'])
    result = {'objects': len(ids), 'projected_dimension': len(observer),
              'vector_order': ['nonlinear_noiseless_observer', 'nonlinear_noiseless_sed', 'difference', 'original_linear_observer'],
              'unshared_gram': info_plain.tolist(), 'shared_gram': info_shared.tolist(),
              'nominal_projection_removed_squared_norm': projection_lost,
              'qr_woodbury_max_error': float(np.max(abs(apply(vectors)-qr(vectors)))),
              'original_input_hashes_verified': inputs_verified,
              'input_sha256': INPUTS, 'design_sha256': sha(OUT/'design.npz'),
              'scope': 'Outcome-free geometry conditional on existing discovery-designed directions and frozen shared covariance.'}
    write(OUT / 'design.json', result)
    print(json.dumps({k: result[k] for k in ('unshared_gram', 'shared_gram', 'nominal_projection_removed_squared_norm')}, indent=2))


def score():
    saved = json.loads((OUT / 'design.json').read_text())
    for p, digest in saved['input_sha256'].items():
        assert sha(ROOT / p) == digest, p
    assert sha(OUT/'design.npz') == saved['design_sha256']
    with np.load(OUT/'design.npz') as d:
        data = dict(d)
    a = data['calibration_design']; covariance = data['posterior_covariance']
    shared = arrays(BASE/'expanded12/validation-projected-modes.npz', ['residual', 'CID'])
    post = arrays(BASE/'expanded12/posterior-models.npz', ['systematics_only_discovery_posterior_mean'])
    assert np.array_equal(shared['CID'], data['CID'])
    centered = shared['residual'] - a @ post['systematics_only_discovery_posterior_mean']
    apply, qr = inverse_functions(a, covariance)
    assert np.max(abs(apply(centered)-qr(centered))) < 1e-10
    result = {}
    for mode in ('observer', 'sed', 'original_observer'):
        v = data[mode]; information = float(v @ apply(v)); matched = float(v @ apply(centered))
        result[mode] = {'I': information, 'M': matched, 'G': matched - .5*information}
    original = json.loads(record(BASE/'expanded12/result.json').read_text())['original_fixed_observer_after_systematic_null_conditioning']
    error = max(abs(result['original_observer'][key] - original[field]) for key, field in
                [('I', 'information'), ('M', 'centered_matched_product'), ('G', 'full_fixed_direction_shift_gain')])
    assert error < 1e-9
    write(OUT/'score.json', {'probes': result, 'sed_minus_observer_G': result['sed']['G']-result['observer']['G'],
                           'original_score_reproduction_max_error': error, 'design_manifest_sha256': sha(OUT/'design.json'),
                           'input_sha256': INPUTS, 'scope': 'Fixed additive probes in the same shared predictive distribution; not jointly refitted physical model evidence.'})
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('phase', choices=['design', 'score'])
    args = parser.parse_args()
    (design if args.phase == 'design' else score)()
