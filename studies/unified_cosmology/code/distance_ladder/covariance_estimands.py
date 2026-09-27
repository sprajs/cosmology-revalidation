"""Independent covariance-estimand audit; no release repair or cosmology fit.

This uses a general-LU/normal-matrix Schur calculation, independently of the
producer's whitened least squares and the prior review's data-space QR/SVD.
The extra-host-noise construction is a counterexample, not an inferred term.
"""
import hashlib
import importlib.metadata
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy import linalg

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
OUT = ROOT / 'studies/unified_cosmology/results/distance_ladder'
LADDER = ROOT / '.work/unified-cosmology/distance-ladder'
CAL = ROOT / '.work/unified-cosmology/calibration-interface'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pin_records(records):
    for name, digest in records.items():
        assert sha(ROOT/name) == digest, name


def covariance(path):
    with path.open() as stream:
        n = int(stream.readline())
        raw = np.loadtxt(stream).reshape(n, n)
    assert np.isfinite(raw).all()
    return raw


def mathematical_controls():
    # Arbitrary positive definite covariances, not fits to observed outcomes.
    H = np.array([[.09, .02], [.02, .04]])
    v = np.arange(1., 8.) / 20
    R = np.diag(np.linspace(.015, .035, 7)) + np.outer(v, v)
    J = np.zeros((7, 2)); J[:2, 0] = 1; J[2:4, 1] = 1
    S = R + J@H@J.T
    one = np.ones(7); r = np.array([.1, .3, -.2, .4, -.1, .5, -.3])
    Ri = linalg.solve(R, np.eye(7)); Hi = linalg.solve(H, np.eye(2))
    B = np.column_stack([one, J]); F = B.T@Ri@B; F[1:, 1:] += Hi
    b = B.T@Ri@r
    joint = r@Ri@r - b@linalg.solve(F, b)
    joint_norm = np.linalg.slogdet(R)[1] + np.linalg.slogdet(H)[1] + np.linalg.slogdet(F)[1]
    Si = linalg.solve(S, np.eye(7)); a = one@Si@one
    P = Si - np.outer(Si@one, Si@one)/a
    embedded = r@P@r
    embedded_norm = np.linalg.slogdet(S)[1] + np.log(a)
    w = Si@one/a; T = np.eye(7)-np.outer(one, w)
    projection_error = float(np.max(abs(T@S@T.T - (S-np.ones((7, 7))/a))))
    # Extra independent host-level SN error gives the observed *type* of pattern.
    # It is not evidence that this term was present in the released analysis.
    D = np.diag([.03, .07]); extra = J@D@J.T
    assert linalg.eigvalsh(extra).min() > -1e-14
    assert extra[0, 1] == .03 and extra[2, 3] == .07 and extra[0, 2] == 0
    K = np.array([[1, -1, 0, 0, 0, 0, 0], [0, 0, 1, -1, 0, 0, 0]])
    assert np.array_equal(K@J, np.zeros((2, 2)))
    cancellation = float(np.max(abs(K@(S+extra)@K.T - K@R@K.T)))
    errors = {'joint_host_and_M_integration_quadratic_error': float(abs(joint-embedded)),
              'joint_host_and_M_integration_logdet_error': float(abs(joint_norm-embedded_norm)),
              'fitted_M_residual_covariance_identity_error': projection_error,
              'sibling_contrast_host_noise_cancellation_error': cancellation}
    assert max(errors.values()) < 1e-11
    assert np.linalg.matrix_rank(P, tol=1e-10) == 6
    return {**errors, 'projected_precision_rank': 6, 'data_dimension': 7,
            'counterexample_extra_host_covariance_psd': True,
            'counterexample_is_empirical_component_estimate': False}


def main():
    factor_path = OUT/'cepheid-factor.json'
    factor = json.loads(factor_path.read_text())
    old_review_path = OUT/'independent-factor-review.json'
    old_review = json.loads(old_review_path.read_text())
    interface_path = OUT/'calibration-interface.json'
    interface = json.loads(interface_path.read_text())
    table_path = OUT/'calibration-table6-review.json'
    table = json.loads(table_path.read_text())
    # Re-identify all prior evidence before using its labels or fit quantities.
    for field in ['source_sha256', 'dependencies_sha256', 'output_sha256']:
        pin_records(factor[field])
    pin_records(old_review['dependency_sha256'])
    assert sha(HERE/'independent_factor_review.py') == old_review['source_sha256']
    for field in ['source_sha256', 'dependencies_sha256', 'output_sha256']:
        pin_records(interface[field])
    for name, entry in interface['input_sources'].items():
        assert sha(ROOT/name) == entry['sha256']
    pin_records(table['input_sha256']); pin_records(table['source_sha256'])
    paths = [LADDER/f'all{c}_shoes_ceph_topantheonwt6.0_112221.fits' for c in 'ylc']
    ledger_path = CAL/'calibrator-host-ledger.csv'
    dat_path = CAL/'Pantheon+SH0ES.dat'
    stat_path = CAL/'Pantheon+SH0ES_STATONLY.cov'
    total_path = CAL/'Pantheon+SH0ES_STAT+SYS.cov'
    inputs = paths + [factor_path, old_review_path, interface_path, table_path,
                      ledger_path, dat_path, stat_path, total_path, OUT/'calibration-source-review.json']
    hashes = {str(p.relative_to(ROOT)): sha(p) for p in inputs}
    y, A, C = [np.asarray(fits.getdata(p), dtype=float) for p in paths]; A = A.T
    assert A.shape == (3492, 47) and C.shape == (3492, 3492)
    keep = A[:, 42] == 0; discard = ~keep
    cols = np.flatnonzero(np.any(A[keep] != 0, axis=0))
    assert keep.sum() == 3138 and discard.sum() == 354
    assert cols.tolist() == [i for i in range(47) if i not in (42, 46)]
    assert np.count_nonzero(C[np.ix_(keep, discard)]) == 0
    X = A[np.ix_(keep, cols)]
    # Full retained covariance solved by general LU; no whitened producer call.
    solved = linalg.solve(C[np.ix_(keep, keep)], np.column_stack([y[keep], X]), assume_a='gen')
    N = X.T@solved[:, 1:]; b = X.T@solved[:, 0]
    Q = N[:37, :37] - N[:37, 37:]@linalg.solve(N[37:, 37:], N[37:, :37])
    t = b[:37] - N[:37, 37:]@linalg.solve(N[37:, 37:], b[37:])
    H = linalg.solve(Q, np.eye(37)); mu = linalg.solve(Q, t)
    Hfixed = linalg.solve(N[:37, :37], np.eye(37))
    savedH = np.asarray(factor['covariance_distance_modulus'])
    errors = {'mean_max_absolute_mag': float(np.max(abs(mu-factor['mean_distance_modulus']))),
              'covariance_max_absolute_mag2': float(np.max(abs(H-savedH)))}
    assert errors['mean_max_absolute_mag'] < 1e-8 and errors['covariance_max_absolute_mag2'] < 1e-10
    nuisance_increase_min = float(linalg.eigvalsh((H-Hfixed + (H-Hfixed).T)/2).min())
    assert nuisance_increase_min > -1e-10
    ledger = pd.read_csv(ledger_path)
    data = pd.read_csv(dat_path, sep=r'\s+')
    S = covariance(stat_path)
    assert len(data) == len(S) == 1701 and len(ledger) == 77
    assert ledger.physical_SN.nunique() == 42 and ledger.host.nunique() == 37
    assert np.array_equal(ledger.original_row_zero_based, np.flatnonzero(data.IS_CALIBRATOR.eq(1)))
    siblings = []
    for original in table['host_comparisons']:
        h = original['host']; i = factor['host_order'].index(h)
        sub = ledger.loc[ledger.host.eq(h)]
        values = [S[int(a.original_row_zero_based), int(b.original_row_zero_based)]
                  for _, a in sub.iterrows() for _, b in sub.iterrows() if a.physical_SN != b.physical_SN]
        assert np.ptp(values) == 0
        value = float(values[0]); excess = value-H[i, i]
        assert excess > 0
        assert abs(value-original['distinct_sibling_STATONLY_covariance_mag2']) < 1e-14
        siblings.append({'host': h, 'physical_SNe': sorted(sub.physical_SN.unique().tolist()),
                         'STATONLY_shared_mag2': value, 'SN_free_host_variance_mag2': float(H[i, i]),
                         'ratio': float(value/H[i, i]), 'positive_unattributed_excess_mag2': float(excess),
                         'extra_host_shared_sigma_if_attributed_entirely_to_independent_SN_host_noise_mag': float(np.sqrt(excess)),
                         'host_variance_with_eight_nuisances_fixed_mag2': float(Hfixed[i, i]),
                         'nuisance_marginalization_increment_mag2': float(H[i, i]-Hfixed[i, i]),
                         'excess_over_fixed_nuisance_host_variance': float(excess/Hfixed[i, i]),
                         'printed_Table6b_sigma_mag': original['table6b_rows'][0]['sigma_b']})
    assert len(siblings) == 4
    off = []
    for i in range(37):
        a = ledger.loc[ledger.host_index.eq(i), 'original_row_zero_based'].to_numpy()
        for j in range(i+1, 37):
            b = ledger.loc[ledger.host_index.eq(j), 'original_row_zero_based'].to_numpy()
            block = S[np.ix_(a, b)]; assert np.ptp(block) == 0
            off.append({'host1': factor['host_order'][i], 'host2': factor['host_order'][j],
                        'released': float(block[0, 0]), 'host_factor': float(H[i, j])})
    assert len(off) == 666
    def compare(rows):
        observed = np.array([v['released'] for v in rows]); base = np.array([v['host_factor'] for v in rows])
        return {'pairs': len(rows), 'descriptive_least_squares_multiplier_through_zero': float(base@observed/(base@base)),
                'rms_minus_one_host_covariance_mag2': float(np.sqrt(np.mean((observed-base)**2))),
                'rms_minus_two_host_covariances_mag2': float(np.sqrt(np.mean((observed-2*base)**2))),
                'ratio_median': float(np.median(observed/base)),
                'released_zero_pairs': int(np.count_nonzero(observed == 0))}
    raw = covariance(total_path); select = data.IS_CALIBRATOR.eq(1) | data.zHD.gt(.01)
    total = (raw[np.ix_(select, select)] + raw.T[np.ix_(select, select)])/2
    one = np.ones(len(total)); a = float(one@linalg.solve(total, one, assume_a='gen'))
    assert abs(a/interface['likelihood']['flat_M_precision']-1) < 1e-11
    result = {'status': 'passed_covariance_estimand_and_counterexample_audit_origin_unidentified',
              'independent_SN_free_factor': {'method': 'General LU on full retained covariance, normal-matrix Schur elimination of all8 nuisances; no producer import.',
                    'full_rows': 3492, 'full_parameters': 47, 'retained_rows': 3138, 'removed_SN_rows': 354,
                    'retained_columns': cols.tolist(), 'host_parameters': 37, 'nuisance_parameters': 8,
                    'retained_removed_covariance_nonzero': 0, **errors,
                    'marginal_minus_fixed_nuisance_covariance_min_eigenvalue_mag2': nuisance_increase_min},
              'siblings': siblings,
              'off_host_comparison': {'all_pairs': compare(off),
                    'excluding_N1365': compare([v for v in off if 'N1365' not in [v['host1'], v['host2']]]),
                    'multiplier_is_descriptive_not_fitted_physical_component': True},
              'flat_M_projection': {'selected_rows': len(total), 'intercept_precision_per_mag2': a,
                    'subtracted_constant_from_every_fitted_residual_covariance_entry_mag2': 1/a,
                    'identity': 'Cov(r-1 Mhat)=S-11^T/(1^T S^-1 1); marginalized precision is S^-1-S^-1 11^T S^-1/(1^T S^-1 1).',
                    'explains_positive_same_host_only_raw_covariance_increment': False},
              'mathematical_controls': mathematical_controls(),
              'conclusions': [
                    'The independently reconstructed H is the marginalized first-two-rung host-distance covariance, not a covariance conditioned on Cepheid nuisance values or on SN data.',
                    'Table6b and H refer to the same kind of SN-free host-distance estimand; the complete released SN-pair residual covariance is a different quantity.',
                    'Under independent host/SN measurements, residual covariance is C_SN+J H J^T. One host variance enters a distinct sibling pair once; duplicating the number of siblings does not multiply it.',
                    'Unknown host-SN cross-covariance adds the two signed cross terms; an independent host-shared SN term J D J^T with diagonal positive D can raise same-host entries while preserving off-host entries.',
                    'The latter is a valid mathematical counterexample, not evidence for an actual released covariance component. The paper prescriptions and missing construction inputs must decide whether such a term is justified.',
                    'Twice the entire H would also double off-host terms and is not this observed pattern. Adding a duplicate host-specific diagonal is compatible with the pattern but is not established as its cause.',
                    'Flat common-M integration changes precision globally and cannot account for positive host-specific increments in the raw released matrix. Sibling differences cancel both Cepheid and additional host-common terms and cannot identify their absolute size.',
                    'No released covariance is repaired or rescaled; no calibrated H0 or cosmological inference follows from this audit.'],
              'dependency_sha256': hashes,
              'source_sha256': {str(Path(__file__).relative_to(ROOT)): sha(__file__)},
              'versions': {k: importlib.metadata.version(k) for k in ['numpy', 'scipy', 'pandas', 'astropy']},
              'new_cosmology_fits': 0, 'released_covariance_edits': 0}
    pin_records(hashes)
    output = OUT/'covariance-estimands.json'
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: result[k] for k in ['status', 'independent_SN_free_factor', 'siblings', 'off_host_comparison', 'flat_M_projection', 'mathematical_controls']}, indent=2))


if __name__ == '__main__':
    main()
