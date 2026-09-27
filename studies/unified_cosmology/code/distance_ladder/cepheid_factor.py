"""SN-free compression of the released Cepheid/anchor Gaussian ladder block."""
import json
from decimal import Decimal
import numpy as np
from scipy import linalg
from common import HERE, WORK, OUT, sha, relative, write
from solve import load, gls, NAMES


def decimal_half_width(token):
    return float(Decimal(10) ** Decimal(Decimal(token).as_tuple().exponent)) / 2


def table_groups():
    groups = {}
    for line in (WORK / 'table2.tex').read_text().splitlines():
        if '&' not in line or line.startswith('\\'):
            continue
        fields = [v.strip().rstrip('\\') for v in line.split('&')]
        groups.setdefault(fields[0], []).append(fields)
    assert sum(map(len, groups.values())) == 3130
    return groups


def identify_hosts(A):
    groups = table_groups()
    mapping = []
    for j in range(37):
        rows = np.flatnonzero((A[:, j] != 0) & (A[:, 38] != 0))
        candidates = []
        for name, data in groups.items():
            if len(data) != len(rows):
                continue
            p = np.array([float(v[4]) for v in data])
            ph = np.array([decimal_half_width(v[4]) for v in data])
            m = np.array([float(v[9]) for v in data])
            mh = np.array([decimal_half_width(v[9]) for v in data])
            lp = A[rows, 41]
            lp_ulp = abs(np.spacing(lp.astype(np.float32)).astype(float))
            metal = A[rows, 43]
            m_ulp = abs(np.spacing(metal.astype(np.float32)).astype(float))
            ok_p = (lp >= np.log10(p-ph)-1-2*lp_ulp) & (lp <= np.log10(p+ph)-1+2*lp_ulp)
            ok_m = abs(metal-m) <= mh+2*m_ulp
            if np.all(ok_p & ok_m):
                candidates.append({'host': name, 'maximum_period_rounding_fraction': float(np.max(abs(10**(lp+1)-p)/ph)),
                                   'maximum_metallicity_rounding_fraction': float(np.max(abs(metal-m)/mh))})
        assert len(candidates) == 1, (j, len(rows), candidates)
        mapping.append({'parameter_index_zero_based': j, 'Cepheid_rows': len(rows),
                        'released_row_indices_zero_based': rows.tolist(), **candidates[0]})
    assert len({v['host'] for v in mapping}) == 37
    return mapping


def main():
    y, A, C = load()
    labels = identify_hosts(A)
    keep = A[:, 42] == 0
    columns = np.flatnonzero(np.any(A[keep] != 0, axis=0))
    assert columns.tolist() == [i for i in range(47) if i not in (42, 46)]
    assert np.count_nonzero(C[np.ix_(keep, ~keep)]) == 0
    data, design, noise = y[keep], A[np.ix_(keep, columns)], C[np.ix_(keep, keep)]
    q, cov, checks = gls(data, design, noise)
    mean, host_cov = q[:37], cov[:37, :37]
    white = linalg.solve_triangular(linalg.cholesky(noise, lower=True), design, lower=True)
    wd = linalg.solve_triangular(linalg.cholesky(noise, lower=True), data, lower=True)
    normal = white.T @ white
    nuisance = white[:, 37:]
    marginal_precision = normal[:37, :37] - normal[:37, 37:] @ linalg.solve(normal[37:, 37:], normal[37:, :37], assume_a='pos')
    identity_error = float(np.max(abs(marginal_precision @ host_cov - np.eye(37))))
    assert identity_error < 1e-8
    rng = np.random.default_rng(273310)
    deviations = rng.normal(size=(12, 37)) @ linalg.cholesky(host_cov, lower=True).T
    deviations *= np.tile([.1, 1., 3.], 4)[:, None]
    profile_errors = []
    for delta in deviations:
        fixed = mean + delta
        target = wd - white[:, :37] @ fixed
        fitted, _, rank, _ = linalg.lstsq(nuisance, target, lapack_driver='gelsd')
        assert rank == 8
        residual = target - nuisance @ fitted
        direct = residual @ residual - checks['chi2']
        compressed = delta @ linalg.solve(host_cov, delta, assume_a='pos')
        profile_errors.append(float(abs(direct-compressed)))
    assert max(profile_errors) < 1e-7
    npz = WORK / 'cepheid-host-factor.npz'
    np.savez_compressed(npz, host=np.array([v['host'] for v in labels]), mean_mu=mean, covariance_mu=host_cov)
    base = json.loads((OUT/'reconstruction.json').read_text())
    baseline = np.array([v['fit'] for v in base['parameters_by_index'][:37]])
    result = {
        'status': 'passed_conditional_SN_free_host_factor',
        'data_rows': len(data), 'parameters': len(columns), 'host_parameters': 37, 'nuisance_parameters': 8,
        'removed_SN_rows': int((~keep).sum()), 'removed_parameters': [NAMES[i] for i in (42, 46)],
        'retained_removed_covariance_nonzero': 0,
        'host_order': [v['host'] for v in labels], 'mean_distance_modulus': mean.tolist(),
        'covariance_distance_modulus': host_cov.tolist(),
        'host_mapping': labels,
        'host_sigma_summary_mag': {'min': float(np.sqrt(np.diag(host_cov)).min()), 'median': float(np.median(np.sqrt(np.diag(host_cov)))), 'max': float(np.sqrt(np.diag(host_cov)).max())},
        'maximum_host_mean_change_after_removing_SN_mag': float(np.max(abs(mean-baseline))),
        'validation': {**checks, 'schur_covariance_identity_max_abs': identity_error, 'profile_perturbations': 12,
                       'profile_delta_chi2_absolute_errors': profile_errors},
        'likelihood_form': '-0.5 (mu-mean)^T covariance^-1 (mu-mean); nuisance profile and flat-nuisance marginalization have identical mu dependence. No evidence normalization is claimed.',
        'assumptions': ['Original selected/corrected Cepheids and fixed covariance, period and metallicity design values, Wesenheit law and compressed geometric/MW constraints.',
                        'All eight nuisance coordinates marginalized jointly. Host errors are correlated.',
                        'No SN magnitudes, redshifts, Hubble-flow cosmography or H0 coordinate remain.',
                        'Not fused with Dovekie: cross-release calibration response and physical calibrator linkage must be supplied explicitly.'],
        'source_sha256': {relative(p): sha(p) for p in [HERE/'common.py', HERE/'solve.py', HERE/'cepheid_factor.py', HERE/'cepheid-design.json']},
        'dependencies_sha256': {relative(p): sha(p) for p in [OUT/'acquisition.json', OUT/'reconstruction.json', WORK/'table2.tex']},
        'output_sha256': {relative(npz): sha(npz)}
    }
    write(OUT/'cepheid-factor.json', result)
    print(json.dumps({k: result[k] for k in ['status','data_rows','parameters','host_sigma_summary_mag','maximum_host_mean_change_after_removing_SN_mag','validation']}, indent=2))


if __name__ == '__main__':
    main()
