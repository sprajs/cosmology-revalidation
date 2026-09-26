#!/usr/bin/env python3
"""Independent algebra/identifiability checks; never reads observed SN residuals."""
from pathlib import Path
import hashlib
import json
import platform

import numpy as np
import scipy
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
RNG = np.random.default_rng(2026092607)


def solve_spd(a, b):
    return cho_solve(cho_factor(a, lower=True), b)


def logdet_spd(a):
    chol = np.linalg.cholesky(a)
    return 2 * np.log(np.diag(chol)).sum()


def normal_logpdf(y, mean, cov):
    e = y - mean
    return float(-0.5 * (len(y) * np.log(2 * np.pi) + logdet_spd(cov)
                        + e @ solve_spd(cov, e)))


def low_rank_R(f, g):
    s = np.eye(len(g)) + f
    return float(0.5 * (g @ solve_spd(s, g) - logdet_spd(s)))


def posterior(x, y):
    v = solve_spd(np.eye(x.shape[1]) + x.T @ x, np.eye(x.shape[1]))
    return v @ x.T @ y, v


def predictive_low_rank(xd, yd, xv, yv):
    fd, fv = xd.T @ xd, xv.T @ xv
    gd, gv = xd.T @ yd, xv.T @ yv
    base = -0.5 * (len(yv) * np.log(2 * np.pi) + yv @ yv)
    return float(base + low_rank_R(fd + fv, gd + gv) - low_rank_R(fd, gd))


def predictive_dense_posterior(xd, yd, xv, yv):
    mean, var = posterior(xd, yd)
    return normal_logpdf(yv, xv @ mean, np.eye(len(yv)) + xv @ var @ xv.T)


def predictive_dense_joint(xd, yd, xv, yv):
    x = np.vstack([xd, xv])
    y = np.concatenate([yd, yv])
    return (normal_logpdf(y, np.zeros(len(y)), np.eye(len(y)) + x @ x.T)
            - normal_logpdf(yd, np.zeros(len(yd)), np.eye(len(yd)) + xd @ xd.T))


def independent_dense_checks():
    nd, nv, k, p = 7, 11, 3, 2
    xd = RNG.normal(size=(nd, k + p))
    xv = RNG.normal(size=(nv, k + p))
    yd, yv = RNG.normal(size=nd), RNG.normal(size=nv)
    methods = [predictive_low_rank, predictive_dense_posterior, predictive_dense_joint]
    values = {method.__name__: [method(xd[:, :n], yd, xv[:, :n], yv)
                                for n in (k, k + p)] for method in methods}
    errors = [abs(values[name][j] - values[methods[0].__name__][j])
              for name in values for j in (0, 1)]
    rotation, _ = np.linalg.qr(RNG.normal(size=(k + p, k + p)))
    rotated = predictive_low_rank(xd @ rotation, yd, xv @ rotation, yv)
    wrongly_rotated = predictive_low_rank(xd, yd, xv @ rotation, yv)
    common_rotation_error = abs(rotated - values[methods[0].__name__][1])
    assert max(errors) < 1e-11
    assert common_rotation_error < 1e-11
    # Proper shared-mode posterior has indispensable calibration-observer blocks.
    _, var = posterior(xd, yd)
    return {
        "dimensions": {"discovery": nd, "validation": nv, "calibration": k, "extra": p},
        "log_predictive_density": values,
        "extension_minus_calibration": values[methods[0].__name__][1] - values[methods[0].__name__][0],
        "maximum_method_disagreement": max(errors),
        "shared_orthogonal_rotation_error": common_rotation_error,
        "validation_only_rotation_logscore_change": wrongly_rotated - values[methods[0].__name__][1],
        "maximum_posterior_calibration_extra_cross_covariance": float(np.abs(var[:k, k:]).max()),
    }


def coherent_calibration_toy():
    # Dimensionless projected data. Means exactly 0.5 chosen deterministically;
    # this is a teaching construction, never an estimate from DES.
    nd, nv = 43, 1020
    yd, yv = np.full(nd, 0.5), np.full(nv, 0.5)
    xd, xv = np.ones((nd, 1)), np.ones((nv, 1))
    mean, var = posterior(xd, yd)
    m, v = float(mean[0]), float(var[0, 0])
    nominal_z = yv.mean() * np.sqrt(nv)
    conditional_sd = np.sqrt(1 / nv + v)
    proper_z = (yv.mean() - m) / conditional_sd
    wrong_shared_sign_z = (yv.mean() + m) / conditional_sd
    joint_score = predictive_low_rank(xd, yd, xv, yv)
    iid_baseline = -0.5 * (nv * np.log(2 * np.pi) + yv @ yv)
    sum_marginals = nv * normal_logpdf(np.array([0.5]), np.array([m]), np.array([[1 + v]]))
    xd2, xv2 = np.repeat(xd, 2, axis=1), np.repeat(xv, 2, axis=1)
    mean2, var2 = posterior(xd2, yd)
    extension_gain = predictive_low_rank(xd2, yd, xv2, yv) - joint_score
    corr = var2[0, 1] / np.sqrt(var2[0, 0] * var2[1, 1])
    # An exact likelihood-null coordinate remains prior dominated even with
    # more SNe: the difference between two indistinguishable unit-prior modes.
    null_contrast = np.array([1.0, -1.0]) / np.sqrt(2)
    null_var = float(null_contrast @ var2 @ null_contrast)
    assert abs(null_var - 1) < 1e-12
    # Independent conditional simulations check the one-sided matched statistic.
    draws = 200_000
    ybar = RNG.normal(m, np.sqrt(v), draws) + RNG.normal(0, 1 / np.sqrt(nv), draws)
    conditional_z = (ybar - m) / conditional_sd
    conditional_tail = float(np.mean(conditional_z > norm.isf(0.05)))
    nominal_tail = float(np.mean(ybar * np.sqrt(nv) > norm.isf(0.05)))
    return {
        "dimensionless_toy_only": True,
        "discovery_rows": nd, "validation_rows": nv, "all_observed_toy_values": 0.5,
        "calibration_posterior_mean": m, "calibration_posterior_variance": v,
        "nominal_zero_calibration_Z": float(nominal_z),
        "proper_shared_calibration_conditioned_Z": float(proper_z),
        "wrong_validation_only_mode_sign_Z": float(wrong_shared_sign_z),
        "joint_validation_log_gain_vs_iid_zero": float(joint_score - iid_baseline),
        "sum_separate_one_row_marginal_log_gain_vs_iid_zero": float(sum_marginals - iid_baseline),
        "indistinguishable_extension_predictive_gain": float(extension_gain),
        "two_mode_posterior_means": mean2.tolist(),
        "two_mode_posterior_covariance": var2.tolist(),
        "two_mode_posterior_correlation": float(corr),
        "unidentified_difference_prior_and_posterior_variance": [1.0, null_var],
        "conditional_simulation": {
            "draws": draws,
            "proper_Z_mean": float(conditional_z.mean()),
            "proper_Z_sd": float(conditional_z.std(ddof=1)),
            "proper_Z_one_sided_5percent_exceedance": conditional_tail,
            "binomial_standard_error_at_5percent": float(np.sqrt(0.05 * 0.95 / draws)),
            "nominal_Z_one_sided_5percent_exceedance_under_same_conditioned_null": nominal_tail,
        },
    }


def finite_realization_toy():
    # Nine draws in ten true latent dimensions necessarily miss a direction.
    response = RNG.normal(size=(10, 9))
    covariance = 0.09 * response @ response.T
    u, s, _ = np.linalg.svd(response, full_matrices=True)
    missed = u[:, -1]
    missing_variance = float(missed @ covariance @ missed)
    assert abs(missing_variance) < 1e-12
    return {
        "true_covariance": "I_10, solely for this synthetic construction",
        "realizations": 9,
        "response_weight": 0.3,
        "finite_covariance_rank": int(np.linalg.matrix_rank(covariance)),
        "covariance_eigenvalues": np.linalg.eigvalsh(covariance).tolist(),
        "missed_unit_direction_finite_covariance_variance": missing_variance,
        "missed_unit_direction_true_variance": 1.0,
        "expected_covariance_factor_if_draws_iid_from_true_covariance": 9 * 0.3**2,
        "relative_sd_of_fixed_direction_second_moment_if_iid_Gaussian_draws": float(np.sqrt(2 / 9)),
        "warning": "The independent-Gaussian-draw premise is not asserted for released DES variants.",
    }


def main():
    result = {
        "status": "PASS; algebra and synthetic-model checks only, no observed-data fits",
        "seed": 2026092607,
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "independent_dense_checks": independent_dense_checks(),
        "coherent_calibration_toy": coherent_calibration_toy(),
        "finite_realization_toy": finite_realization_toy(),
    }
    (HERE / "synthetic-checks.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
