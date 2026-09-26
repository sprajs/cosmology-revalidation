"""Gaussian sign/censor likelihood for C=diag(d)+v v^T; no physical noise claim."""
from functools import lru_cache
import numpy as np
from scipy.special import erfcx, log_ndtr, logsumexp, roots_hermitenorm


def inverse_mills(t):
    """phi(t)/Phi(t), with stable negative-tail evaluation."""
    t = np.asarray(t, float)
    value = np.zeros_like(t)
    negative = t < 0
    value[negative] = np.sqrt(2 / np.pi) / erfcx(-t[negative] / np.sqrt(2))
    moderate = (t >= 0) & (t < 38)
    u = t[moderate]
    value[moderate] = np.exp(-.5*u*u - .5*np.log(2*np.pi) - log_ndtr(u))
    return value


@lru_cache(None)
def nodes(order):
    x, w = roots_hermitenorm(order)
    return x, np.log(w) - 0.5 * np.log(2 * np.pi)


def log_orthant(mean, variance, loading, signs, order=64):
    """P(signs * Y > 0) in the explicitly supported weak-factor domain.

    Scalar and independent cases are analytic. The multirow quadrature is
    restricted to total factor precision <= .5 and latent log-density mode
    within [-4,4]. Other cases require a separately validated integrator.
    Agreement between quadrature orders alone does not establish accuracy.
    """
    mean, variance, loading, signs = [np.asarray(x, float) for x in (mean, variance, loading, signs)]
    assert mean.shape == variance.shape == loading.shape == signs.shape
    assert np.all(variance > 0) and np.all(np.abs(signs) == 1)
    assert all(np.all(np.isfinite(x)) for x in (mean, variance, loading, signs))
    if len(mean) == 0:
        return 0.0
    if len(mean) == 1:
        return float(log_ndtr(signs[0] * mean[0] / np.sqrt(variance[0] + loading[0]**2)))
    if np.all(loading == 0):
        return float(log_ndtr(signs * mean / np.sqrt(variance)).sum())
    factor_precision = np.sum(loading**2 / variance)
    if factor_precision > .5:
        raise ValueError('Strong-factor orthants are outside this validated quadrature domain')
    a = signs * mean / np.sqrt(variance); b = signs * loading / np.sqrt(variance)
    def score(z):
        t = a + b*z
        mills = inverse_mills(t)
        return -z + np.sum(b*mills)
    if score(-4.) < 0 or score(4.) > 0:
        raise ValueError('Latent mode outside [-4,4]; require a tail-adapted validated integrator')
    z, lw = nodes(order)
    t = signs[:, None] * (mean[:, None] + loading[:, None] * z) / np.sqrt(variance[:, None])
    return float(logsumexp(lw + log_ndtr(t).sum(axis=0)))


def log_gaussian(y, mean, variance, loading):
    y, mean, variance, loading = [np.asarray(x, float) for x in (y, mean, variance, loading)]
    assert y.shape == mean.shape == variance.shape == loading.shape
    assert np.all(variance > 0)
    residual = y - mean
    precision_u = 1 + np.sum(loading * loading / variance)
    projection = np.sum(loading * residual / variance)
    q = np.sum(residual * residual / variance) - projection * projection / precision_u
    logdet = np.log(variance).sum() + np.log(precision_u)
    return float(-0.5 * (len(y) * np.log(2 * np.pi) + logdet + q))


def log_censored(y_observed, mean, variance, loading, observed, censored_signs, order=64):
    """Joint density of observed values and signs of all other eligible rows.

    observed is a boolean mask in the full schedule. censored_signs is ordered
    by ~observed. This includes sign-pattern information; no extra truncation
    normalizer is multiplied in. For negative omitted rows use signs=-1.
    """
    observed = np.asarray(observed, bool)
    mean, variance, loading = [np.asarray(x, float) for x in (mean, variance, loading)]
    y_observed = np.asarray(y_observed, float)
    assert mean.shape == variance.shape == loading.shape == observed.shape
    assert len(y_observed) == int(observed.sum())
    assert np.all(variance > 0)
    residual = y_observed - mean[observed]
    vu = 1 / (1 + np.sum(loading[observed] ** 2 / variance[observed]))
    mu = vu * np.sum(loading[observed] * residual / variance[observed])
    lp = log_gaussian(y_observed, mean[observed], variance[observed], loading[observed])
    return lp + log_orthant(mean[~observed] + loading[~observed] * mu,
                            variance[~observed], loading[~observed] * np.sqrt(vu), censored_signs, order)


def log_positive_conditioned(y, mean, variance, loading, order=64):
    """Retained-value density conditional on all these rows being positive."""
    if np.any(np.asarray(y) <= 0):
        return -np.inf
    return log_gaussian(y, mean, variance, loading) - log_orthant(
        mean, variance, loading, np.ones(len(y)), order)


def decompose_covariance(covariance, relative_tolerance=1e-12):
    """Verify rather than assume the diagonal-plus-one-factor representation."""
    C = np.asarray(covariance, float)
    assert C.ndim == 2 and C.shape[0] == C.shape[1]
    assert np.linalg.norm(C-C.T) <= relative_tolerance * np.linalg.norm(C)
    off = C - np.diag(np.diag(C))
    if np.max(np.abs(off)) == 0:
        return np.diag(C).copy(), np.zeros(len(C))
    i, j = np.unravel_index(np.argmax(np.abs(off)), C.shape)
    candidates = [k for k in range(len(C)) if k not in (i, j) and C[j, k] != 0 and C[i,k] != 0]
    if not candidates:
        raise ValueError('Covariance factor reconstruction needs a nonzero off-diagonal triangle; provide a verified factor explicitly for sparse/two-row cases')
    k = candidates[0]
    square = C[i, j] * C[i, k] / C[j, k]
    assert square > 0
    vi = np.sqrt(square)
    v = C[:, i] / vi; v[i] = vi
    d = np.diag(C) - v*v
    assert np.all(d > 0)
    residual = C - np.diag(d) - np.outer(v, v)
    assert np.linalg.norm(residual) <= relative_tolerance * np.linalg.norm(C)
    return d, v
