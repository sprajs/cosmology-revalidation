"""Absorption-only slab attenuation; a constructed geometry, not a fitted law."""

import numpy as np

K = 2.5 / np.log(10)


def slab(tau):
    """Absorption-only uniform-slab attenuation, including the transparent limit."""
    tau = np.asarray(tau, dtype=float)
    if np.any(~np.isfinite(tau)) or np.any(tau < 0):
        raise ValueError("Slab optical depth must be finite and nonnegative")
    transmission = np.ones_like(tau)
    np.divide(-np.expm1(-tau), tau, out=transmission, where=tau != 0)
    attenuation = np.asarray(-2.5 * np.log10(transmission))
    # The logarithm loses relative accuracy when transmission approaches one.
    # -ln[(1-exp(-tau))/tau] = tau/2 - tau^2/24 + tau^4/2880 + O(tau^6).
    thin = tau < 1e-4
    t = tau[thin]
    attenuation[thin] = K * (t / 2 - t * t / 24 + t**4 / 2880)
    return attenuation


def result(rv, tau, limit=0.5):
    av = float(slab(tau))
    ab = float(slab(tau * (1 + 1 / rv)))
    f = min(1, limit / (K * tau * (1 + 1 / rv)))
    # Uniform SN depth u~U[0,1], extinction=K*tau*u; deterministic A_B<limit.
    return {
        "microscopic_RV": rv,
        "tau_V": tau,
        "galaxy_A_V": av,
        "galaxy_A_B": ab,
        "galaxy_effective_RV": av / (ab - av),
        "SN_screen_RV": rv,
        "selected_fraction_AB_lt_0p5": f,
        "unselected_mean_SN_E": K * tau / (2 * rv),
        "selected_mean_SN_E": K * tau * f / (2 * rv),
    }
