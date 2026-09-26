"""Reproduce absorption geometry and exact latent-colour reparameterization."""

import numpy as np
from scipy.integrate import quad
from lib.dust import slab, result
from lib.records import write_rows

DEFAULTS = {}


def run(out, cfg):
    rows = [
        result(rv, tau)
        for rv in [2.0, 3.1, 4.0]
        for tau in [0.01, 0.1, 0.3, 1.0, 3.0, 10.0]
    ]
    write_rows(out / "geometry.csv", rows)
    errors = [
        abs(quad(lambda depth: np.exp(-tau * depth), 0, 1)[0] - (-np.expm1(-tau) / tau))
        for tau in [0.01, 0.1, 1.0, 3.0, 10.0]
    ]
    if float(slab(0)) != 0 or max(errors) > 1e-12:
        raise RuntimeError("Slab integration or transparent-limit failure")
    age = np.linspace(0, 10, 101)
    e, ci = 0.03 + 0.004 * age, -0.02 + 0.001 * age
    beta, rb, b, k = 2.0, 4.0, -0.03, 0.005
    colour, brightness = ci + e, beta * ci + rb * e + b * age
    ep, cip, bp = e + k * age, ci - k * age, b - (rb - beta) * k
    gap = float(np.max(abs(brightness - (beta * cip + rb * ep + bp * age))))
    if gap > 1e-12:
        raise RuntimeError("Exact colour/brightness reparameterization failed")
    return {
        "integration_error": max(errors),
        "grid": rows,
        "original_age_slope": b,
        "alternative_age_slope": bp,
        "brightness_identity_error": gap,
        "colour_identity_error": float(np.max(abs(colour - cip - ep))),
        "scope": "Constructed absorption-only slab and exact latent-mean degeneracy. Equal galaxy attenuation and SN sightline extinction are not assumed. This is not an observed dust population or a fitted cosmological correction.",
    }
