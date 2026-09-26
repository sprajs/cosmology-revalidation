"""Photon-count integrals on the union of piecewise-linear table knots."""

import numpy as np
from numpy.polynomial.legendre import leggauss


def counts(wave, sed, band, order):
    lo, hi = band[0, 0], band[-1, 0]
    assert wave[0] <= lo and wave[-1] >= hi
    if order == 0:
        x = np.linspace(lo, hi, int(np.ceil((hi - lo) / 0.5)) + 1)
        return np.trapezoid(
            x * np.interp(x, band[:, 0], band[:, 1]) * np.interp(x, wave, sed), x
        )
    knots = np.unique(np.r_[band[:, 0], wave[(wave > lo) & (wave < hi)]])
    mid = (knots[1:] + knots[:-1]) / 2
    half = (knots[1:] - knots[:-1]) / 2
    nodes, weights = leggauss(order)
    x = mid[:, None] + half[:, None] * nodes
    value = x * np.interp(x, band[:, 0], band[:, 1]) * np.interp(x, wave, sed)
    return float(np.sum(half * (value @ weights)))
