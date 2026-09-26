"""Project released BAO onto the flat nonaccelerating shape cone."""

import numpy as np
from lib.bao_shape import load_data, cone_matrices, Cone, tail
from lib.records import write_rows

DEFAULTS = {"draws": 5000, "seed": 260926}


def run(out, config):
    if config["draws"] < 1:
        raise ValueError("draws must be positive")
    z, y, covariance, _, order, _, _ = load_data()
    inequalities, rays, labels = cone_matrices(np.log1p(z))
    cone = Cone(covariance, rays, inequalities)
    statistic, fitted, kkt = cone.fit(y)
    second, _, success = cone.independent(y)
    if not success or abs(statistic - second) > 1e-6 or kkt > 1e-6:
        raise RuntimeError("Independent constrained solvers disagree")
    rng = np.random.default_rng(config["seed"])
    apex = np.array(
        [cone.fit(cone.L @ rng.normal(size=len(y)))[0] for _ in range(config["draws"])]
    )
    np.savez_compressed(out / "null_draws.npz", statistic=apex)
    write_rows(
        out / "projection.csv",
        [
            {"component": i, "observed": a, "projected": b}
            for i, (a, b) in enumerate(zip(y, fitted))
        ],
    )
    return {
        "statistic": statistic,
        "independent_statistic": second,
        "kkt_error": kkt,
        "redshifts": z.tolist(),
        "original_row_order": order,
        "apex_dominating_tail": tail(apex, statistic),
        "scope": "Gaussian released BAO, flat geometry, common free ruler, nonaccelerating expansion across sampled intervals. The apex tail dominates feasible cone-mean tails; it is not a point estimate of present q(0). BGS isotropic DV is omitted from this 12-dimensional anisotropic cone.",
    }
