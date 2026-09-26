"""Finite Gaussian shared-mode inference and distance-response geometry."""

import numpy as np
from lib.paths import DATA
from lib.records import write_rows

DEFAULTS = {}


def run(out, cfg):
    with np.load(
        DATA / "calibration/discovery-projected-modes.npz", allow_pickle=False
    ) as f:
        fd, ud = f["joint_F"].sum(0), f["joint_u"].sum(0)
    with np.load(
        DATA / "calibration/validation-projected-modes.npz", allow_pickle=False
    ) as f:
        fv, uv = f["joint_F"].sum(0), f["joint_u"].sum(0)
    with np.load(DATA / "calibration/distance-responses.npz", allow_pickle=False) as f:
        h = f["contrast"].copy()
    h[12:] /= 0.02
    bandmap = np.array([[1.0, 0, 0], [-1, -1, -1], [0, 1, 0], [0, 0, 1]])
    iso = np.eye(15)
    iso[12:, 12:] = np.linalg.cholesky(2 * 0.02**2 * np.linalg.inv(bandmap.T @ bandmap))
    transforms = {
        "systematics_only": np.eye(15)[:, :12],
        "inherited_observer": np.diag(np.r_[np.ones(12), [0.02] * 3]),
        "isotropic_observer": iso,
    }
    rows = []
    for name, t in transforms.items():
        target = h @ t
        identity = np.eye(t.shape[1])
        f0 = t.T @ fd @ t
        f1 = t.T @ (fd + fv) @ t
        c0 = np.linalg.inv(identity + f0)
        c1 = np.linalg.inv(identity + f1)
        variance0 = float(target @ c0 @ target)
        variance1 = float(target @ c1 @ target)
        if not 0 <= variance1 <= variance0 + 1e-12 <= target @ target + 1e-10:
            raise RuntimeError("Gaussian information monotonicity failed")
        u0 = t.T @ ud
        u1 = t.T @ (ud + uv)
        mean0 = c0 @ u0
        mean1 = c1 @ u1
        # Integrated likelihood increment relative to zero shared response;
        # the unchanged per-object base residual constant cancels in comparisons.
        increment = 0.5 * (u1 @ mean1 - u0 @ mean0) - 0.5 * (
            np.linalg.slogdet(identity + f1)[1] - np.linalg.slogdet(identity + f0)[1]
        )
        rows.append(
            {
                "family": name,
                "dimensions": len(target),
                "prior_sd_mag": float(np.linalg.norm(target)),
                "discovery_sd_mag": np.sqrt(variance0),
                "combined_sd_mag": np.sqrt(variance1),
                "discovery_mean_response_mag": float(target @ mean0),
                "combined_mean_response_mag": float(target @ mean1),
                "validation_log_evidence_increment": float(increment),
            }
        )
    write_rows(out / "shared_mode_geometry.csv", rows)
    return {
        "calibration_modes": 12,
        "discovery_objects": 43,
        "validation_objects": 1020,
        "records": rows,
        "scope": "Exact Gaussian calculation for the bundled projected modes and priors. The matrices are derived inputs, not raw detector measurements. Shared-mode predictions and shrinking conditional variances do not identify an actual calibration bias or arbitrary grey luminosity evolution.",
    }
