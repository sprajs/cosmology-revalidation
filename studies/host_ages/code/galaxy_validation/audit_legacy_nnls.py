#!/usr/bin/env python3
"""Read-only numerical audit of older BAO/RAISIN NNLS calls on real inputs.

Current first-party sources are loaded with only their old ROOT assignment
redirected to a user-supplied frozen input archive. No source or old result is
modified. The NNLS residual is checked against the returned coefficients and
against bounded-variable least squares on the same original design.
"""

from pathlib import Path
import argparse, datetime, json
import numpy as np, scipy
from scipy.optimize import nnls, lsq_linear
from scipy.linalg import solve_triangular
from astropy.cosmology import FlatLambdaCDM
from acquire import ROOT, OUT, sha
from match import ARCHIVE


def load(path, archive):
    source = path.read_text()
    old = "ROOT = Path(__file__).resolve().parents[2]"
    assert old in source
    source = source.replace(old, "ROOT = ARCHIVE_OVERRIDE", 1)
    ns = {
        "__file__": str(path),
        "__name__": "read_only_numerical_audit",
        "ARCHIVE_OVERRIDE": archive,
    }
    exec(compile(source, str(path), "exec"), ns)
    return ns


def check(a, y):
    x, r = nnls(a, y, maxiter=10000)
    actual = float(np.sum((a @ x - y) ** 2))
    reported = float(r * r)
    bv = lsq_linear(a, y, bounds=(0, np.inf), method="bvls", tol=1e-12, max_iter=3000)
    alt = float(np.sum((a @ bv.x - y) ** 2))
    norms = np.linalg.norm(a, axis=0)
    keep = norms > 0
    if keep.any():
        scaled = lsq_linear(
            a[:, keep] / norms[keep],
            y,
            bounds=(0, np.inf),
            method="bvls",
            tol=1e-13,
            max_iter=3000,
        )
        scaled_alt = float(np.sum(((a[:, keep] / norms[keep]) @ scaled.x - y) ** 2))
    else:
        scaled_alt = float(y @ y)
    gradient = a.T @ (a @ x - y)
    scale = max(np.linalg.norm(a) * np.linalg.norm(y), 1.0)
    kkt = (
        max(
            np.max(abs(gradient[x > 1e-10]), initial=0),
            np.max(-gradient[x <= 1e-10], initial=0),
        )
        / scale
    )
    return {
        "shape": list(a.shape),
        "reported_residual_squared": reported,
        "reconstructed_residual_squared": actual,
        "BVLS_residual_squared": alt,
        "column_scaled_BVLS_residual_squared": scaled_alt,
        "column_scaled_BVLS_objective_difference": abs(actual - scaled_alt),
        "column_norm_range": [float(norms.min()), float(norms.max())],
        "maximum_NNLS_weight": float(np.max(x)),
        "norm_identity_absolute_error": abs(reported - actual),
        "objective_difference_absolute": abs(actual - alt),
        "normalized_KKT_error": float(kkt),
        "BVLS_success": bool(bv.success),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", type=Path, default=ARCHIVE)
    args = ap.parse_args()
    paths = [
        ROOT / "studies/expansion/code/bao_shape.py",
        ROOT / "studies/infrared/code/raisin_differential.py",
    ]
    b = load(paths[0], args.archive)
    z, y, c, *_ = b["load_data"]()
    a, r, _ = b["cone_matrices"](np.log1p(z))
    cone = b["Cone"](c, r, a)
    rng = np.random.default_rng(260926)
    records = []
    for name, cc, yy in [
        ("full_cone", cone, y),
        (
            "radial_cone",
            b["Cone"](
                c[len(z) :, len(z) :],
                np.triu(np.ones((len(z), len(z)))),
                np.vstack([np.eye(len(z)), np.eye(len(z))[:-1] - np.eye(len(z))[1:]]),
            ),
            y[len(z) :],
        ),
    ]:
        tests = [("observed", yy)] + [
            (
                f"draw{k}",
                (np.zeros(len(yy)) if k < 500 else yy)
                + cc.L @ rng.normal(size=len(yy)),
            )
            for k in range(1000)
        ]
        for label, value in tests:
            v = solve_triangular(cc.L, value, lower=True)
            records.append({"model": name, "case": label, **check(cc.W, v)})
    r = load(paths[1], args.archive)
    raisin = []
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    for branch in r["BRANCHES"]:
        d = r["table"](branch)
        n = len(d)
        h = np.eye(n) - np.ones((n, n)) / n
        mu = cosmo.distmod(d.zHD.to_numpy()).value
        for group, ks in r["GROUPS"].items():
            c = r["syscov"](branch, group)
            target = (h @ c @ h).ravel()
            for flavor in ["DLMAG", "HD", "without_mass"]:
                vs = []
                for k in ks:
                    s = r["table"](branch, k).loc[d.index]
                    v = (s.DLMAG - d.DLMAG).to_numpy()
                    if flavor == "HD":
                        v -= cosmo.distmod(s.zHD.to_numpy()).value - mu
                    if flavor == "without_mass":
                        v -= (s.MASS_CORR - d.MASS_CORR).to_numpy()
                    vs.append(h @ v)
                x = np.stack([np.outer(v, v).ravel() for v in vs], axis=1)
                record = {
                    "branch": branch,
                    "group": group,
                    "flavor": flavor,
                    **check(x, target),
                }
                record["original_relative_covariance_residual"] = (
                    float(
                        np.sqrt(record["reconstructed_residual_squared"])
                        / np.linalg.norm(target)
                    )
                    if np.linalg.norm(target) > 0
                    else 0.0
                )
                if flavor == "without_mass" and any(k in [4, 5] for k in ks):
                    removed = [j for j, k in enumerate(ks) if k in [4, 5]]
                    assert max(np.max(abs(vs[j])) for j in removed) < 1e-12
                    keep = [j for j, k in enumerate(ks) if k not in [4, 5]]
                    xx = x[:, keep]
                    norms = np.linalg.norm(xx, axis=0)
                    nonzero = norms > 0
                    if np.any(nonzero):
                        scaled = xx[:, nonzero] / norms[nonzero]
                        fit = lsq_linear(
                            scaled,
                            target,
                            bounds=(0, np.inf),
                            method="bvls",
                            tol=1e-13,
                            max_iter=3000,
                        )
                        residual = scaled @ fit.x - target
                    else:
                        residual = -target
                    record["remove_algebraically_cancelled_massstep"] = {
                        "removed_variants": [ks[j] for j in removed],
                        "removed_max_abs_centered_magnitude_response": max(
                            float(np.max(abs(vs[j]))) for j in removed
                        ),
                        "remaining_variants": [ks[j] for j in keep],
                        "relative_covariance_residual": (
                            float(np.linalg.norm(residual) / np.linalg.norm(target))
                            if np.linalg.norm(target) > 0
                            else 0.0
                        ),
                        "interpretation": "Known massstep-only responses cancel in DLMAG-minus-MASS_CORR. Their ~1e-15mag floating residuals must be zero; no arbitrary coefficient may amplify them into physical covariance.",
                    }
                raisin.append(record)

    def summary(rows):
        return {
            "cases": len(rows),
            "shapes": sorted(set(tuple(x["shape"]) for x in rows)),
            "maximum_norm_identity_absolute_error": max(
                x["norm_identity_absolute_error"] for x in rows
            ),
            "maximum_BVLS_objective_difference_absolute": max(
                x["objective_difference_absolute"] for x in rows
            ),
            "maximum_normalized_KKT_error": max(
                x["normalized_KKT_error"] for x in rows
            ),
            "maximum_column_scaled_BVLS_objective_difference": max(
                x["column_scaled_BVLS_objective_difference"] for x in rows
            ),
            "all_BVLS_success": all(x["BVLS_success"] for x in rows),
        }

    result = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "scipy_version": scipy.__version__,
        "input_sources_sha256": {str(p.relative_to(ROOT)): sha(p) for p in paths},
        "audit_code_sha256": sha(Path(__file__)),
        "adaptation": "Only old source ROOT redirected in memory to --archive, no algorithms modified; no writes to archive or older results.",
        "BAO": summary(records),
        "BAO_observed": [x for x in records if x["case"] == "observed"],
        "RAISIN": summary(raisin),
        "RAISIN_records": raisin,
        "limits": "These actual designs and finite draws pass/fail only this audit; this is not a universal validation of scipy.optimize.nnls. BAO uses square normalized cone matrices; RAISIN overdetermined vectorized covariance designs and independently reconstructs residuals in original code.",
    }
    (OUT / "legacy-nnls-audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {k: result[k] for k in ["scipy_version", "BAO", "BAO_observed", "RAISIN"]},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
