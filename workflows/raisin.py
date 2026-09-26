"""Paired branch accounting, covariance closure and released-flux inventory."""

import numpy as np
from scipy.optimize import nnls
from lib.paths import DATA
from lib.records import fitres, write_rows
from lib.photometry import read_raisin

DEFAULTS = {}
BASE = DATA / "raisin"


def distance(branch, option=0):
    return fitres(
        BASE / f"distances/w/{branch}_dist/RAISIN_combined_FITOPT{option:03d}.FITRES"
    )


def covariance(branch, group):
    folder = BASE / f"distances/w/{branch}_syst"
    rows = np.loadtxt(folder / f"RAISIN_{group}_lcparams_cosmosis.txt")
    stat = np.loadtxt(folder / "RAISIN_stat_lcparams_cosmosis.txt")
    raw = np.loadtxt(folder / f"RAISIN_{group}.covmat")
    n = int(raw[0])
    if len(rows) != n or len(raw) != 1 + n * n:
        raise ValueError("Covariance dimensions do not match table")
    matrix = raw[1:].reshape(n, n) + np.diag(rows[:, 5] ** 2 - stat[:, 5] ** 2)
    if not np.allclose(matrix, matrix.T, atol=1e-12):
        raise ValueError("Asymmetric exported covariance")
    return matrix


def run(out, cfg):
    tabs = {name: distance(name) for name in ["nir", "optical", "opticalnir"]}
    nir = tabs["nir"]
    low = nir.zHD.to_numpy() < 0.1
    high = nir.zHD.to_numpy() > 0.2
    if len(nir) != 79 or low.sum() != 42 or high.sum() != 37:
        raise ValueError("Frozen 42-low/37-high common cohort changed")
    for name, tab in tabs.items():
        if not np.array_equal(tab.index, nir.index) or not np.array_equal(
            tab.zHD, nir.zHD
        ):
            raise ValueError("Branch identity/redshift order differs")
    weight = high / high.sum() - low / low.sum()
    contrasts = []
    closure = []
    variants = {
        "released": (0, 0),
        "without_mass": (0, -1),
        "without_bias": (1, 0),
        "without_bias_or_mass": (1, -1),
    }
    for branch in ["optical", "opticalnir"]:
        tab = tabs[branch]
        for label, (bias, mass) in variants.items():
            delta = (
                tab.DLMAG
                - nir.DLMAG
                + bias * (tab.DLMAG_biascor - nir.DLMAG_biascor)
                + mass * (tab.MASS_CORR - nir.MASS_CORR)
            ).to_numpy()
            se = np.sqrt(
                np.var(delta[high], ddof=1) / high.sum()
                + np.var(delta[low], ddof=1) / low.sum()
            )
            a, b = tab.DLMAGERR.to_numpy(), nir.DLMAGERR.to_numpy()
            contrasts.append(
                {
                    "branch_minus_nir": branch,
                    "variant": label,
                    "high_minus_low_mag": float(weight @ delta),
                    "descriptive_paired_se_mag": float(se),
                    "unknown_pair_correlation_min_stat_se": float(
                        np.sqrt(np.sum(weight**2 * (a - b) ** 2))
                    ),
                    "unknown_pair_correlation_max_stat_se": float(
                        np.sqrt(np.sum(weight**2 * (a + b) ** 2))
                    ),
                }
            )
        C = covariance(branch, "all")
        Cn = covariance("nir", "all")
        a, b = np.sqrt(weight @ C @ weight), np.sqrt(weight @ Cn @ weight)
        closure.append(
            {
                "branch": branch,
                "unknown_cross_systematic_min_se": float(abs(a - b)),
                "unknown_cross_systematic_max_se": float(a + b),
            }
        )
    write_rows(out / "paired_contrasts.csv", contrasts)
    # Identity accounting: do not merge alternate photometry by brightness.
    by_cid = {}
    alternates = []
    photrows = []
    for path in sorted((BASE / "photometry/RAISIN").rglob("*")):
        if not path.is_file() or path.suffix.lower() != ".dat":
            continue
        header, rows = read_raisin(path)
        cid = header["SNID"].split()[0]
        if cid in by_cid:
            if cid in nir.index:
                raise ValueError("Ambiguous analysis-cohort photometry")
            alternates.append(str(path.relative_to(DATA)))
            continue
        by_cid[cid] = header
        flux = np.array([float(r["FLUXCAL"]) for r in rows])
        error = np.array([float(r["FLUXCALERR"]) for r in rows])
        photrows.append(
            {
                "CID": cid,
                "survey": header["SURVEY"].split()[0],
                "rows": len(rows),
                "negative": int((flux < 0).sum()),
                "zero": int((flux == 0).sum()),
                "nonfinite": int((~np.isfinite(flux)).sum()),
                "nonpositive_error": int((error <= 0).sum()),
            }
        )
    write_rows(out / "photometry_inventory.csv", photrows)
    nir_header_gap = max(
        abs(float(by_cid[c]["PEAKMJD"].split()[0]) - nir.loc[c, "PKMJD"])
        for c in nir.index
    )
    # Test whether nonnegative sums of the released signed variant vectors
    # reconstruct the covariance after eliminating a common magnitude offset.
    reconstruct = []
    H = np.eye(len(nir)) - np.ones((len(nir), len(nir))) / len(nir)
    for branch in tabs:
        baseline = tabs[branch]
        target = H @ covariance(branch, "all") @ H
        vectors = [
            H
            @ (
                distance(branch, k).loc[baseline.index].DLMAG - baseline.DLMAG
            ).to_numpy()
            for k in range(1, 29)
        ]
        X = np.stack([np.outer(v, v).ravel() for v in vectors], axis=1)
        coefficients, _ = nnls(X, target.ravel(), maxiter=10000)
        gap = float(
            np.linalg.norm(X @ coefficients - target.ravel()) / np.linalg.norm(target)
        )
        reconstruct.append(
            {
                "branch": branch,
                "relative_covariance_shape_residual": gap,
                "nonnegative_weights": coefficients.tolist(),
            }
        )
    return {
        "objects": 79,
        "low": 42,
        "high": 37,
        "contrasts": contrasts,
        "systematic_bounds": closure,
        "photometry_totals": {
            k: sum(r[k] for r in photrows)
            for k in ["rows", "negative", "zero", "nonfinite", "nonpositive_error"]
        },
        "alternate_files_outside_distance_sample": alternates,
        "nir_peak_header_max_gap_days": float(nir_header_gap),
        "covariance_reconstruction": reconstruct,
        "scope": "Same-object descriptive distance differences. Undoing additive exported terms leaves fitted light curves and selection fixed. Cross-branch covariance is unavailable; no independent combination, repaired covariance or cosmological correction is inferred.",
    }
