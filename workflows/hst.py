"""Repeat signed photometry directly from calibrated pixels and fixed operators."""

import json
import numpy as np
from astropy.io import fits
from lib.paths import DATA
from lib.records import write_rows, sha256

DEFAULTS = {"bootstrap_draws": 2000, "seed": 260926}


def run(out, cfg):
    if cfg["bootstrap_draws"] < 1:
        raise ValueError("bootstrap_draws must be positive")
    state = json.loads((DATA / "hst/stage_a.json").read_text())
    if (
        sha256(DATA / "hst/operators.npz") != state["operators_sha256"]
        or sha256(DATA / "hst/source_masks.npz") != state["source_masks_sha256"]
    ):
        raise ValueError("Frozen geometry file identity mismatch")
    with np.load(DATA / "hst/operators.npz", allow_pickle=False) as f:
        ops = {k: f[k] for k in ["flat", "aperture", "annulus", "pam"]}
    offsets = {
        (r["candidate"], r["exposure"]): (r["start"], r["end"])
        for r in state["operator_rows"]
    }
    strict = [r for r in state["rows"] if r.get("strict")]
    secondary = [r for r in state["rows"] if r.get("secondary")]
    if len(secondary) < 30:
        raise RuntimeError("Insufficient masked secondary support")
    exposures = state["exposures"]
    ledger = []
    measurements = {}
    for i, e in enumerate(exposures):
        if e["ordinal"] not in [2, 4]:
            continue
        path = DATA / "hst/exposures" / e["filename"]
        if sha256(path) != e["sha256"]:
            raise ValueError("Exposure differs from frozen geometry")
        with fits.open(path, memmap=True, do_not_scale_image_data=True) as f:
            for candidate in secondary:
                lo, hi = offsets[candidate["id"], i]
                flat = ops["flat"][lo:hi]
                aperture, annulus, pam = [
                    ops[k][lo:hi].astype(float) for k in ["aperture", "annulus", "pam"]
                ]
                y, error = [
                    f[k].data.ravel()[flat].astype(float) for k in ["SCI", "ERR"]
                ]
                dq = f["DQ"].data.ravel()[flat]
                valid = np.isfinite(y) & np.isfinite(error) & (error > 0) & (dq == 0)
                acov = np.sum(aperture[valid] * pam[valid]) / np.sum(aperture * pam)
                bcov = np.sum(annulus[valid]) / np.sum(annulus)
                if acov < 0.90 - 1e-7 or bcov < 0.75 - 1e-7:
                    raise RuntimeError("Frozen masked coverage changed")
                a, b = aperture[valid] * pam[valid], annulus[valid]
                k = e["photflam"] / exposures[0]["photflam"]
                weight = k * (a - a.sum() * b / b.sum())
                flux = float(weight @ y[valid])
                variance = float(np.sum((weight * error[valid]) ** 2))
                if abs(weight.sum()) > 1e-8 or not variance > 0:
                    raise RuntimeError(
                        "Aperture constant-image or positive variance gate failed"
                    )
                row = {
                    "candidate": candidate["id"],
                    "tile": candidate["tile"],
                    "visit": e["visit"],
                    "ordinal": e["ordinal"],
                    "filename": e["filename"],
                    "flux": flux,
                    "variance": variance,
                    "aperture_coverage": acov,
                    "annulus_coverage": bcov,
                }
                ledger.append(row)
                measurements[candidate["id"], e["visit"], e["ordinal"]] = row
    write_rows(out / "signed_exposures.csv", ledger)
    pairs = []
    summaries = {}
    rng = np.random.default_rng(cfg["seed"])
    for visit in ["search", "template"]:
        rows = []
        for c in secondary:
            a, b = [measurements[c["id"], visit, j] for j in [2, 4]]
            rows.append(
                {
                    "candidate": c["id"],
                    "tile": c["tile"],
                    "visit": visit,
                    "difference": a["flux"] - b["flux"],
                    "variance": a["variance"] + b["variance"],
                }
            )
        pairs.extend(rows)
        d = np.array([r["difference"] for r in rows])
        v = np.array([r["variance"] for r in rows])
        tiles = np.array([r["tile"] for r in rows])

        def metric(ix):
            dd, vv = d[ix], v[ix]
            intercept = np.sum(dd / vv) / np.sum(1 / vv)
            return float(np.sum((dd - intercept) ** 2 / vv) / (len(ix) - 1))

        unique = np.unique(tiles)
        counts = [int((tiles == t).sum()) for t in unique]
        boot = []
        if len(unique) >= 12 and min(counts) >= 3:
            for _ in range(cfg["bootstrap_draws"]):
                draw = rng.choice(unique, len(unique), replace=True)
                ix = np.concatenate([np.flatnonzero(tiles == t) for t in draw])
                boot.append(metric(ix))
        summaries[visit] = {
            "n": len(d),
            "tiles": len(unique),
            "sum_d2_over_sum_V": float((d @ d) / v.sum()),
            "centered_variance_ratio": metric(np.arange(len(d))),
            "signed_standardized_mean": float(np.mean(d / np.sqrt(v))),
            "tile_bootstrap_95_interval": np.quantile(boot, [0.025, 0.975]).tolist()
            if boot
            else None,
            "tile_bootstrap_draws": len(boot),
        }
    write_rows(out / "signed_pairs.csv", pairs)
    return {
        "strict_support": len(strict),
        "strict_status": "no support; no estimate"
        if len(strict) < 30
        else "not scored by this secondary workflow",
        "masked_secondary_support": len(secondary),
        "visits": summaries,
        "scope": "Held-out signed repeat aperture sums on current calibrated HST products and frozen design-only operators. The variance uses quoted diagonal ERR and independent exposures as a reference. Correlated pixels, detector reads, reference terms, source masking and historical reduction remain separate questions; no errors are rescaled.",
    }
