"""Extract fixed paired read moments and signed slopes from original dark RAWs."""

import numpy as np
from astropy.io import fits
from lib.dark import BASE, QUADRANTS, design, support, blocks, time_weights
from lib.records import write_rows

DEFAULTS = {}


def read_pair(earlier, later, times):
    delta = np.empty((7, 1024, 1024), dtype=np.float64)
    endpoints, flagged = np.zeros((1024, 1024), bool), np.zeros((1024, 1024), bool)
    with (
        fits.open(
            BASE / "raw" / (earlier + "_raw.fits"), memmap=False, uint=True
        ) as old,
        fits.open(BASE / "raw" / (later + "_raw.fits"), memmap=False, uint=True) as new,
    ):
        for f in (old, new):
            if f[0].header["NSAMP"] != 16 or f[0].header["SAMP_SEQ"] != "SPARS50":
                raise ValueError("Unexpected RAW sampling")
        for i, extver in enumerate(range(15, 8, -1)):
            a, b = old["SCI", extver].data, new["SCI", extver].data
            if (
                a.dtype != np.uint16
                or b.dtype != np.uint16
                or a.shape != (1024, 1024)
                or b.shape != a.shape
            ):
                raise ValueError("Physical unsigned 16-bit RAW encoding required")
            # Promote before subtracting: uint16 subtraction can wrap.
            delta[i] = b.astype(np.float64) - a.astype(np.float64)
            endpoints |= (a == 0) | (a == 65535) | (b == 0) | (b == 65535)
            for f in (old, new):
                if abs(f["SCI", extver].header["SAMPTIME"] - times[i]) > 1e-6:
                    raise ValueError("Nominal sampling time differs")
                dq = f["DQ", extver]
                flagged |= (
                    bool(dq.header["PIXVALUE"])
                    if dq.header["NAXIS"] == 0
                    else dq.data != 0
                )
    return delta, endpoints, flagged


def moments(v, h):
    n = v.shape[1]
    if not n:
        raise ValueError("Empty fixed region")
    mean = v.mean(axis=1)
    gamma = v @ v.T / (2 * n)
    projected = h @ v
    direct = np.mean(projected**2, axis=1) / 2
    matrix = np.einsum("pi,ij,pj->p", h, gamma, h)
    if not np.allclose(direct, matrix, rtol=1e-10, atol=1e-12):
        raise ValueError("Direct/matrix contraction failed")
    return mean, gamma, direct, float(np.max(abs(direct - matrix)))


def run(out, cfg):
    active, sites, coverage = support()
    times, powers, h = time_weights()
    plan = design()
    with fits.open(BASE / "references/t2c16200i_ccd.fits") as f:
        matches = [
            r
            for r in f["CCD"].data
            if str(r["CCDAMP"]).strip() == "ABCD"
            and r["CCDCHIP"] == 1
            and r["CCDGAIN"] == 2.5
            and r["BINAXIS1"] == r["BINAXIS2"] == 1
        ]
        if len(matches) != 1 or matches[0]["AMPX"] != 512 or matches[0]["AMPY"] != 512:
            raise ValueError("Ambiguous gain reference")
        gains = {q: float(matches[0]["ATODGN" + q]) for q in QUADRANTS}
    qg = np.empty((8, 2, 4, 7, 7))
    qm = np.empty((8, 2, 4, 7))
    bg = np.empty((8, 16, 16, 7, 7))
    bm = np.empty((8, 16, 16, 7))
    quadrants, blockrows, apertures, categories, decomposition = [], [], [], [], []
    max_gap = 0.0
    for ip, pair in enumerate(plan["pairs"]):
        d, endpoints, flagged = read_pair(pair["earlier"], pair["later"], times)
        for iq, (q, (yy, xx)) in enumerate(QUADRANTS.items()):
            for j, (scope, mask) in enumerate(
                (
                    ("masked", active[yy, xx]),
                    ("geometry_only", np.ones(active[yy, xx].shape, bool)),
                )
            ):
                mean, gamma, scores, gap = moments(d[:, yy, xx][:, mask], h)
                max_gap = max(max_gap, gap)
                qm[ip, j, iq] = mean
                qg[ip, j, iq] = gamma
                for power, score in zip(powers, scores):
                    quadrants.append(
                        {
                            "pair": pair["pair"],
                            "kind": pair["kind"],
                            "quadrant": q,
                            "scope": scope,
                            "n": int(mask.sum()),
                            "power": power,
                            "dn2_per_nominal_s2": float(score),
                            "electron2_per_nominal_s2": float(score * gains[q] ** 2),
                        }
                    )
                if j == 0:
                    coherent = float((h[0] @ mean) ** 2 / 2)
                    decomposition.append(
                        {
                            "pair": pair["pair"],
                            "quadrant": q,
                            "total": float(scores[0]),
                            "coherent": coherent,
                            "centered": float(scores[0] - coherent),
                            "coherent_fraction": coherent / float(scores[0]),
                        }
                    )
        for by, bx, yy, xx in blocks():
            mean, gamma, scores, gap = moments(d[:, yy, xx][:, active[yy, xx]], h)
            max_gap = max(max_gap, gap)
            bm[ip, by, bx] = mean
            bg[ip, by, bx] = gamma
            blockrows.append(
                {
                    "pair": pair["pair"],
                    "block_y": by,
                    "block_x": bx,
                    "n": int(active[yy, xx].sum()),
                    "primary_dn2_per_nominal_s2": float(scores[0]),
                }
            )
        for s in sites:
            z = np.einsum("tyx,yx->t", d[:, s["yy"], s["xx"]], s["weight"])
            for power, slope in zip(powers, h @ z):
                apertures.append(
                    {
                        "pair": pair["pair"],
                        "x": s["x"],
                        "y": s["y"],
                        "power": power,
                        "signed_dn_per_nominal_s": float(slope),
                        "half_squared_dn2_per_nominal_s2": float(slope**2 / 2),
                    }
                )
        primary = np.einsum("t,tyx->yx", h[0], d)
        for label, condition in (("digital_endpoint", endpoints), ("raw_DQ", flagged)):
            categories.append(
                {
                    "pair": pair["pair"],
                    "category": label,
                    "selected_pixels": int((active & condition).sum()),
                    "half_squared_slope_contribution": float(
                        np.sum(primary[active & condition] ** 2) / 2
                    ),
                }
            )
        del d
    np.savez_compressed(
        out / "read-matrix.npz",
        quadrant_gamma=qg,
        quadrant_mean=qm,
        block_gamma=bg,
        block_mean=bm,
        quadrant_names=np.array(list(QUADRANTS)),
        powers=np.array(powers),
    )
    for name, rows in (
        ("coverage", coverage),
        ("quadrants", quadrants),
        ("blocks", blockrows),
        ("apertures", apertures),
        ("categories", categories),
        ("decomposition", decomposition),
        ("cohort", plan["cohort"]),
    ):
        write_rows(out / (name + ".csv"), rows)
    return {
        "raw_exposures": len(plan["cohort"]),
        "primary_pairs": 7,
        "overlapping_sensitivity_pairs": 1,
        "fixed_pixels": int(active.sum()),
        "eligible_apertures": len(sites),
        "powers": powers,
        "indeterminate_exposures": sum(
            r["EXPFLAG"] == "INDETERMINATE" for r in plan["cohort"]
        ),
        "max_direct_matrix_gap": max_gap,
        "gains": gains,
        "coherent_fraction_range": [
            min(r["coherent_fraction"] for r in decomposition),
            max(r["coherent_fraction"] for r in decomposition),
        ],
        "scope": "Total digitized RAW repeat moments in nominal time coordinates. Unflagged events, shot noise and detector-state differences remain mixed; no electronic-noise isolation, independent-pixel significance, or error rescaling.",
    }
