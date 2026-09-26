"""Compare frozen native dark slopes and quoted variance on reference-only masks."""

import numpy as np
from astropy.io import fits
from lib.dark import BASE, QUADRANTS, support, blocks
from lib.records import write_rows

DEFAULTS = {}
PAIR = ("idbx41onq", "idbx43p7q")


def run(out, cfg):
    active, sites, coverage = support()
    arrays = []
    ledger = []
    for root in PAIR:
        with fits.open(BASE / "calibrated" / (root + "_flt.fits"), memmap=False) as f:
            if f[0].header["NSAMP"] != 8 or f[0].header["EXPFLAG"] != "NORMAL":
                raise ValueError("Unexpected calibrated input")
            a = {k: f[k, 1].data.copy() for k in ("SCI", "ERR", "DQ", "SAMP", "TIME")}
            if (
                any(v.shape != (1014, 1014) for v in a.values())
                or f["SCI", 1].header["BUNIT"] != "COUNTS/S"
                or f["ERR", 1].header["BUNIT"] != "COUNTS/S"
            ):
                raise ValueError("Calibrated units/geometry differ")
            for key in ("SCI", "ERR", "TIME"):
                if not np.isfinite(a[key][active[5:1019, 5:1019]]).all():
                    raise ValueError("Nonfinite value on fixed support")
            if np.any(a["ERR"][active[5:1019, 5:1019]] < 0):
                raise ValueError("Negative quoted error")
            # Translate trimmed FLTs to the zero-based RAW detector coordinates.
            arrays.append({k: np.pad(v, 5) for k, v in a.items()})
            for bit in (1 << i for i in range(16)):
                ledger.append(
                    {
                        "root": root,
                        "bit": bit,
                        "pixels": int(
                            np.sum(
                                (arrays[-1]["DQ"][active].astype(np.uint16) & bit) != 0
                            )
                        ),
                    }
                )
    old, new = arrays
    d = new["SCI"].astype(float) - old["SCI"].astype(float)
    v = old["ERR"].astype(float) ** 2 + new["ERR"].astype(float) ** 2
    if not np.all(v[active] > 0):
        raise ValueError("Nonpositive pair variance")

    def summarize(mask, label):
        x, q = d[mask], v[mask]
        ss = float(x @ x)
        vv = float(q.sum())
        return {
            "region": label,
            "pixels": len(x),
            "signed_mean": float(x.mean()),
            "sum_squared_difference": ss,
            "sum_quoted_pair_variance": vv,
            "repeat_to_quoted_ratio": ss / vv,
        }

    regions = [summarize(active, "all_active_fixed_mask")]
    for name, (yy, xx) in QUADRANTS.items():
        mask = np.zeros_like(active)
        mask[yy, xx] = active[yy, xx]
        regions.append(summarize(mask, name))
    blockrows = []
    for by, bx, yy, xx in blocks():
        mask = np.zeros_like(active)
        mask[yy, xx] = active[yy, xx]
        blockrows.append(summarize(mask, f"{by},{bx}"))
    apertures = []
    all_weight = np.zeros_like(d)
    for s in sites:
        yy, xx, w = s["yy"], s["xx"], s["weight"]
        value = float(np.sum(w * d[yy, xx]))
        variance = float(np.sum(w * w * v[yy, xx]))
        apertures.append(
            {
                "raw_x": s["x"],
                "raw_y": s["y"],
                "difference_dn_per_nominal_s": value,
                "quoted_diagonal_pair_variance": variance,
                "squared_difference": value**2,
            }
        )
        all_weight[yy, xx] += np.abs(w)
    # Descriptive influence decomposition; flags/ranks never change eligibility.
    flat = np.flatnonzero(active)
    ranking = flat[np.argsort(d.ravel()[flat] ** 2)[-10:][::-1]]
    total = regions[0]["sum_squared_difference"]
    influence = []
    for pos in ranking:
        y, x = np.unravel_index(pos, active.shape)
        influence.append(
            {
                "raw_x": int(x),
                "raw_y": int(y),
                "squared_difference": float(d[y, x] ** 2),
                "fraction_of_total": float(d[y, x] ** 2 / total),
                "absolute_aperture_weight": float(all_weight[y, x]),
                "earlier_DQ": int(old["DQ"][y, x]),
                "later_DQ": int(new["DQ"][y, x]),
            }
        )
    bit32 = active & (
        ((old["DQ"].astype(np.uint16) | new["DQ"].astype(np.uint16)) & 32) != 0
    )
    for name, rows in (
        ("regions", regions),
        ("blocks", blockrows),
        ("apertures", apertures),
        ("influence", influence),
        ("dq-bits", ledger),
        ("coverage", coverage),
    ):
        write_rows(out / (name + ".csv"), rows)
    return {
        "pair": list(PAIR),
        "regions": regions,
        "apertures": {
            "n": len(apertures),
            "ratio": sum(r["squared_difference"] for r in apertures)
            / sum(r["quoted_diagonal_pair_variance"] for r in apertures),
        },
        "largest_pixel": influence[0],
        "DQ32": {
            "pixels": int(bit32.sum()),
            "sum_squared_difference": float(np.sum(d[bit32] ** 2)),
        },
        "starting_point": "Frozen local CALWF3 3.7.3 FLTs; native calibration is not rerun by this workflow.",
        "scope": "Single conditional NORMAL dark pair, DN per nominal second. No output-DQ or residual cuts. Pixel and aperture ratios are distinct spatial estimands; no detector-noise attribution, historical photometry reconstruction or error rescaling.",
    }
