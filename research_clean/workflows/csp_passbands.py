"""Frozen spectral/passband photon integrals with independent quadrature."""

import numpy as np
from astropy.io import fits
from lib.paths import DATA
from lib.passbands import counts
from lib.records import write_rows

DEFAULTS = {
    "phase_days": [-7, 0, 7, 10, 15, 20, 30],
    "redshifts": [0, 0.01, 0.03, 0.05, 0.08],
}
PAIRS = {
    "J_WIRC_minus_RC1": ("Jrc1_SWO_TAM_scan_atm.dat", "J_DUP_TAM_scan_atm.dat"),
    "H_WIRC_minus_RetroCam": ("H_SWO_TAM_scan_atm.dat", "H_DUP_TAM_scan_atm.dat"),
    "Y_WIRC_minus_RetroCam": ("Y_SWO_TAM_scan_atm.dat", "Y_DUP_TAM_scan_atm.dat"),
    "J_RC2_minus_RC1": ("Jrc1_SWO_TAM_scan_atm.dat", "Jrc2_SWO_TAM_scan_atm.dat"),
}


def run(out, cfg):
    if (
        0 not in cfg["phase_days"]
        or 0 not in cfg["redshifts"]
        or any(z < 0 for z in cfg["redshifts"])
    ):
        raise ValueError(
            "Grid must include phase zero and redshift zero, with nonnegative redshifts"
        )
    with fits.open(DATA / "raisin/kcor/kcor_CSPDR3_BD17.fits") as f:
        h = f["SN SED"].header
        wave = h["LMIN"] + np.arange(h["NBL"]) * h["LBIN"]
        times = h["TMIN"] + np.arange(h["NBT"]) * h["TBIN"]
        spectra = np.asarray(f["SN SED"].data.field(0), float).reshape(
            len(times), len(wave)
        )
        rw = np.asarray(f["PrimarySED"].data.field(0), float)
        rf = np.asarray(f["PrimarySED"].data["BD17"], float)
    bands = {
        n: np.loadtxt(DATA / "csp/filters" / n) for pair in PAIRS.values() for n in pair
    }
    refs = {
        order: {n: counts(rw, rf, b, order) for n, b in bands.items()}
        for order in [2, 4, 0]
    }
    rows = []
    for phase in cfg["phase_days"]:
        ix = np.flatnonzero(times == phase)
        if len(ix) != 1:
            raise ValueError("Phase not present in archived SED")
        for z in cfg["redshifts"]:
            sc = {
                order: {
                    n: counts(wave * (1 + z), spectra[ix[0]] / (1 + z), b, order)
                    for n, b in bands.items()
                }
                for order in [2, 4, 0]
            }
            if not all(
                np.isfinite(v) and v > 0
                for group in sc.values()
                for v in group.values()
            ):
                raise ValueError("Nonpositive photon integral")
            for label, (a, b) in PAIRS.items():
                dm = {
                    o: float(
                        -2.5
                        * np.log10((sc[o][b] / refs[o][b]) / (sc[o][a] / refs[o][a]))
                    )
                    for o in [2, 4, 0]
                }
                rows.append(
                    {
                        "pair": label,
                        "phase": phase,
                        "z": z,
                        "delta_mag_equal_BD17": dm[2],
                        "gauss2_4_gap": abs(dm[2] - dm[4]),
                        "trapezoid_halfA_gap": abs(dm[2] - dm[0]),
                    }
                )
    for r in rows:
        phase0 = next(
            x
            for x in rows
            if x["pair"] == r["pair"] and x["z"] == r["z"] and x["phase"] == 0
        )
        z0 = next(
            x
            for x in rows
            if x["pair"] == r["pair"] and x["phase"] == r["phase"] and x["z"] == 0
        )
        r["phase_contrast_vs_phase0"] = (
            r["delta_mag_equal_BD17"] - phase0["delta_mag_equal_BD17"]
        )
        r["redshift_contrast_vs_z0"] = (
            r["delta_mag_equal_BD17"] - z0["delta_mag_equal_BD17"]
        )
    g = max(x["gauss2_4_gap"] for x in rows)
    t = max(x["trapezoid_halfA_gap"] for x in rows)
    if g > 1e-10 or t > 1e-4:
        raise RuntimeError("Independent integration gate failed")
    write_rows(out / "passband_grid.csv", rows)
    return {
        "grid_points": len(rows),
        "max_gauss2_4_gap": g,
        "max_halfA_trapezoid_gap": t,
        "scope": "Conditional passband-shape sensitivity for the archived KCOR SN SED and BD17 reference. No fitted zero-point offset, observed brightness calibration or cosmological correction. Signed transmission tails are retained.",
    }
