"""Compute progenitor-delay grids with independent clock and resolution checks."""

import numpy as np
import pandas as pd
from astropy.cosmology import Flatw0waCDM
from lib.populations import COSMOS, DTDS, clock, calculate

DEFAULTS = {
    "cosmology": "son_cpl_H63p6",
    "sfh": "B13",
    "delay_points": 12001,
    "clock_points": 60001,
}


def run(out, cfg):
    c = COSMOS[cfg["cosmology"]]
    if (
        cfg["sfh"] not in {"B13", "MD14"}
        or cfg["delay_points"] < 1000
        or cfg["clock_points"] < 1000
    ):
        raise ValueError("Unknown star formation law or inadequate integration grid")
    ages, _ = clock(c, cfg["clock_points"])
    z = np.array([0.0, 0.1, 0.5, 1.0, 2.0, 10.0])
    error = float(
        np.max(abs(ages(-np.log1p(z)) - Flatw0waCDM(**c, Tcmb0=0).age(z).value))
    )
    frames = []
    for name in DTDS:
        frame = calculate(c, cfg["sfh"], name, cfg["delay_points"], cfg["clock_points"])
        frame["dtd"] = name
        frames.append(frame)
    table = pd.concat(frames, ignore_index=True)
    fine = calculate(
        c,
        cfg["sfh"],
        "C14_smooth",
        2 * cfg["delay_points"] - 1,
        2 * cfg["clock_points"] - 1,
    )
    coarse = frames[0]
    gap = max(
        float(np.max(abs(fine[k] - coarse[k])))
        for k in ["mean_delay_gyr", "median_delay_gyr", "formed_mass_mean_gyr"]
    )
    if error >= 2e-6 or gap >= 1e-4:
        raise RuntimeError(f"Clock/resolution gate failed: {error}, {gap}")
    table.to_csv(out / "delay_curves.csv", index=False)
    for label, column in [("mean", "mean_delay_gyr"), ("median", "median_delay_gyr")]:
        pd.DataFrame(
            {"z": coarse.z, "delta_mu": coarse["correction_subtracted_mag_" + column]}
        ).to_csv(out / (label + "_correction.csv"), index=False)
    return {
        "clock_max_error_gyr": error,
        "double_resolution_max_error_gyr": gap,
        "rows": len(table),
        "at_z1": table[np.isclose(table.z, 1)].to_dict("records"),
        "scope": "Cosmic-volume SFH/DTD convolution. Galaxy ages, progenitor delays and survey-selected SN populations are distinct. The fixed 0.030 mag/Gyr template coefficient is an imposed sensitivity, not an identified extra correction.",
    }
