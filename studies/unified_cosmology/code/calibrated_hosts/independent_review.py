"""Independent DESI band arithmetic, units, gradient and masked-continuum audit.

Reconstruct weights from observed-frame interval intersections. Test index
Jacobians against directional finite differences and fixed-anchor masking
against synthetic continua. No physical calibration or age validation follows.
"""

import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"
import sys, json, hashlib
from pathlib import Path
import numpy as np
import astropy.units as u
from bands import measure
from acquire import ROOT, WORK, RESULT

base = WORK
source = Path(__file__).with_name("bands.py")
rows = json.loads((base / "band-records.json").read_text())
rng = np.random.default_rng(927090)
windows = [
    (3850, 3950),
    (4000, 4100),
    (4041.6, 4079.75),
    (4083.5, 4122.25),
    (4128.5, 4161),
]
errors = []
derivatives = []
linear = []
constant = []
index_errors = []
for row in rows:
    with np.load(base / "spectra" / f"{row['targetid']}.npz") as ar:
        vs = []
        for camera, lo, hi in [
            ("B", -np.inf, 5780),
            ("R", 5780, 7570),
            ("Z", 7570, np.inf),
        ]:
            keep = (ar[camera + "_WAVELENGTH"] >= lo) & (
                ar[camera + "_WAVELENGTH"] < hi
            )
            vs.append(
                [
                    ar[camera + "_" + s][keep].astype(float)
                    for s in ["WAVELENGTH", "FLUX", "IVAR", "MASK"]
                ]
                + [np.repeat(camera, keep.sum())]
            )
        wave, flux, ivar, mask, arm = [
            np.concatenate([v[i] for v in vs]) for i in range(5)
        ]
    z = row["z"]
    good = np.isfinite(flux) & np.isfinite(ivar) & (ivar > 0) & (mask == 0)
    var = np.zeros(len(wave))
    var[good] = 1 / ivar[good]
    W = []
    widths = []
    for i, (lo, hi) in enumerate(windows):
        overlap = (
            np.clip(
                np.minimum(wave + 0.4, hi * (1 + z))
                - np.maximum(wave - 0.4, lo * (1 + z)),
                0,
                None,
            )
            * good
        )
        widths.append(overlap / (1 + z))
        weights = overlap / overlap.sum()
        if i < 2:
            weights *= 1e-17 * wave**2 / 2997924580000000000 / 1e-29
        W.append(weights)
    W = np.array(W)
    y = W @ np.where(good, flux, 0)
    C = (W * var) @ W.T
    se = np.sqrt(np.diag(C))
    saved = np.array(row["band_covariance_formal_diagonal_ivar"])
    errors.append(
        [
            float(np.max(abs(y - row["band_flux"]) / (1 + abs(y)))),
            float(np.max(abs(C - saved) / (se[:, None] * se[None, :]))),
        ]
    )
    got = measure(wave, flux, ivar, mask, arm, z)
    for k in ["Dn4000", "Hdelta_native_A"]:
        if got[k] is not None:
            index_errors.append(abs(got[k] - row[k]))
    for rep in range(3):
        delta = np.sqrt(var) * rng.normal(size=len(wave))
        h = 1e-5
        plus = measure(wave, flux + h * delta, ivar, mask, arm, z)
        minus = measure(wave, flux - h * delta, ivar, mask, arm, z)
        for k, g in zip(["Dn4000", "Hdelta_native_A"], got["gradient_arrays"]):
            if g is not None:
                numeric = (plus[k] - minus[k]) / (2 * h)
                pred = float(g @ delta)
                derivatives.append(abs(numeric - pred) / (1 + abs(pred)))
    const = measure(wave, np.ones(len(wave)), ivar, mask, arm, z)
    if const["Hdelta_native_A"] is not None:
        constant.append(abs(const["Hdelta_native_A"]))
    trend = 1 + (wave / (1 + z) - 4100) / 500
    tr = measure(wave, trend, ivar, mask, arm, z)
    if tr["Hdelta_native_A"] is not None:
        linear.append(
            {
                "targetid": row["targetid"],
                "spurious_EW_A": tr["Hdelta_native_A"],
                "actual_EW_sigma_A": row["Hdelta_formal_sigma_A"],
                "coverage": row["band_coverage"][2:],
            }
        )
unit = (1e-17 * u.erg / u.s / u.cm**2 / u.AA).to_value(
    u.uJy, equivalencies=u.spectral_density(4000 * u.AA)
)
assert max(x[0] for x in errors) < 1e-10 and max(x[1] for x in errors) < 1e-10
assert max(derivatives) < 1e-6 and max(index_errors) < 1e-10
r = {
    "status": "passed",
    "source_sha256": {
        str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in [Path(__file__), source, source.with_name("design.json")]
    },
    "input_sha256": {
        str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in [base / "band-records.json"]
        + [base / "spectra" / f"{row['targetid']}.npz" for row in rows]
    },
    "bands": len(rows),
    "maximum_mean_relative_error": max(x[0] for x in errors),
    "maximum_covariance_normalized_error": max(x[1] for x in errors),
    "gradient_direction_checks": len(derivatives),
    "maximum_gradient_relative_error": max(derivatives),
    "recomputed_indices_maximum_absolute_error": max(index_errors),
    "astropy_microJy_for_native_unit_at_4000A": unit,
    "coded_conversion": 1e12 * 4000**2 / 2.99792458e18,
    "constant_spectrum_maximum_spurious_EW_A": max(constant),
    "linear_spectrum_actual_masks": linear,
    "interpretation": "Formal diagonal-IVAR arithmetic review, not DESI uncertainty/calibration validation. Linear-continuum test slope0.2%/A is a diagnostic of fixed-anchor/masking response, not an empirical spectral systematic.",
}
(RESULT / "independent-review.json").write_text(json.dumps(r, indent=2) + "\n")
print(json.dumps({k: v for k, v in r.items() if k != "linear_spectrum_actual_masks"}))
print(
    "worstlinear",
    sorted(linear, key=lambda r: abs(r["spurious_EW_A"]), reverse=True)[:3],
)
