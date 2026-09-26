"""Fixed detector support shared by raw and calibrated dark comparisons."""

import json
import numpy as np
from lib.paths import DATA
from lib.dark_geometry import map_bpixtab, pixel_circle_fraction

BASE = DATA / "hst-dark"
QUADRANTS = {
    "B": (slice(5, 512), slice(5, 512)),
    "C": (slice(5, 512), slice(512, 1019)),
    "A": (slice(512, 1019), slice(5, 512)),
    "D": (slice(512, 1019), slice(512, 1019)),
}


def design():
    return json.loads((BASE / "design.json").read_text())


def support():
    """Rebuild eligibility from reference tables without opening outcome pixels."""
    bad = np.zeros((1024, 1024), bool)
    for name in ("3562029fi_bpx.fits", "3562028ni_bpx.fits"):
        bits, _ = map_bpixtab(BASE / "references" / name)
        bad |= bits != 0
    active = np.zeros_like(bad)
    active[5:1019, 5:1019] = True
    active &= ~bad
    if active.sum() != 993750:
        raise ValueError("Reference mask support changed")
    r0, r1, r2 = design()["spatial_operator"]["radii_pixels"]
    extent = int(np.ceil(r2 + 0.5))
    offsets = range(-extent, extent + 1)
    aperture, inner, outer = [
        np.array([[pixel_circle_fraction(x, y, r) for x in offsets] for y in offsets])
        for r in (r0, r1, r2)
    ]
    for weights, radius in zip((aperture, inner, outer), (r0, r1, r2)):
        if abs(weights.sum() - np.pi * radius**2) > 1e-8:
            raise ValueError("Circle area closure failed")
    annulus = outer - inner
    sites, coverage = [], []
    for y in range(32, 1024, 64):
        for x in range(32, 1024, 64):
            yy, xx = (
                slice(y - extent, y + extent + 1),
                slice(x - extent, x + extent + 1),
            )
            a, b = aperture * ~bad[yy, xx], annulus * ~bad[yy, xx]
            ac, bc = float(a.sum() / aperture.sum()), float(b.sum() / annulus.sum())
            row = {
                "x": x,
                "y": y,
                "aperture_coverage": ac,
                "annulus_coverage": bc,
                "eligible": ac >= 0.9 and bc >= 0.75,
            }
            coverage.append(row)
            if row["eligible"]:
                w = a - b * a.sum() / b.sum()
                if abs(w.sum()) > 1e-8:
                    raise ValueError("Constant-image cancellation failed")
                sites.append({**row, "yy": yy, "xx": xx, "weight": w})
    if len(sites) != 247:
        raise ValueError("Fixed aperture support changed")
    return active, sites, coverage


def blocks():
    for by in range(16):
        for bx in range(16):
            yield (
                by,
                bx,
                slice(max(5, 64 * by), min(1019, 64 * (by + 1))),
                slice(max(5, 64 * bx), min(1019, 64 * (bx + 1))),
            )


def time_weights():
    spec = design()["time_operator"]
    t, powers = np.array(spec["t_seconds_nominal"]), spec["powers"]
    h = []
    for power in powers:
        u = np.abs(np.arange(-3, 4) / 3.0) ** power
        centered = t - (u @ t) / u.sum()
        row = u * centered / (u @ centered**2)
        if abs(row.sum()) > 1e-16 or abs(row @ t - 1) > 1e-12:
            raise ValueError("Slope operator moment identities failed")
        h.append(row)
    return t, powers, np.array(h)
