"""Frozen tangent-plane and fractional-pixel operators for WFC3/IR FLTs.

This module reads FITS headers and the official PAM. It never reads a science
SCI/ERR/DQ plane. FITS pixel coordinates are one-based throughout.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from astropy.wcs import WCS

RAD = math.pi / 180
ASEC = 3600 / RAD


def tangent(ra_deg, dec_deg, ra0_deg, dec0_deg):
    ra = np.asarray(ra_deg) * RAD
    dec = np.asarray(dec_deg) * RAD
    ra0, dec0 = ra0_deg * RAD, dec0_deg * RAD
    d = ra - ra0
    den = math.sin(dec0) * np.sin(dec) + math.cos(dec0) * np.cos(dec) * np.cos(d)
    return (
        np.cos(dec) * np.sin(d) / den * ASEC,
        (math.cos(dec0) * np.sin(dec) - math.sin(dec0) * np.cos(dec) * np.cos(d))
        / den
        * ASEC,
    )


def untangent(east_arcsec, north_arcsec, ra0_deg, dec0_deg):
    xi, eta = np.asarray(east_arcsec) / ASEC, np.asarray(north_arcsec) / ASEC
    ra0, dec0 = ra0_deg * RAD, dec0_deg * RAD
    den = math.cos(dec0) - eta * math.sin(dec0)
    ra = ra0 + np.arctan2(xi, den)
    dec = np.arctan2(math.sin(dec0) + eta * math.cos(dec0), np.hypot(den, xi))
    return np.rad2deg(ra) % 360, np.rad2deg(dec)


def world_tangent(wcs: WCS, pix_xy, ra0, dec0):
    xy = np.asarray(pix_xy, float)
    rd = wcs.all_pix2world(xy, 1)
    e, n = tangent(rd[:, 0], rd[:, 1], ra0, dec0)
    return np.column_stack((e, n))


def curvature_gate(wcs: WCS, ra0: float, dec0: float, grid: list[int]) -> float:
    """Maximum direct-WCS minus corner-bilinear offset in native pixels."""
    uv = np.array([(u, v) for u in (0.25, 0.5, 0.75) for v in (0.25, 0.5, 0.75)], float)
    maxpix = 0.0
    for y in grid:
        for x in grid:
            cornerpix = np.array(
                [
                    (x + dx, y + dy)
                    for dx, dy in [(-0.5, -0.5), (0.5, -0.5), (-0.5, 0.5), (0.5, 0.5)]
                ]
            )
            c = world_tangent(wcs, cornerpix, ra0, dec0)
            a = uv[:, 0, None]
            b = uv[:, 1, None]
            bilinear = (
                c[0] * (1 - a) * (1 - b)
                + c[1] * a * (1 - b)
                + c[2] * (1 - a) * b
                + c[3] * a * b
            )
            directpix = np.column_stack((x - 0.5 + uv[:, 0], y - 0.5 + uv[:, 1]))
            direct = world_tangent(wcs, directpix, ra0, dec0)
            # Conservative scale: smallest native IR axis is roughly .12 arcsec/pixel.
            maxpix = max(
                maxpix, float(np.max(np.linalg.norm(direct - bilinear, axis=1) / 0.12))
            )
    return maxpix


@dataclass
class Operator:
    flat: np.ndarray
    aperture: np.ndarray
    annulus: np.ndarray
    pam: np.ndarray
    resolution: int
    closure: float


def _disk_overlap(corners: np.ndarray, radius: float, resolution: int) -> np.ndarray:
    """Fractional native-pixel overlap, using exact WCS pixel corners.

    Corners have shape (npixels,4,2), ordered lower-left, lower-right,
    upper-left, upper-right. Subpixel interpolation is independently bounded
    by curvature_gate before image values are used.
    """
    rr = np.linalg.norm(corners, axis=2)
    fully_in = rr.max(axis=1) <= radius
    xmin, xmax = corners[:, :, 0].min(axis=1), corners[:, :, 0].max(axis=1)
    ymin, ymax = corners[:, :, 1].min(axis=1), corners[:, :, 1].max(axis=1)
    nearest_x = np.maximum.reduce((xmin, np.zeros_like(xmin), -xmax))
    nearest_y = np.maximum.reduce((ymin, np.zeros_like(ymin), -ymax))
    fully_out = np.hypot(nearest_x, nearest_y) >= radius
    frac = fully_in.astype(float)
    edge = np.flatnonzero(~fully_in & ~fully_out)
    if len(edge):
        uv = (np.arange(resolution) + 0.5) / resolution
        u, v = np.meshgrid(uv, uv)
        u = u.ravel()[None, :, None]
        v = v.ravel()[None, :, None]
        # Bound peak transient allocation at the 128x128 refinement. Each
        # eight-pixel chunk's two-coordinate sample is 2 MiB at that level;
        # intermediate arithmetic arrays remain within the array budget.
        for start in range(0, len(edge), 8):
            part = edge[start : start + 8]
            c = corners[part]
            sample = (
                c[:, 0, None, :] * (1 - u) * (1 - v)
                + c[:, 1, None, :] * u * (1 - v)
                + c[:, 2, None, :] * (1 - u) * v
                + c[:, 3, None, :] * u * v
            )
            frac[part] = np.mean(
                np.sum(sample * sample, axis=2) < radius * radius, axis=1
            )
    return frac


def _one_resolution(corners, flat, pam, resolution, aperture_radius, annulus):
    a = _disk_overlap(corners, aperture_radius, resolution)
    inner = _disk_overlap(corners, annulus[0], resolution)
    outer = _disk_overlap(corners, annulus[1], resolution)
    b = np.maximum(0, outer - inner)
    use = (a > 0) | (b > 0)
    return flat[use], a[use], b[use], pam[use]


def make_operator(
    wcs: WCS,
    center_east: float,
    center_north: float,
    ra0: float,
    dec0: float,
    pam_image: np.ndarray,
    aperture_radius: float,
    annulus: tuple[float, float],
    levels=(16, 32, 64, 128),
) -> Operator:
    ra, dec = untangent(center_east, center_north, ra0, dec0)
    xy0 = wcs.all_world2pix([[ra, dec]], 1)[0]
    # 2 arcsec annulus plus conservative native-scale/boundary margin.
    bound = math.ceil(annulus[1] / 0.11) + 4
    xx = np.arange(
        max(1, math.floor(xy0[0]) - bound), min(1014, math.ceil(xy0[0]) + bound) + 1
    )
    yy = np.arange(
        max(1, math.floor(xy0[1]) - bound), min(1014, math.ceil(xy0[1]) + bound) + 1
    )
    xg, yg = np.meshgrid(xx, yy)
    x, y = xg.ravel(), yg.ravel()
    assert len(x) > 0
    xy = np.column_stack((x, y))
    c = []
    for dx, dy in [(-0.5, -0.5), (0.5, -0.5), (-0.5, 0.5), (0.5, 0.5)]:
        c.append(
            world_tangent(wcs, xy + [dx, dy], ra0, dec0) - [center_east, center_north]
        )
    corners = np.stack(c, axis=1)
    flat = ((y - 1) * 1014 + x - 1).astype(np.int32)
    p = pam_image[y - 1, x - 1].astype(float)
    if not np.all(np.isfinite(p)) or not np.all(p > 0):
        raise ValueError("invalid PAM in operator footprint")
    prior = None
    for resolution in levels:
        now = _one_resolution(corners, flat, p, resolution, aperture_radius, annulus)
        if prior is not None:
            # Compare true fixed full-length vectors, including boundary pixels
            # that enter/leave with the doubled integration grid.
            def dense(o):
                aa = np.zeros(len(flat))
                bb = np.zeros(len(flat))
                loc = np.searchsorted(flat, o[0]) if np.all(np.diff(flat) > 0) else None
                if loc is None:
                    look = {int(f): i for i, f in enumerate(flat)}
                    loc = [look[int(f)] for f in o[0]]
                aa[loc], bb[loc] = o[1], o[2]
                return aa, bb

            pa, pb = dense(prior)
            na, nb = dense(now)
            # Synthetic signed scene and heteroscedastic ERR fixed in detector
            # coordinates; compares both flux and assembled-weight variance.
            mid = corners.mean(axis=1)
            synthetic = 1 + np.exp(-0.5 * np.sum(mid * mid, axis=1) / (0.3**2))
            err = 1 + 0.0002 * x + 0.0003 * y

            def score(a, b):
                ap = a * p
                weight = ap - ap.sum() * b / b.sum()
                return (
                    np.dot(weight, synthetic),
                    np.dot(weight * weight, err * err),
                    ap.sum(),
                    b.sum(),
                )

            old, new = score(pa, pb), score(na, nb)
            relative = [abs(q - z) / max(abs(z), 1e-12) for q, z in zip(old, new)]
            closure = max(relative)
            if closure <= 0.001:
                return Operator(now[0], now[1], now[2], now[3], resolution, closure)
        prior = now
    raise ValueError(
        "doubled-resolution operator closure exceeds 0.1% at 128 subpixels"
    )
