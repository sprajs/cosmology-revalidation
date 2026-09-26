#!/usr/bin/env python3
"""Source-mapped BPIXTAB union and frozen circle support; no RAW SCI access."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
from astropy.io import fits

HERE = Path(__file__).resolve().parent
N = 1024
ACTIVE = slice(5, 1019)
R0, R1, R2 = (3.118908382066277, 9.35672514619883, 15.594541910331383)


def hash_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def map_bpixtab(path: Path) -> tuple[np.ndarray, dict]:
    """Exact IR doDQIIR/DQINormal map: one-index PIX, then raw LTV=+5."""
    mask = np.zeros((N, N), dtype=np.uint16)
    rows_used = 0
    clamped = 0
    with fits.open(path, memmap=True) as hdus:
        table = hdus["BPIX"]
        if table.header["SIZAXIS1"] != N or table.header["SIZAXIS2"] != N:
            raise ValueError("BPIXTAB detector dimension mismatch")
        for row in table.data:
            amp = str(row["CCDAMP"]).strip()
            chip = int(row["CCDCHIP"])
            gain = float(row["CCDGAIN"])
            if amp not in ("N/A", "ABCD") or chip not in (-999, 1) or gain not in (-999.0, 2.5):
                continue
            x = int(row["PIX1"]) - 1 + 5  # raw SCI LTV1=5, LTM1_1=1
            y = int(row["PIX2"]) - 1 + 5
            length = int(row["LENGTH"])
            axis = int(row["AXIS"])
            flag_float = float(row["VALUE"])
            flag = int(flag_float)
            if flag_float != flag or not 0 <= flag <= 65535 or length <= 0 or axis not in (1, 2):
                raise ValueError("invalid BPIXTAB row")
            rows_used += 1
            if axis == 1:
                if y < 0 or y >= N or x + length - 1 < 0 or x >= N:
                    clamped += 1
                    continue
                lo, hi = max(x, 0), min(x + length, N)
                clamped += int(lo != x or hi != x + length)
                mask[y, lo:hi] |= np.uint16(flag)
            else:
                if x < 0 or x >= N or y + length - 1 < 0 or y >= N:
                    clamped += 1
                    continue
                lo, hi = max(y, 0), min(y + length, N)
                clamped += int(lo != y or hi != y + length)
                mask[lo:hi, x] |= np.uint16(flag)
    return mask, {"rows_used": rows_used, "out_of_frame_or_clamped_rows": clamped,
                  "nonzero_pixels_all_raw": int(np.count_nonzero(mask)),
                  "nonzero_pixels_active": int(np.count_nonzero(mask[ACTIVE, ACTIVE]))}


def quadrant_area(x: float, y: float, radius: float) -> float:
    """Area of circle in [0,x]×[0,y] for x,y>=0."""
    x, y = min(max(x, 0.0), radius), min(max(y, 0.0), radius)
    if x == 0 or y == 0:
        return 0.0
    if x * x + y * y <= radius * radius:
        return x * y
    x0 = min(x, math.sqrt(max(radius * radius - y * y, 0.0)))

    def integral(t: float) -> float:
        return 0.5 * (t * math.sqrt(max(radius * radius - t * t, 0.0))
                      + radius * radius * math.asin(min(1.0, t / radius)))

    return y * x0 + integral(x) - integral(x0)


def signed_area(x: float, y: float, radius: float) -> float:
    return math.copysign(1.0, x) * math.copysign(1.0, y) * quadrant_area(abs(x), abs(y), radius)


def pixel_circle_fraction(x: int, y: int, radius: float) -> float:
    x0, x1, y0, y1 = x - 0.5, x + 0.5, y - 0.5, y + 0.5
    return (signed_area(x1, y1, radius) - signed_area(x0, y1, radius)
            - signed_area(x1, y0, radius) + signed_area(x0, y0, radius))


def footprint(radius: float) -> np.ndarray:
    extent = math.ceil(radius + 0.5)
    grid = np.array([[pixel_circle_fraction(x, y, radius)
                      for x in range(-extent, extent + 1)]
                     for y in range(-extent, extent + 1)], dtype=np.float64)
    if np.min(grid) < -1e-10 or np.max(grid) > 1 + 1e-10:
        raise AssertionError("circle fraction out of [0,1]")
    return grid


def synthetic_checks() -> dict:
    checks = {}
    for label, r in (("aperture", R0), ("inner", R1), ("outer", R2)):
        a = footprint(r)
        err = float(np.sum(a) - math.pi * r * r)
        if abs(err) > 1e-8:
            raise AssertionError(f"{label} area closure {err}")
        checks[label] = {"sum": float(np.sum(a)), "pi_r2": math.pi * r * r,
                         "area_error": err, "shape": list(a.shape)}
    # Independent geometric special cases and reflection symmetry.
    assert abs(pixel_circle_fraction(0, 0, 0.5) - math.pi * 0.25) < 1e-12
    for r in (R0, R1, R2):
        for x, y in ((1, 2), (3, 0), (10, 7)):
            assert abs(pixel_circle_fraction(x, y, r) - pixel_circle_fraction(-x, y, r)) < 1e-12
    # Four-row source-coordinate synthetic: 1-index PIX plus LTV5, X/Y runs, OR.
    grid = np.zeros((N, N), dtype=np.uint16)
    grid[5, 5:8] |= 4
    grid[5:8, 5] |= 8
    assert [int(grid[5, 5]), int(grid[5, 6]), int(grid[6, 5])] == [12, 4, 8]
    checks["row_geometry"] = "pass"
    return checks


def main() -> None:
    protocol = HERE / "mask-protocol.json"
    plan = json.loads(protocol.read_text())
    if plan["source_sha256"] != hash_file(Path(__file__)):
        raise SystemExit("source changed since pre-mask freeze")
    checks = synthetic_checks()
    files = [HERE / "references" / name for name in plan["reference_names"]]
    for path in files:
        if hash_file(path) != plan["reference_sha256"][path.name]:
            raise SystemExit(f"reference hash changed: {path}")
    masks = []
    table_stats = {}
    for path in files:
        mask, stats = map_bpixtab(path)
        masks.append(mask)
        table_stats[path.name] = stats
    union = np.bitwise_or.reduce(masks)
    bad = union != 0
    np.save(HERE / "fixed-union-bad.npy", bad)
    n_rejected = 0
    coverage = []
    outer_extent = math.ceil(R2 + 0.5)
    side = 2 * outer_extent + 1
    a = np.zeros((side, side), dtype=np.float64)
    b = np.zeros_like(a)
    for iy, dy in enumerate(range(-outer_extent, outer_extent + 1)):
        for ix, dx in enumerate(range(-outer_extent, outer_extent + 1)):
            a[iy, ix] = pixel_circle_fraction(dx, dy, R0)
            b[iy, ix] = pixel_circle_fraction(dx, dy, R2) - pixel_circle_fraction(dx, dy, R1)
    area_a, area_b = float(np.sum(a)), float(np.sum(b))
    for j in range(16):
        cy = 32 + 64 * j
        for i in range(16):
            cx = 32 + 64 * i
            if cx - outer_extent < 5 or cx + outer_extent > 1018 or cy - outer_extent < 5 or cy + outer_extent > 1018:
                raise AssertionError("fixed centre outside active support")
            good = ~bad[cy - outer_extent:cy + outer_extent + 1,
                        cx - outer_extent:cx + outer_extent + 1]
            frac_a = float(np.sum(good * a) / area_a)
            frac_b = float(np.sum(good * b) / area_b)
            eligible = frac_a >= 0.90 and frac_b >= 0.75
            n_rejected += not eligible
            weight = good * a - (np.sum(good * a) / np.sum(good * b)) * good * b if eligible else None
            if weight is not None and abs(float(np.sum(weight))) > 1e-8:
                raise AssertionError("background-balanced weight sum")
            coverage.append({"i": i, "j": j, "x": cx, "y": cy,
                             "aperture_coverage": frac_a, "annulus_coverage": frac_b,
                             "eligible": bool(eligible)})
    result = {"source_sha256": hash_file(Path(__file__)),
              "protocol_sha256": hash_file(protocol), "reference_sha256": plan["reference_sha256"],
              "synthetic": checks, "tables": table_stats,
              "union_bad_all_raw": int(np.count_nonzero(bad)),
              "union_bad_active": int(np.count_nonzero(bad[ACTIVE, ACTIVE])),
              "fixed_union_mask_sha256": hash_file(HERE / "fixed-union-bad.npy"),
              "support_total": 256, "support_eligible": 256 - n_rejected,
              "support_rejected": n_rejected, "coverage": coverage,
              "aperture_area": area_a, "annulus_area": area_b}
    (HERE / "fixed-mask-result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: result[k] for k in ("union_bad_active", "support_total", "support_eligible", "support_rejected")}, indent=2))


if __name__ == "__main__":
    main()
