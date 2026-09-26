#!/usr/bin/env python3
"""Independent physical-uint16 RAW extraction of the fixed dark estimands.

This verifies arithmetic, not detector-noise attribution. It imports only the
previously independently checked circle/square geometry, never the scorer.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "runs/research_2026_09_26/raisin_hst_pixel_feasibility"
D = BASE / "dark_ramp_execution"
OUT = D / "root-raw-review"


def sha(path):
    return hashlib.file_digest(Path(path).open("rb"), "sha256").hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def rows(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream))


def main():
    OUT.mkdir(exist_ok=False)
    started = time.monotonic()
    design_path = BASE / "pixel_design/variance-diagnostic-review/calwf3_read_matrix/dark-quality-design/experiment-protocol.json"
    design = read_json(design_path)
    freeze = read_json(D / "execution-freeze.json")
    header = read_json(D / "header-gate.json")
    paths = {x["root"]: Path(x["path"]) for x in header["records"]}
    for root, path in paths.items():
        assert sha(path) == freeze["raw_sha256"][root]
    assert sha(design_path) == freeze["design_sha256"]
    assert sha(D / "fixed-union-bad.npy") == freeze["mask_sha256"]
    spec = importlib.util.spec_from_file_location("checked_geometry", D / "fixed_mask.py")
    geometry = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(geometry)
    assert sha(D / "fixed_mask.py") == freeze["source_files"][str(D / "fixed_mask.py")]
    inputs = [design_path, D / "execution-freeze.json", D / "header-gate.json",
              D / "fixed-union-bad.npy", D / "fixed-mask-result.json", D / "fixed_mask.py",
              D / "read-matrix.npz", D / "pair-quadrant-power.csv",
              D / "block-pair-power.csv", D / "aperture-pair-power.csv"]
    protocol = {
        "source_sha256": sha(__file__), "input_sha256": {str(p): sha(p) for p in inputs},
        "method": "Astropy physical uint16; float64 before difference; independent scalar means, dot products and slope sums; all saved fixed regions and powers",
        "shared_geometry": "Only circle/square helper, already independently quadrature-validated in root-mask-review; weights assembled independently here",
        "tolerance": "rtol=1e-10 and atol=1e-10 for slopes/moments; matrices atol=1e-8; no rescaling, new cuts or science attribution",
        "pairs": 8, "apertures_per_pair": 247, "time_cap_seconds": 180,
    }
    (OUT / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    bad = np.load(D / "fixed-union-bad.npy")
    t = np.array(design["time_operator"]["t_seconds_nominal"])
    powers = design["time_operator"]["powers"]
    h = []
    for power in powers:
        u = np.abs(np.arange(-3, 4) / 3.0) ** power
        center = np.dot(t, u) / u.sum()
        row = u * (t - center) / np.dot(u, (t - center) ** 2)
        assert abs(row.sum()) < 1e-16 and abs(np.dot(row, t) - 1) < 1e-12
        h.append(row)
    h = np.array(h)
    saved = np.load(D / "read-matrix.npz")
    quadrants = [("B", slice(5, 512), slice(5, 512)),
                 ("C", slice(5, 512), slice(512, 1019)),
                 ("A", slice(512, 1019), slice(5, 512)),
                 ("D", slice(512, 1019), slice(512, 1019))]
    qr = {(r["pair"], r["scope"], r["quadrant"], float(r["power"])): r
          for r in rows(D / "pair-quadrant-power.csv")}
    br = {(r["pair"], int(r["block_y"]), int(r["block_x"])): r
          for r in rows(D / "block-pair-power.csv")}
    ar = {(r["pair"], int(r["y"]), int(r["x"]), float(r["power"])): r
          for r in rows(D / "aperture-pair-power.csv")}
    r0, r1, r2 = design["spatial_operator"]["radii_pixels"]
    radius = int(np.ceil(r2 + 0.5))
    offsets = range(-radius, radius + 1)
    a = np.array([[geometry.pixel_circle_fraction(x, y, r0) for x in offsets] for y in offsets])
    b = np.array([[geometry.pixel_circle_fraction(x, y, r2) - geometry.pixel_circle_fraction(x, y, r1)
                   for x in offsets] for y in offsets])
    sites = []
    for entry in read_json(D / "fixed-mask-result.json")["coverage"]:
        y, x = entry["y"], entry["x"]
        yy, xx = slice(y-radius, y+radius+1), slice(x-radius, x+radius+1)
        ma, mb = a * ~bad[yy, xx], b * ~bad[yy, xx]
        eligible = ma.sum()/a.sum() >= .9 and mb.sum()/b.sum() >= .75
        assert eligible == entry["eligible"]
        if eligible:
            w = ma - mb * (ma.sum()/mb.sum())
            assert abs(w.sum()) < 1e-8
            sites.append((y, x, yy, xx, w))
    assert len(sites) == 247
    gaps, counts = {}, {}

    def check(label, actual, expected, atol=1e-10):
        actual, expected = np.asarray(actual), np.asarray(expected)
        assert np.allclose(actual, expected, rtol=1e-10, atol=atol), label
        gaps[label] = max(gaps.get(label, 0.0), float(np.max(np.abs(actual-expected))))
        counts[label] = counts.get(label, 0) + actual.size

    for ip, pair in enumerate(design["pairs"]):
        assert time.monotonic()-started < 180
        delta = np.empty((7, 1024, 1024), dtype=np.float64)
        with fits.open(paths[pair["earlier"]], memmap=False, uint=True) as earlier, \
             fits.open(paths[pair["later"]], memmap=False, uint=True) as later:
            for it, extver in enumerate(range(15, 8, -1)):
                raw_a, raw_b = earlier["SCI", extver].data, later["SCI", extver].data
                assert raw_a.dtype == np.uint16 and raw_b.dtype == np.uint16
                delta[it] = raw_b.astype(np.float64) - raw_a.astype(np.float64)

        def region(ys, xs, good):
            v = delta[:, ys, xs][:, good]
            mean = np.array([np.sum(row)/v.shape[1] for row in v])
            gamma = np.array([[np.sum(left*right)/(2*v.shape[1]) for right in v] for left in v])
            return v, mean, gamma

        for iq, (name, yy, xx) in enumerate(quadrants):
            for scope_id, (scope, mask) in enumerate((
                    ("masked", ~bad[yy, xx]), ("geometry_only", np.ones(bad[yy, xx].shape, bool)))):
                v, mean, gamma = region(yy, xx, mask)
                check("quadrant_mean", mean, saved["quadrant_mean"][ip, scope_id, iq])
                check("quadrant_gamma", gamma, saved["quadrant_gamma"][ip, scope_id, iq], 1e-8)
                for power, hh in zip(powers, h):
                    slopes = np.sum(hh[:, None]*v, axis=0)
                    score = np.dot(slopes, slopes)/(2*len(slopes))
                    expected = qr[pair["pair"], scope, name, float(power)]
                    assert int(expected["n"]) == v.shape[1]
                    check("quadrant_score", score, float(expected["dn2_per_nominal_s2"]))
        for by in range(16):
            yy = slice(max(5, by*64), min(1019, (by+1)*64))
            for bx in range(16):
                xx = slice(max(5, bx*64), min(1019, (bx+1)*64))
                v, mean, gamma = region(yy, xx, ~bad[yy, xx])
                check("block_mean", mean, saved["block_mean"][ip, by, bx])
                check("block_gamma", gamma, saved["block_gamma"][ip, by, bx], 1e-8)
                slopes = np.sum(h[0, :, None]*v, axis=0)
                check("block_score", np.mean(slopes**2)/2,
                      float(br[pair["pair"], by, bx]["primary_dn2_per_nominal_s2"]))
        for y, x, yy, xx, w in sites:
            for power, hh in zip(powers, h):
                field = np.sum(delta[:, yy, xx] * hh[:, None, None], axis=0)
                scalar = float(np.sum(w * field))
                expected = ar[pair["pair"], y, x, float(power)]
                check("aperture_slope", scalar, float(expected["signed_dn_per_nominal_s"]))
                check("aperture_half_square", scalar**2/2, float(expected["half_squared_dn2_per_nominal_s2"]))
        print(json.dumps({"pair": pair["pair"], "pass": True}), flush=True)
    result = {"pass": True, "elapsed_seconds": time.monotonic()-started,
              "protocol_sha256": sha(OUT / "protocol.json"), "max_abs_gap": gaps,
              "scalar_comparison_counts": counts,
              "scope": "All 8 fixed pair raw re-extractions; arithmetic only, no new noise attribution"}
    (OUT / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
