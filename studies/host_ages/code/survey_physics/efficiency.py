#!/usr/bin/env python3
"""Independent bounded-efficiency and interpolation sensitivity checks."""
import gzip
import json
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from common import WORK, RESULTS, DATA, sha


def parse(text):
    maps = []
    current = None
    metadata = {}
    for line in text.splitlines():
        if line.startswith("PEAKMJD_RANGE:"):
            metadata["peak_range"] = [float(v) for v in line.split()[1:3]]
        elif line.startswith("FIELDLIST:"):
            metadata["fields"] = line.split()[1]
        elif line.startswith("VARNAMES:"):
            current = dict(metadata, names=line.split()[1:], rows=[])
            maps.append(current)
        elif line.startswith("HOSTEFF:"):
            current["rows"].append([float(v) for v in line.split()[1:4]])
    return maps


def main():
    path = DATA / "models/searcheff/SEARCHEFF_zHOST_DES-SN5YR_OBS.DAT.gz"
    maps = parse(gzip.open(path, "rt").read())
    rng = np.random.default_rng(2026092703)
    records = []
    for m in maps:
        a = np.array(m["rows"])
        x, y = np.unique(a[:, 0]), np.unique(a[:, 1])
        v = np.full((len(x), len(y)), np.nan)
        v[np.searchsorted(x, a[:, 0]), np.searchsorted(y, a[:, 1])] = a[:, 2]
        assert np.isfinite(v).all() and len(np.unique(a[:, :2], axis=0)) == len(a)
        bounded = np.clip(v, 0, 1)
        points = rng.uniform(
            [x.min() - 1, y.min() - 1], [x.max() + 1, y.max() + 1], (10000, 2)
        )
        points = np.clip(points, [x.min(), y.min()], [x.max(), y.max()])
        raw = RegularGridInterpolator((x, y), v)(points)
        before = RegularGridInterpolator((x, y), bounded)(points)
        after = np.clip(raw, 0, 1)
        # Independent explicit bilinear interpolation, including exact endpoints.
        ix = np.clip(np.searchsorted(x, points[:, 0], side="right") - 1, 0, len(x) - 2)
        iy = np.clip(np.searchsorted(y, points[:, 1], side="right") - 1, 0, len(y) - 2)
        tx = (points[:, 0] - x[ix]) / (x[ix + 1] - x[ix])
        ty = (points[:, 1] - y[iy]) / (y[iy + 1] - y[iy])
        manual = (
            (1 - tx) * (1 - ty) * bounded[ix, iy]
            + tx * (1 - ty) * bounded[ix + 1, iy]
            + (1 - tx) * ty * bounded[ix, iy + 1]
            + tx * ty * bounded[ix + 1, iy + 1]
        )
        delta = float(np.max(abs(manual - before)))
        assert delta < 1e-12 and before.min() >= -1e-12 and before.max() <= 1 + 1e-12
        # Held-out grid knots test interpolation, NOT new detection observations.
        predicted = 0.5 * (bounded[:-2] + bounded[2:])
        held = bounded[1:-1]
        u = rng.random(len(raw))
        assert np.array_equal(u < raw, u < after)
        records.append(
            {
                "fields": m["fields"],
                "peak_range": m["peak_range"],
                "grid_nodes": len(a),
                "invalid_original_nodes": int(((v < 0) | (v > 1)).sum()),
                "raw_probability_max": float(v.max()),
                "manual_vs_scipy_max": delta,
                "bounded_grid_vs_legacy_saturated_interpolation_max": float(
                    np.max(abs(before - after))
                ),
                "bounded_grid_vs_legacy_rms": float(
                    np.sqrt(np.mean((before - after) ** 2))
                ),
                "held_out_grid_knot_rmse": float(
                    np.sqrt(np.mean((held - predicted) ** 2))
                ),
                "legacy_threshold_equals_saturated_probability": True,
            }
        )
    result = {
        "code_sha256": sha(__file__),
        "input": {"path": str(path), "sha256": sha(path)},
        "current_native_original_map_failure": {
            "event": 37,
            "efficiency": 1.004,
            "status": "fatal guard; preserved .work/survey-physics/pilot.log",
        },
        "legacy_source": {
            "commit": "8ccf23a17cd98ac667dbc9e9666970c1e36246af",
            "date": "2025-11-12",
            "url": "https://github.com/RickKessler/SNANA/blob/8ccf23a17cd98ac667dbc9e9666970c1e36246af/src/sntools_trigger.c",
            "sha256": sha(WORK / "release-era-sntools_trigger.c"),
            "rule": "Uniform variate RAN>EFF, with no host-efficiency bound check. Mathematically identical to saturating the interpolated probability at0/1.",
        },
        "maps": records,
        "maps_tested": len(records),
        "random_interpolation_points": 10000 * len(records),
        "invalid_nodes": sum(r["invalid_original_nodes"] for r in records),
        "maximum_probability_change": max(
            r["bounded_grid_vs_legacy_saturated_interpolation_max"] for r in records
        ),
        "limitations": [
            "Held-out knots are released map values, not independent discovery/redshift follow-up events.",
            "Interpolation is validated numerically; the empirical selection model is not re-estimated here.",
            "Bounded-node interpolation differs slightly from legacy post-interpolation saturation; recorded response cannot assume equality.",
        ],
    }
    (RESULTS / "efficiency.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "maps"}, indent=2))


if __name__ == "__main__":
    main()
