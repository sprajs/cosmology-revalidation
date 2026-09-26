"""Independent ledger-only check of fixed pixel-noise summaries and resamples.

This script never opens image files or imports the scoring/recovery code.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def csv_rows(name):
    with (HERE / name).open(newline="") as f:
        return list(csv.DictReader(f))


def direct_metrics(records, idx):
    d = np.array([float(records[i]["d"]) for i in idx], dtype=np.float64)
    v = np.array([float(records[i]["V"]) for i in idx], dtype=np.float64)
    assert len(d) and np.all(np.isfinite(d)) and np.all(np.isfinite(v)) and np.all(v > 0)
    n = len(d)
    mu = np.dot(d, 1 / v) / (1 / v).sum()
    z = d / np.sqrt(v)
    return {"n": n, "weighted_intercept": mu, "signed_standardized_mean": z.mean(),
            "uncentered_mean_z2": np.dot(z, z) / n,
            "centered_variance_ratio": np.dot(d - mu, (d - mu) / v) / (n - 1) if n > 1 else None,
            "sum_d": d.sum(), "sum_d2": np.dot(d, d), "sum_V": v.sum(),
            "sum_d2_minus_sum_V": np.dot(d, d) - v.sum(),
            "sum_d2_over_sum_V": np.dot(d, d) / v.sum()}


def close(actual, expected, label, atol=1e-11):
    if not np.isclose(actual, expected, rtol=1e-12, atol=atol):
        raise AssertionError(f"{label}: {actual} != {expected}")


def main():
    frozen = json.loads((HERE / "recovery-protocol.json").read_text())
    for record in frozen["inputs"]:
        p = Path(record["path"])
        assert sha(p) == record["sha256"] and p.stat().st_size == record["bytes"]
    assert frozen["independent_script_sha256"] == sha(Path(__file__))
    recovered = json.loads((HERE / "stage_b_recovered_primary.json").read_text())
    design = json.loads((HERE.parent / "pixel_design" / "protocol.json").read_text())
    pairs = csv_rows("signed_pair_ledger.csv")
    exposure = csv_rows("signed_exposure_ledger.csv")
    template = csv_rows("common_template_ledger.csv")
    assert len(pairs) == 340 and len(exposure) == 1360 and len(template) == 170
    assert set(x["branch"] for x in pairs) == {"secondary"}
    checks = {}
    for visit in ("search", "template"):
        data = [r for r in pairs if r["visit"] == visit]
        target = recovered["branches"]["secondary"]["visits"][visit]
        assert len(data) == target["n"] == 170
        stats = direct_metrics(data, range(len(data)))
        for name, value in stats.items():
            if value is not None:
                close(value, target[name], visit + "." + name)
        tiles = sorted({int(r["tile"]) for r in data})
        lookup = {t: np.array([i for i, r in enumerate(data) if int(r["tile"]) == t], int)
                  for t in tiles}
        counts = {str(t): len(lookup[t]) for t in tiles}
        assert counts == target["tile_counts"]
        assert len(tiles) == 16 and min(counts.values()) >= 3
        delete = [direct_metrics(data, np.concatenate([lookup[q] for q in tiles if q != t]))
                  for t in tiles]
        for t, stats_delete in zip(tiles, delete):
            expected = next(r for r in target["tile_delete_one"] if r["removed_tile"] == t)
            for name, value in stats_delete.items():
                if value is not None:
                    close(value, expected[name], f"{visit}.delete{t}.{name}")
        rng = np.random.default_rng(design["bootstrap"]["seed"])
        boot = []
        for i in range(design["bootstrap"]["spatial_tile_resamples"]):
            draw = rng.choice(tiles, size=len(tiles), replace=True)
            subset = np.concatenate([lookup[int(t)] for t in draw])
            m = direct_metrics(data, subset)
            expected = target["tile_bootstrap_replicates"][i]
            for name, value in m.items():
                if value is not None:
                    close(value, expected[name], f"{visit}.bootstrap{i}.{name}")
            boot.append(m["centered_variance_ratio"])
        quantiles = np.quantile(boot, [.025, .975])
        for i, value in enumerate(quantiles):
            close(value, target["tile_bootstrap_centered_ratio_2p5_97p5"][i],
                  f"{visit}.bootstrap_quantile{i}")
        checks[visit] = {"n": len(data), "unique_tiles": len(tiles),
                         "min_per_tile": min(counts.values()), "bootstrap_replicates": len(boot),
                         "raw_moments": {k: float(stats[k]) for k in
                             ("sum_d", "sum_d2", "sum_V", "sum_d2_minus_sum_V",
                              "sum_d2_over_sum_V")},
                         "signed_standardized_mean": float(stats["signed_standardized_mean"]),
                         "uncentered_mean_z2": float(stats["uncentered_mean_z2"]),
                         "centered_variance_ratio": float(stats["centered_variance_ratio"])}
    keyed = {(r["branch"], r["candidate"], r["visit"], r["ordinal"]): r for r in exposure}
    for r in template:
        cid = r["candidate"]
        s2,s4,t2,t4 = [keyed[("secondary",cid,visit,ordinal)] for visit,ordinal in
                       (("search","2"),("search","4"),("template","2"),("template","4"))]
        x = np.array([float(q["flux"]) for q in (s2,s4,t2,t4)])
        v = np.array([float(q["V"]) for q in (s2,s4,t2,t4)])
        common = (x[0] - .5 * (x[2] + x[3]), x[1] - .5 * (x[2] + x[3]))
        separate = (x[0]-x[2], x[1]-x[3])
        close(common[0], float(r["common1"]), cid + ".common1")
        close(common[1], float(r["common2"]), cid + ".common2")
        close(separate[0], float(r["separate1"]), cid + ".separate1")
        close(separate[1], float(r["separate2"]), cid + ".separate2")
        close(.25 * (v[2] + v[3]), float(r["predicted_common_cov_offdiag"]), cid + ".cov")
        close(common[0]*common[1]-separate[0]*separate[1],
              float(r["crossproduct_difference"]), cid + ".crossproduct")
    result = {"status": "PASS ledger-only independent primary arithmetic and all tile resamples",
              "no_fits_reopen": True, "strict_gate": recovered["branches"]["strict"]["primary_adequacy"],
              "recovered_sha256": sha(HERE / "stage_b_recovered_primary.json"),
              "recovery_protocol_sha256": sha(HERE / "recovery-protocol.json"),
              "visits": checks, "common_template_objects_checked": len(template)}
    dest = HERE / "independent_statistics.json"
    dest.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
