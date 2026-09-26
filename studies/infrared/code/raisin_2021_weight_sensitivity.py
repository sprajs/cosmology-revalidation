#!/usr/bin/env python3
"""Frozen source-weight sensitivity on the coherent 2021 RAISIN fit cohort."""
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "runs/research_2026_09_26/raisin_2021_selection_accounting"
OUT = BASE / "weight_sensitivity"
DATA = ROOT / "runs/research_2026_09_26/astra_design/raisin_timing_assets/author_20211111"
BINS = (0., .2, .3, .4, .5, .6, 1.)
SIGMAS = tuple(.005 * k for k in range(60))
NEEDED_NIR = ("CID", "zHD", "DLMAG", "DLMAGERR", "SIM_DLMAG")
NEEDED_OPT = ("CID", "AV", "RV", "STRETCH")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path, needed):
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("VARNAMES:"):
                names = line.split()[1:]
                positions = tuple(names.index(c) + 1 for c in needed)
                out = {c: [] for c in needed}
            elif line.startswith("SN:"):
                values = line.split()
                if len(values) != len(names) + 1:
                    raise ValueError((path, "column mismatch"))
                for col, pos in zip(needed, positions):
                    out[col].append(values[pos])
    return {c: np.asarray(v, dtype=int if c == "CID" else float) for c, v in out.items()}


def main():
    parent = json.loads((BASE / "manifest.json").read_text())
    for rel, expected in parent.items():
        if sha(ROOT / rel) != expected:
            raise ValueError(("parent input changed", rel))
    nir = read(DATA / "nir.FITRES.gz", NEEDED_NIR)
    opt = read(DATA / "optnir.FITRES.gz", NEEDED_OPT)
    if len(nir["CID"]) != 30000 or len(opt["CID"]) != 29995:
        raise ValueError("archived FITRES size changed")
    if len(set(nir["CID"])) != len(nir["CID"]) or len(set(opt["CID"])) != len(opt["CID"]):
        raise ValueError("duplicate CID")
    omap = {int(c): i for i, c in enumerate(opt["CID"])}
    selected = []
    for j, cid in enumerate(nir["CID"]):
        i = omap.get(int(cid))
        if i is not None and opt["AV"][i] < .3 * opt["RV"][i] and .75 < opt["STRETCH"][i] < 1.185:
            selected.append(j)
    if len(selected) != 25606:
        raise ValueError(("parent cohort closure", len(selected)))
    selected = np.asarray(selected, int)
    z = nir["zHD"][selected]
    resid = nir["DLMAG"][selected] - nir["SIM_DLMAG"][selected]
    err = nir["DLMAGERR"][selected]
    if not np.all(np.isfinite(z)) or not np.all(np.isfinite(resid)) or not np.all(np.isfinite(err)) or np.any(err <= 0):
        raise ValueError("nonfinite or nonpositive fit inputs")
    regions = {"global": np.ones(len(selected), bool)}
    for lo, hi in zip(BINS[:-1], BINS[1:]):
        name = f"[{lo:g},{hi:g}{']' if hi == 1. else ')'}"
        regions[name] = (z >= lo) & ((z <= hi) if hi == 1. else (z < hi))
    if sum(int(v.sum()) for k, v in regions.items() if k != "global") != len(selected):
        raise ValueError("fixed z-bin partition failure")
    rows = []
    summary = {"selected_n": len(selected), "selected_cids_sha256": hashlib.sha256("\n".join(map(str, nir["CID"][selected])).encode()).hexdigest(), "regions": {}}
    for name, mask in regions.items():
        n = int(mask.sum())
        region = {"n": n, "mean_ranges": {}, "p2_minus_p1_range": None}
        paired_diffs = []
        for s in SIGMAS:
            means = {}
            for p in (1, 2):
                if n:
                    weights = (err[mask] ** 2 + s ** 2) ** (-p / 2)
                    mean = float(np.dot(weights, resid[mask]) / weights.sum())
                    neff = float(weights.sum() ** 2 / np.dot(weights, weights))
                else:
                    mean = neff = None
                rows.append({"region": name, "n": n, "sigma": format(s, ".3f"), "p": p, "weighted_mean": mean, "effective_n": neff})
                means[p] = mean
            if n:
                paired_diffs.append(means[2] - means[1])
        for p in (1, 2):
            values = [r["weighted_mean"] for r in rows if r["region"] == name and r["p"] == p and r["weighted_mean"] is not None]
            region["mean_ranges"][str(p)] = {"min": min(values), "max": max(values)} if values else None
        region["p2_minus_p1_range"] = {"min": min(paired_diffs), "max": max(paired_diffs)} if paired_diffs else None
        summary["regions"][name] = region
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "curves.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=("region", "n", "sigma", "p", "weighted_mean", "effective_n"))
        writer.writeheader()
        writer.writerows(rows)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    paths = [OUT / "protocol.json", Path(__file__), BASE / "manifest.json", DATA / "nir.FITRES.gz", DATA / "optnir.FITRES.gz"]
    (OUT / "manifest.json").write_text(json.dumps({str(p.relative_to(ROOT)): sha(p) for p in paths}, indent=2) + "\n")


if __name__ == "__main__":
    main()
