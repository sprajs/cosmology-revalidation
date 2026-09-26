#!/usr/bin/env python3
"""Frozen descriptive accounting for coherent 2021 RAISIN simulated FITRES."""
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "runs/research_2026_09_26/raisin_2021_selection_accounting"
IN = ROOT / "runs/research_2026_09_26/astra_design/raisin_timing_assets/author_20211111"
BINS = np.array([0, .2, .3, .4, .5, .6, 1.])
Q = [0, .05, .16, .5, .84, .95, 1]
COLS = ["CID", "zHD", "PKMJDINI", "PKMJD", "SIM_PKMJD", "SIM_AV", "SIM_STRETCH", "SIM_RV", "AV", "STRETCH", "RV", "DLMAG", "SIM_DLMAG"]


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read_fitres(path):
    cols = None
    out = {c: [] for c in COLS}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("VARNAMES:"):
                cols = line.split()[1:]
                idx = {c: cols.index(c) + 1 for c in COLS}
            elif line.startswith("SN:"):
                a = line.split()
                if len(a) != len(cols) + 1:
                    raise ValueError((path, "field count", len(a), len(cols)))
                for c in COLS:
                    out[c].append(a[idx[c]])
    if cols is None:
        raise ValueError((path, "missing VARNAMES"))
    return {c: np.array(out[c], dtype=int if c == "CID" else float) for c in COLS}


def stats(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if not len(x):
        return {"n": 0, "mean": None, "sd": None, "q": None}
    return {"n": int(len(x)), "mean": float(x.mean()), "sd": float(x.std()),
            "q": [float(v) for v in np.quantile(x, Q)]}


def main():
    nirpath = IN / "nir.FITRES.gz"
    optpath = IN / "optnir.FITRES.gz"
    n = read_fitres(nirpath)
    o = read_fitres(optpath)
    if len(set(n["CID"])) != len(n["CID"]) or len(set(o["CID"])) != len(o["CID"]):
        raise ValueError("duplicate CID")
    if not np.all(np.isfinite(n["SIM_DLMAG"])) or not np.all(np.isfinite(o["SIM_DLMAG"])) or np.min(n["SIM_DLMAG"]) <= 0 or np.min(o["SIM_DLMAG"]) <= 0:
        raise ValueError("SIM_DLMAG missing/sentinel; refuse residual summaries")
    oi = {int(v): i for i, v in enumerate(o["CID"])}
    oi_for_n = np.array([oi.get(int(v), -1) for v in n["CID"]])
    joint = oi_for_n >= 0
    avpass = np.zeros(len(n["CID"]), bool)
    stpass = np.zeros(len(n["CID"]), bool)
    avpass[joint] = o["AV"][oi_for_n[joint]] < .3 * o["RV"][oi_for_n[joint]]
    stpass[joint] = (o["STRETCH"][oi_for_n[joint]] > .75) & (o["STRETCH"][oi_for_n[joint]] < 1.185)
    selected = joint & avpass & stpass
    if len(n["CID"]) != 30000 or len(o["CID"]) != 29995 or int(joint.sum()) != 29995:
        raise ValueError("expected archived cohort count failed")
    keys = {
        "sim_av": n["SIM_AV"], "sim_stretch": n["SIM_STRETCH"], "sim_rv": n["SIM_RV"],
        "nir_fit_av": n["AV"], "nir_fit_stretch": n["STRETCH"],
        "nir_peak_init_minus_truth": n["PKMJDINI"] - n["SIM_PKMJD"],
        "nir_peak_fit_minus_truth": n["PKMJD"] - n["SIM_PKMJD"],
        "nir_d_minus_sim_d": n["DLMAG"] - n["SIM_DLMAG"],
    }
    jo = oi_for_n[joint]
    jointkeys = {
        "joint_fit_av": o["AV"][jo], "joint_fit_stretch": o["STRETCH"][jo], "joint_fit_rv": o["RV"][jo],
        "joint_peak_fit_minus_truth": o["PKMJD"][jo] - n["SIM_PKMJD"][joint],
        "joint_minus_nir_peak": o["PKMJD"][jo] - n["PKMJD"][joint],
        "joint_d_minus_sim_d": o["DLMAG"][jo] - o["SIM_DLMAG"][jo],
    }
    summary = {
        "counts": {"nir": len(n["CID"]), "joint": int(joint.sum()), "missing_joint": int((~joint).sum()),
                   "fail_av_only": int((joint & ~avpass & stpass).sum()),
                   "fail_stretch_only": int((joint & avpass & ~stpass).sum()),
                   "fail_both": int((joint & ~avpass & ~stpass).sum()),
                   "pass_both": int(selected.sum())},
        "missing_joint_cids": n["CID"][~joint].astype(int).tolist(),
        "all_nir": {k: stats(v) for k, v in keys.items()},
        "joint_present": {k: stats(v) for k, v in jointkeys.items()},
        "selected_nir": {k: stats(v[selected]) for k, v in keys.items()},
        "selected_joint": {k: stats(v[selected[joint]]) for k, v in jointkeys.items()},
        "z_bins": [],
    }
    for lo, hi in zip(BINS[:-1], BINS[1:]):
        zb = (n["zHD"] >= lo) & ((n["zHD"] < hi) if hi < 1 else (n["zHD"] <= hi))
        zj = zb & joint
        zs = zb & selected
        summary["z_bins"].append({"lo": float(lo), "hi": float(hi), "nir": int(zb.sum()),
                                  "joint": int(zj.sum()), "selected": int(zs.sum()),
                                  "sim_av_before": stats(n["SIM_AV"][zb]),
                                  "sim_av_after": stats(n["SIM_AV"][zs]),
                                  "sim_stretch_before": stats(n["SIM_STRETCH"][zb]),
                                  "sim_stretch_after": stats(n["SIM_STRETCH"][zs]),
                                  "sim_rv_before": stats(n["SIM_RV"][zb]),
                                  "sim_rv_after": stats(n["SIM_RV"][zs]),
                                  "timing_nir_before": stats(keys["nir_peak_fit_minus_truth"][zb]),
                                  "timing_nir_after": stats(keys["nir_peak_fit_minus_truth"][zs]),
                                  "nir_d_minus_sim_d_before": stats(keys["nir_d_minus_sim_d"][zb]),
                                  "nir_d_minus_sim_d_after": stats(keys["nir_d_minus_sim_d"][zs])})
    RUN.mkdir(parents=True, exist_ok=True)
    (RUN / "result.json").write_text(json.dumps(summary, indent=2) + "\n")
    paths = [nirpath, optpath, RUN / "protocol.json", RUN / "source/cosmo_sys-20211107.py", Path(__file__)]
    (RUN / "manifest.json").write_text(json.dumps({str(p.relative_to(ROOT)): digest(p) for p in paths}, indent=2) + "\n")


if __name__ == "__main__":
    main()
