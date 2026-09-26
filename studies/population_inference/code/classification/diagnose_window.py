#!/usr/bin/env python3
"""Separate source-motivated DES clump-window diagnostic; preserves primary failure."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reconstruct import MODEL, OUT, RAW, REF, file_sha, infer, load_curves

CLUMP = OUT / "clump_run/DES-SN5YR_DES.SNANA.TEXT"
PLAN = OUT / "clump_correction_plan.md"


def main():
    plan_hash = file_sha(PLAN)
    clump = pd.read_csv(CLUMP, comment="#", delimiter=" ", skipinitialspace=True,
                        dtype={"CID": str})
    assert clump.CID.is_unique
    peaks = dict(zip(clump.CID, clump.PKMJDINI.astype(float)))
    ref = pd.read_csv(REF, dtype={"CID": str})
    ids = set(ref.CID)
    assert len(ids) == 1635 and ids <= peaks.keys()
    assert np.isfinite([peaks[cid] for cid in ids]).all()
    df, meta = load_curves(RAW / "DES-SN5YR_DES", ids, peaks)
    print(f"clump diagnostic: {len(meta)} HEAD; {len(df)} retained epochs; "
          f"{meta.n_window_excluded.sum()} window exclusions; "
          f"{meta.n_flag_excluded.sum()} flag exclusions", flush=True)
    pred = infer(df)
    out = meta.merge(pred, on="CID", validate="one_to_one").merge(
        ref[["CID", "PROB_SNNV19"]], on="CID", validate="one_to_one")
    prior = pd.read_csv(OUT / "des_probabilities.csv", dtype={"CID": str})
    out = out.merge(prior[["CID", "pIa"]].rename(columns={"pIa": "primary_pIa"}),
                    on="CID", validate="one_to_one")
    out["diff"] = out.pIa - out.PROB_SNNV19
    out["abs_diff"] = out["diff"].abs()
    out["window_effect"] = out.pIa - out.primary_pIa
    out.sort_values("CID").to_csv(OUT / "des_clump_diagnostic.csv", index=False)
    err = out["diff"].to_numpy()
    ab = np.abs(err)
    metrics = {"status": "diagnostic current-SNANA clump recreation, not exact historical clump",
               "plan_sha256": plan_hash, "clump_sha256": file_sha(CLUMP),
               "model_sha256": file_sha(MODEL), "n": len(out),
               "n_raw": int(meta.n_raw.sum()),
               "n_window_excluded": int(meta.n_window_excluded.sum()),
               "n_flag_excluded": int(meta.n_flag_excluded.sum()),
               "n_retained": int(meta.n_retained.sum()),
               "median_abs": float(np.median(ab)), "p95_abs": float(np.quantile(ab, .95)),
               "p99_abs": float(np.quantile(ab, .99)), "max_abs": float(ab.max()),
               "rms": float(np.sqrt(np.mean(err * err))), "bias": float(np.mean(err)),
               "spearman": float(spearmanr(out.pIa, out.PROB_SNNV19).statistic),
               "n_window_effect_gt_0p1": int((out.window_effect.abs() > .1).sum())}
    for label, cut in (("gt999", .999), ("gt5", .5)):
        p = out.pIa.to_numpy() > cut
        r = out.PROB_SNNV19.to_numpy() > cut
        metrics[label] = {"agreement": float(np.mean(p == r)), "both": int(np.sum(p & r)),
                          "pred_only": int(np.sum(p & ~r)), "released_only": int(np.sum(~p & r)),
                          "neither": int(np.sum(~p & ~r))}
    metrics["adequate_same_registered_gate"] = bool(
        metrics["rms"] <= .005 and metrics["p99_abs"] <= .02 and
        metrics["spearman"] >= .995 and metrics["gt5"]["agreement"] >= .99 and
        metrics["gt999"]["agreement"] >= .98)
    (OUT / "des_clump_diagnostic_metrics.json").write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    print(json.dumps(metrics, indent=2), flush=True)


if __name__ == "__main__":
    main()
