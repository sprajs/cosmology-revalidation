#!/usr/bin/env python3
"""Reconstruct DES SNNV19 probabilities with pinned official SuperNNova inference."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / "phase2/classification/sources/SuperNNova"
sys.path.insert(0, str(SOURCE))
from supernnova.validation.validate_onthefly import classify_lcs  # noqa: E402

OUT = ROOT / "phase2/classification/reconstruction_20260926"
RAW = ROOT / "sources/repos/des-science__DES-SN5YR@1.3/0_DATA/DES-SN5YR_DES"
REF = ROOT / "sources/repos/des-science__DES-SN5YR@1.3/3_CLASSIFICATION/DES_classification.csv"
MODEL = ROOT / "phase2/official/inputs/SNDATA_ROOT/models/classifiers/DES-SN5YR/SNNTRAINV19_z_TRAINDES_V19/model.pt"
MASK = sum((8, 16, 32, 64, 128, 256, 512))
VARIANTS = ("P21", "BS21", "G10", "P21_dmplus010", "P21_dmminus010", "P21_dmz020", "P21_rho000", "P21_rho090", "P21_noisetrue120")


def file_sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_curves(stem: Path, ids: set[str] | None = None,
                peak_mjd: dict[str, float] | None = None):
    with fits.open(str(stem) + "_HEAD.FITS" + (".gz" if (Path(str(stem) + "_HEAD.FITS.gz")).exists() else ""), memmap=False) as f:
        head = f[1].data.copy()
    with fits.open(str(stem) + "_PHOT.FITS" + (".gz" if (Path(str(stem) + "_PHOT.FITS.gz")).exists() else ""), memmap=False) as f:
        phot = f[1].data.copy()
    names = [str(x).strip() for x in head["SNID"]]
    assert len(names) == len(set(names)), "duplicate HEAD SNID"
    rows = []
    meta = []
    for h, name in zip(head, names):
        if ids is not None and name not in ids:
            continue
        lo, hi = int(h["PTROBS_MIN"]), int(h["PTROBS_MAX"])
        nobs = int(h["NOBS"])
        assert 1 <= lo <= hi <= len(phot) and hi - lo + 1 == nobs, (name, lo, hi, nobs)
        q = phot[lo - 1:hi]
        assert len(q) == nobs and np.all(np.isfinite(q["MJD"]))
        assert np.all(np.isfinite(q["FLUXCAL"])) and np.all(np.isfinite(q["FLUXCALERR"]))
        # SNANA's delimiter row, if present, is outside inclusive PTROBS bounds.
        assert np.all(q["MJD"] > 0), name
        if peak_mjd is not None:
            assert name in peak_mjd, name
            dt = np.asarray(q["MJD"], dtype=float) - peak_mjd[name]
            window_keep = (dt > -30) & (dt < 100)
            n_window_excluded = int(np.sum(~window_keep))
            q = q[window_keep]
        flag = np.asarray(q["PHOTFLAG"], dtype=np.int64)
        keep = (flag & MASK) == 0
        n_flag_excluded = int(np.sum(~keep))
        q = q[keep]
        assert len(q) > 0, name
        bands = [str(v).strip() for v in q["BAND"]]
        assert set(bands).issubset({"g", "i", "r", "z"}), (name, set(bands))
        specz = float(h["REDSHIFT_FINAL"])
        specz_err = float(h["REDSHIFT_FINAL_ERR"])
        assert np.isfinite(specz) and np.isfinite(specz_err), name
        snr = np.asarray(q["FLUXCAL"] / q["FLUXCALERR"], dtype=float)
        assert np.all(np.isfinite(snr)), name
        rows.extend({"SNID": name, "MJD": float(p["MJD"]), "FLT": band,
                     "FLUXCAL": float(p["FLUXCAL"]), "FLUXCALERR": float(p["FLUXCALERR"]),
                     "HOSTGAL_SPECZ": specz, "HOSTGAL_SPECZ_ERR": specz_err}
                    for p, band in zip(q, bands))
        item = {"CID": name, "redshift_final": specz, "n_raw": nobs,
                "n_rejected": nobs - len(q), "n_retained": len(q), "peak_snr": float(np.max(snr))}
        if peak_mjd is not None:
            item.update({"clump_pkmjdini": peak_mjd[name], "n_window_excluded": n_window_excluded,
                         "n_flag_excluded": n_flag_excluded})
        meta.append(item)
    if ids is not None:
        assert {m["CID"] for m in meta} == ids, "classification IDs absent from HEAD"
    return pd.DataFrame(rows), pd.DataFrame(meta)


def infer(df: pd.DataFrame, chunk_size: int = 128,
          model_path: Path = MODEL) -> pd.DataFrame:
    ids = sorted(df.SNID.unique())
    res = []
    for start in range(0, len(ids), chunk_size):
        chunk = set(ids[start:start + chunk_size])
        idx, pred = classify_lcs(df.loc[df.SNID.isin(chunk)].copy(), str(model_path), "cpu")
        assert pred.shape == (len(chunk), 1, 2), pred.shape
        res.extend((str(cid), float(p[0, 0])) for cid, p in zip(idx, pred))
        print(f"inferred {min(start + chunk_size, len(ids))}/{len(ids)}", flush=True)
    out = pd.DataFrame(res, columns=["CID", "pIa"])
    assert len(out) == len(ids) and out.CID.is_unique
    assert np.isfinite(out.pIa).all() and out.pIa.between(0, 1).all()
    return out


def run_des():
    prereg = OUT / "preregistration.md"
    expected = "925c07d183ae9a4b73cad5337e2210320301133aafa9be8f3145cde1f44591a8"
    assert file_sha(prereg) == expected, "preregistration changed"
    ref = pd.read_csv(REF, dtype={"CID": str})
    ids = set(ref.CID)
    assert len(ref) == 1635 and len(ids) == 1635
    df, meta = load_curves(RAW / "DES-SN5YR_DES", ids)
    print(f"loaded {len(meta)} heads; {len(df)} retained rows; {meta.n_rejected.sum()} rejected", flush=True)
    subids = [sorted(ids)[i] for i in (0, 1, 2, 10, 20, 50, 100, 200, 500, 1000)]
    subset = df.loc[df.SNID.isin(subids)]
    batch = infer(subset, 128).set_index("CID").pIa
    singles = infer(subset, 1).set_index("CID").pIa
    max_batch_delta = float(np.max(np.abs(batch.sort_index().values - singles.sort_index().values)))
    assert max_batch_delta <= 1e-6, max_batch_delta
    pred = infer(df)
    full = meta.merge(pred, on="CID", validate="one_to_one").merge(
        ref[["CID", "PROB_SNNV19"]], on="CID", validate="one_to_one")
    full["diff"] = full.pIa - full.PROB_SNNV19
    full["abs_diff"] = full["diff"].abs()
    full["pred_gt999"] = full.pIa > 0.999
    full["released_gt999"] = full.PROB_SNNV19 > 0.999
    full["pred_gt5"] = full.pIa > 0.5
    full["released_gt5"] = full.PROB_SNNV19 > 0.5
    full.sort_values("CID").to_csv(OUT / "des_probabilities.csv", index=False)
    err = full["diff"].to_numpy()
    ab = np.abs(err)
    metrics = {"n": len(full), "preregistration_sha256": expected,
               "model_sha256": file_sha(MODEL), "max_batch_delta": max_batch_delta,
               "n_raw": int(meta.n_raw.sum()), "n_rejected": int(meta.n_rejected.sum()),
               "n_retained": int(meta.n_retained.sum()), "median_abs": float(np.median(ab)),
               "p95_abs": float(np.quantile(ab, .95)), "p99_abs": float(np.quantile(ab, .99)),
               "max_abs": float(np.max(ab)), "rms": float(np.sqrt(np.mean(err * err))),
               "bias": float(np.mean(err)),
               "spearman": float(spearmanr(full.pIa, full.PROB_SNNV19).statistic),
               "rounded_release_matches": int(np.sum(ab <= .00005))}
    for label in ("gt999", "gt5"):
        a = full[f"pred_{label}"].to_numpy()
        b = full[f"released_{label}"].to_numpy()
        metrics[label] = {"agreement": float(np.mean(a == b)), "both": int(np.sum(a & b)),
                          "pred_only": int(np.sum(a & ~b)), "released_only": int(np.sum(~a & b)),
                          "neither": int(np.sum(~a & ~b))}
    for col in ("redshift_final", "peak_snr", "n_retained"):
        full[f"{col}_quartile"] = pd.qcut(full[col], 4, labels=False, duplicates="drop")
        metrics[f"by_{col}"] = full.groupby(f"{col}_quartile").agg(
            n=("CID", "size"), mean_diff=("diff", "mean"), mean_abs=("abs_diff", "mean")).to_dict("index")
    metrics["adequate"] = bool(metrics["rms"] <= .005 and metrics["p99_abs"] <= .02 and
                               metrics["spearman"] >= .995 and metrics["gt5"]["agreement"] >= .99 and
                               metrics["gt999"]["agreement"] >= .98)
    (OUT / "des_metrics.json").write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in metrics.items() if not k.startswith("by_")}, indent=2), flush=True)


def run_sim(model: str):
    metrics = json.loads((OUT / "des_metrics.json").read_text())
    assert metrics["adequate"], "DES reproduction adequacy gate failed"
    stem = ROOT / f"phase2/literature/simulations/outputs/PH2_pilot02_{model}/PH2_pilot02_{model}"
    df, meta = load_curves(stem)
    pred = infer(df)
    out = meta.merge(pred, on="CID", validate="one_to_one")
    out["gt999"] = out.pIa > .999
    out["gt5"] = out.pIa > .5
    out.sort_values("CID").to_csv(OUT / f"{model}_probabilities.csv", index=False)
    stats = {"variant": model, "n_detected_head": len(out),
             "n_raw": int(out.n_raw.sum()), "n_rejected": int(out.n_rejected.sum()),
             "n_retained": int(out.n_retained.sum()), "gt999": int(out.gt999.sum()),
             "gt5": int(out.gt5.sum())}
    out["zquartile"] = pd.qcut(out.redshift_final, 4, labels=False, duplicates="drop")
    stats["by_zquartile"] = out.groupby("zquartile").agg(n=("CID", "size"),
        gt999=("gt999", "sum"), gt5=("gt5", "sum")).to_dict("index")
    (OUT / f"{model}_metrics.json").write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    print(json.dumps(stats, indent=2), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("des", *VARIANTS))
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.mode == "des":
        run_des()
    else:
        run_sim(args.mode)
