#!/usr/bin/env python3
"""Separate DES validation using Pippin-style test database normalization."""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reconstruct import MODEL, OUT, RAW, REF, file_sha, infer, load_curves
from supernnova.validation.validate_onthefly import format_data, get_settings

CLUMP = OUT / "clump_run/DES-SN5YR_DES.SNANA.TEXT"
HDF5 = OUT / "test_database/processed/database.h5"
PLAN = OUT / "test_norm_plan.md"
TEN = ("1770702", "1929034", "1814082", "1315795", "1296161",
       "1291580", "1256426", "1257112", "1297365", "1262899")


def setup_model() -> Path:
    model_dir = OUT / "model_with_DES_test_norm"
    model_dir.mkdir(exist_ok=True)
    shutil.copyfile(MODEL.parent / "cli_args.json", model_dir / "cli_args.json")
    target = model_dir / "model.pt"
    if not target.exists():
        target.symlink_to(MODEL)
    assert target.resolve() == MODEL
    with h5py.File(HDF5, "r") as f:
        norm = {}
        for name in ("FLUXCAL_g", "FLUXCAL_i", "FLUXCAL_r", "FLUXCAL_z",
                     "FLUXCALERR_g", "FLUXCALERR_i", "FLUXCALERR_r", "FLUXCALERR_z",
                     "delta_time"):
            prefix = name.split("_")[0]
            path = f"normalizations_global/{prefix}" if name != "delta_time" else "normalizations/delta_time"
            norm[name] = {key: float(f[f"{path}/{key}"][()]) for key in ("min", "mean", "std")}
    (model_dir / "data_norm.json").write_text(json.dumps(norm, indent=2, sort_keys=True) + "\n")
    return target


def validate_hdf5(df: pd.DataFrame, model: Path):
    s = get_settings(str(model))
    fmt = format_data(df[df.SNID.isin(TEN)].copy(), s)
    selected = [f for f in s.all_features if f in s.training_features]
    result = {}
    with h5py.File(HDF5, "r") as f:
        ids = f["SNID"][:].astype(str)
        positions = {cid: i for i, cid in enumerate(ids)}
        features = f["features"][:].astype(str).tolist()
        nfeat = int(f["data"].attrs["n_features"])
        assert features == s.all_features and nfeat == len(features)
        for cid in TEN:
            assert cid in positions
            h = f["data"][positions[cid]].reshape(-1, nfeat)
            a = h[:, [features.index(name) for name in selected]]
            b = fmt.loc[cid, selected].to_numpy(dtype=np.float32).reshape(-1, len(selected))
            assert a.shape == b.shape
            result[cid] = {"shape": list(a.shape), "max_abs_raw_feature_difference": float(np.max(np.abs(a-b)))}
    return result


def main():
    plan_hash = file_sha(PLAN)
    model = setup_model()
    peaks_df = pd.read_csv(CLUMP, comment="#", delimiter=" ", skipinitialspace=True,
                           dtype={"CID": str})
    peaks = dict(zip(peaks_df.CID, peaks_df.PKMJDINI.astype(float)))
    ref = pd.read_csv(REF, dtype={"CID": str})
    ids = set(ref.CID)
    df, meta = load_curves(RAW / "DES-SN5YR_DES", ids, peaks)
    fixture = validate_hdf5(df, model)
    assert all(v["max_abs_raw_feature_difference"] <= 1e-5 for v in fixture.values()), fixture
    pred = infer(df, model_path=model)
    out = meta.merge(pred, on="CID", validate="one_to_one").merge(
        ref[["CID", "PROB_SNNV19"]], on="CID", validate="one_to_one")
    prev = pd.read_csv(OUT / "des_clump_diagnostic.csv", dtype={"CID": str})
    out = out.merge(prev[["CID", "pIa"]].rename(columns={"pIa": "clump_training_norm_pIa"}),
                    on="CID", validate="one_to_one")
    out["diff"] = out.pIa - out.PROB_SNNV19
    out["abs_diff"] = out["diff"].abs()
    out.sort_values("CID").to_csv(OUT / "des_test_norm_diagnostic.csv", index=False)
    err = out["diff"].to_numpy()
    ab = np.abs(err)
    metrics = {"status": "diagnostic source-derived DES test normalization; historical executable not exact",
               "plan_sha256": plan_hash, "hdf5_sha256": file_sha(HDF5),
               "model_original_sha256": file_sha(MODEL),
               "model_test_norm_sha256": file_sha(model.parent / "data_norm.json"),
               "n": len(out), "n_retained": int(meta.n_retained.sum()),
               "median_abs": float(np.median(ab)), "p95_abs": float(np.quantile(ab, .95)),
               "p99_abs": float(np.quantile(ab, .99)), "max_abs": float(ab.max()),
               "rms": float(np.sqrt(np.mean(err*err))), "bias": float(np.mean(err)),
               "spearman": float(spearmanr(out.pIa, out.PROB_SNNV19).statistic),
               "ten_hdf5_feature_checks": fixture}
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
    (OUT / "des_test_norm_diagnostic_metrics.json").write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k:v for k,v in metrics.items() if k != "ten_hdf5_feature_checks"}, indent=2))


if __name__ == "__main__":
    main()
