#!/usr/bin/env python3
"""Compare independent official HDF5 validation with windowed onthefly DES output."""
from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reconstruct import OUT, file_sha


def main():
    paths = list((OUT / "test_database/models").glob("*/PRED_*.pickle"))
    assert len(paths) == 1, paths
    with paths[0].open("rb") as f:
        official = pickle.load(f)
    official = official[["SNID", "all_class0"]].rename(
        columns={"SNID": "CID", "all_class0": "official_validate_pIa"})
    official.CID = official.CID.astype(str)
    assert official.CID.is_unique

    windowed = pd.read_csv(OUT / "des_clump_diagnostic.csv", dtype={"CID": str})
    out = windowed[["CID", "pIa", "PROB_SNNV19"]].merge(
        official, on="CID", how="left", validate="one_to_one")
    assert len(out) == 1635 and out.official_validate_pIa.notna().all()
    out["validate_minus_onthefly"] = out.official_validate_pIa - out.pIa
    out["validate_minus_release"] = out.official_validate_pIa - out.PROB_SNNV19
    out.sort_values("CID").to_csv(OUT / "official_validate_comparison.csv", index=False)
    delta = out.validate_minus_onthefly.to_numpy()
    release_delta = out.validate_minus_release.to_numpy()
    metrics = {
        "status": "independent pinned SuperNNova HDF5 validate_rnn diagnostic",
        "official_pickle_sha256": file_sha(paths[0]),
        "n_official_total": len(official),
        "n_des_classification": len(out),
        "max_abs_validate_minus_onthefly": float(np.max(np.abs(delta))),
        "rms_validate_minus_onthefly": float(np.sqrt(np.mean(delta**2))),
        "n_validate_minus_onthefly_gt_1e_6": int(np.sum(np.abs(delta) > 1e-6)),
        "rms_validate_minus_release": float(np.sqrt(np.mean(release_delta**2))),
    }
    (OUT / "official_validate_comparison.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    print(json.dumps(metrics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
