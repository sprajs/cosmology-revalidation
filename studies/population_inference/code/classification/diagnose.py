#!/usr/bin/env python3
"""Bounded feature-path comparison after the preserved failed DES gate."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reconstruct import MODEL, OUT, RAW, load_curves

from supernnova.validation.validate_onthefly import classify_lcs, format_data, get_settings
from supernnova.utils.training_utils import get_model, normalize_arr
from supernnova.data.make_dataset import pivot_dataframe_single_from_df
from supernnova.utils.data_utils import compute_delta_time

IDS = ("1343533", "1919017", "1335323", "1344639", "1265165",
       "1402076", "1297365", "1256426", "1257112", "1380034")


def feature_arrays(df, settings):
    fmt = format_data(df.copy(), settings)
    idx_norm = [i for i, f in enumerate(settings.training_features)
                if f in settings.training_features_to_normalize]
    idx_selected = [i for i, f in enumerate(settings.all_features)
                    if f in settings.training_features]
    settings.idx_features_to_normalize = idx_norm
    out = {}
    for cid in IDS:
        arr = fmt.loc[cid, settings.all_features].to_numpy().reshape(-1, len(settings.all_features))
        arr = np.asarray(arr, dtype=np.float64)
        arr[:, idx_norm] = np.clip(arr[:, idx_norm], settings.arr_norm[:, 0], np.inf)
        out[cid] = normalize_arr(arr.copy(), settings)[:, idx_selected]
    return fmt, out, idx_norm, idx_selected


def dataset_path_arrays(df, settings, idx_norm, idx_selected):
    # The official make_dataset sequence, separate from onthefly.format_data:
    # compute_delta_time -> pivot -> save_to_HDF5 onehot -> normalize/select.
    pivot = pivot_dataframe_single_from_df(compute_delta_time(df.copy()), settings)
    base = settings.all_features[:13]
    tmp = pd.Series(settings.list_filters_combination).append(pivot["FLT"])
    dummies = pd.get_dummies(tmp)[len(settings.list_filters_combination):]
    combined = pd.concat([pivot[base], dummies], axis=1)[settings.all_features]
    out = {}
    for cid in IDS:
        arr = combined.loc[cid, settings.all_features].to_numpy().reshape(-1, len(settings.all_features))
        arr = np.asarray(arr, dtype=np.float64)
        arr[:, idx_norm] = np.clip(arr[:, idx_norm], settings.arr_norm[:, 0], np.inf)
        out[cid] = normalize_arr(arr.copy(), settings)[:, idx_selected]
    return out


def main():
    minimal, meta = load_curves(RAW / "DES-SN5YR_DES", set(IDS))
    with fits.open(str(RAW / "DES-SN5YR_DES_HEAD.FITS.gz"), memmap=False) as f:
        head = f[1].data.copy()
    hh = {str(h["SNID"]).strip(): h for h in head if str(h["SNID"]).strip() in IDS}
    full = minimal.copy()
    for col in ("HOSTGAL_PHOTOZ", "HOSTGAL_PHOTOZ_ERR", "PEAKMJD", "SNTYPE"):
        full[col] = [hh[cid][col] for cid in full.SNID]
    full["SIM_REDSHIFT_CMB"] = full["HOSTGAL_SPECZ"]
    s = get_settings(str(MODEL))
    mfmt, marr, inorm, iselect = feature_arrays(minimal, s)
    ffmt, farr, _, _ = feature_arrays(full, s)
    darr = dataset_path_arrays(full, s, inorm, iselect)
    # Independent direct model evaluation: pack normalized feature tensors here,
    # bypassing onthefly's batch construction and prediction extraction.
    import torch
    order = sorted(IDS, key=lambda cid: (-len(farr[cid]), cid))
    packed = torch.nn.utils.rnn.pack_sequence(
        [torch.as_tensor(farr[cid], dtype=torch.float32) for cid in order],
        enforce_sorted=True)
    s.use_cuda = False
    s.device = "cpu"
    torch.manual_seed(s.seed)
    rnn = get_model(s, len(iselect))
    rnn.load_state_dict(torch.load(str(MODEL), map_location="cpu"))
    rnn.eval()
    with torch.no_grad():
        direct = torch.nn.functional.softmax(rnn(packed), dim=-1)[:, 0].numpy()
    direct = dict(zip(order, map(float, direct)))
    official_ids, official_pred = classify_lcs(minimal.copy(), str(MODEL), "cpu")
    official = dict(zip(map(str, official_ids), map(float, official_pred[:, 0, 0])))
    original_z = minimal.copy()
    original_z["HOSTGAL_SPECZ"] = [float(hh[cid]["HOSTGAL_SPECZ"]) for cid in original_z.SNID]
    original_z["HOSTGAL_SPECZ_ERR"] = [float(hh[cid]["HOSTGAL_SPECZ_ERR"]) for cid in original_z.SNID]
    original_ids, original_pred = classify_lcs(original_z, str(MODEL), "cpu")
    original = dict(zip(map(str, original_ids), map(float, original_pred[:, 0, 0])))
    per = {}
    for cid in IDS:
        a, b = marr[cid], farr[cid]
        assert a.shape == b.shape, (cid, a.shape, b.shape)
        diff = np.abs(a - b)
        dataset_diff = np.abs(darr[cid] - b)
        per[cid] = {"kind": "largest_mismatch" if IDS.index(cid) < 5 else "near_agreement",
                    "raw_rows": int(meta.set_index("CID").loc[cid, "n_raw"]),
                    "retained_rows": int(meta.set_index("CID").loc[cid, "n_retained"]),
                    "pivot_rows_minimal": int(len(mfmt.loc[cid])),
                    "pivot_rows_full": int(len(ffmt.loc[cid])),
                    "tensor_shape": list(a.shape),
                    "max_abs_feature_difference": float(diff.max()),
                    "max_abs_dataset_onthefly_difference": float(dataset_diff.max()),
                    "direct_pIa": direct[cid], "onthefly_pIa": official[cid],
                    "original_HOSTGAL_SPECZ_pIa_diagnostic": original[cid],
                    "redshift_final": float(hh[cid]["REDSHIFT_FINAL"]),
                    "original_HOSTGAL_SPECZ": float(hh[cid]["HOSTGAL_SPECZ"]),
                    "direct_onthefly_abs_difference": abs(direct[cid] - official[cid]),
                    "first_differing_feature": None if not np.any(diff) else
                        s.all_features[iselect[np.where(np.any(diff != 0, axis=0))[0][0]]]}
    state = torch.load(str(MODEL), map_location="cpu")
    info = {"ids": list(IDS), "per_object": per,
            "all_features": s.all_features, "training_features_saved": s.training_features,
            "actual_selected_features": [s.all_features[i] for i in iselect],
            "idx_features_to_normalize": inorm, "idx_selected": iselect,
            "norm": s.norm, "arr_norm": s.arr_norm.tolist(),
            "model_weight_ih_l0_shape": list(state["rnn_layer.weight_ih_l0"].shape),
            "class_0": "Ia (tag_type in official data_utils.py)",
            "redshift": "zspe; HOSTGAL_SPECZ[_ERR] set from REDSHIFT_FINAL[_ERR]",
            "time_origin": "official compute_delta_time first retained MJD per SNID; pivot recomputes grouped deltas"}
    (OUT / "diagnostic_features.json").write_text(json.dumps(info, indent=2) + "\n")
    print(json.dumps({"per_object": per, "actual_selected_features": info["actual_selected_features"],
                      "model_weight_ih_l0_shape": info["model_weight_ih_l0_shape"]}, indent=2))


if __name__ == "__main__":
    main()
