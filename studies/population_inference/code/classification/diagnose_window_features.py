#!/usr/bin/env python3
"""Compare corrected wrapper with official FITS-reader preprocessing on ten CIDs."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reconstruct import MASK, MODEL, OUT, RAW, load_curves
from supernnova.utils.data_utils import load_pandas_from_fit
from supernnova.validation.validate_onthefly import format_data, get_settings

IDS = ("1770702", "1929034", "1814082", "1315795", "1296161",
       "1291580", "1256426", "1257112", "1297365", "1262899")
FIELDS = ("MJD", "FLUXCAL", "FLUXCALERR", "HOSTGAL_SPECZ", "HOSTGAL_SPECZ_ERR")


def main():
    clump = pd.read_csv(OUT / "clump_run/DES-SN5YR_DES.SNANA.TEXT",
                        comment="#", delimiter=" ", skipinitialspace=True, dtype={"CID": str})
    peaks = dict(zip(clump.CID, clump.PKMJDINI.astype(float)))
    wrapper, _ = load_curves(RAW / "DES-SN5YR_DES", set(IDS), peaks)
    head = load_pandas_from_fit(RAW / "DES-SN5YR_DES_HEAD.FITS.gz")
    phot = load_pandas_from_fit(RAW / "DES-SN5YR_DES_PHOT.FITS.gz")
    if isinstance(head.SNID.iloc[0], bytes):
        head["SNID"] = head.SNID.str.decode("utf8").str.strip()
    else:
        head["SNID"] = head.SNID.astype(str).str.strip()
    head = head.set_index("SNID")
    frames = []
    for cid in IDS:
        h = head.loc[cid]
        q = phot.iloc[int(h.PTROBS_MIN) - 1:int(h.PTROBS_MAX)].copy()
        q["SNID"] = cid
        q["FLT"] = q.BAND.apply(lambda x: x.rstrip()).values.astype(str)
        dt = q.MJD - peaks[cid]
        q = q[(dt > -30) & (dt < 100)]
        q = q[(q.PHOTFLAG.astype(np.int64) & MASK) == 0]
        q["HOSTGAL_SPECZ"] = h.REDSHIFT_FINAL
        q["HOSTGAL_SPECZ_ERR"] = h.REDSHIFT_FINAL_ERR
        q["HOSTGAL_PHOTOZ"] = h.HOSTGAL_PHOTOZ
        q["HOSTGAL_PHOTOZ_ERR"] = h.HOSTGAL_PHOTOZ_ERR
        q["PEAKMJD"] = h.PEAKMJD
        q["SNTYPE"] = h.SNTYPE
        q["SIM_REDSHIFT_CMB"] = h.REDSHIFT_FINAL
        frames.append(q)
    official = pd.concat(frames, ignore_index=True)
    settings = get_settings(str(MODEL))
    a = format_data(wrapper.copy(), settings)
    b = format_data(official.copy(), settings)
    model_features = [f for f in settings.all_features if f in settings.training_features]
    result = {}
    for cid in IDS:
        w = wrapper.loc[wrapper.SNID == cid].reset_index(drop=True)
        o = official.loc[official.SNID == cid].reset_index(drop=True)
        assert len(w) == len(o)
        assert (w.FLT.to_numpy() == o.FLT.to_numpy()).all()
        rowdiff = max(float(np.max(np.abs(w[f].to_numpy(dtype=float) - o[f].to_numpy(dtype=float))))
                      for f in FIELDS)
        x = a.loc[cid, model_features].to_numpy(dtype=float).reshape(-1, len(model_features))
        y = b.loc[cid, model_features].to_numpy(dtype=float).reshape(-1, len(model_features))
        assert x.shape == y.shape
        result[cid] = {"kind": "largest_corrected_mismatch" if IDS.index(cid) < 5 else "near_corrected_agreement",
                       "retained_rows": len(w), "pivot_rows": len(x),
                       "max_abs_raw_row_difference": rowdiff,
                       "max_abs_model_feature_difference": float(np.max(np.abs(x-y))),
                       "excluded_photoz_difference": float(abs(a.loc[cid, "HOSTGAL_PHOTOZ"].iloc[0] -
                                                                b.loc[cid, "HOSTGAL_PHOTOZ"].iloc[0]))}
    (OUT / "diagnostic_window_features.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
