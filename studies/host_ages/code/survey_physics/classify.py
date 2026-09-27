#!/usr/bin/env python3
"""Reconstructed SNNV19 classification with measured peak/redshift inputs only.

Run under the preserved Python3.10/SuperNNova environment documented in README.
This is conditional inference, not training or proof of original probabilities.
"""
import argparse
import json
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from astropy.io import fits
from supernnova.validation.validate_onthefly import classify_lcs
from common import ROOT, WORK, RESULTS, DATA, ARCHIVE, sha

sys.path.insert(0, str(ROOT))
from lib.records import fitres


def load_curves(stem, ids, peaks):
    with fits.open(str(stem) + "_HEAD.FITS") as h:
        head = h[1].data.copy()
    with fits.open(str(stem) + "_PHOT.FITS") as h:
        phot = h[1].data.copy()
    names = [str(v).strip() for v in head["SNID"]]
    assert len(names) == len(set(names))
    aliases = {"DES-" + b: b for b in "griz"}
    aliases.update({b: b for b in "griz"})
    rows, metadata = [], []
    for h, name in zip(head, names):
        if name not in ids:
            continue
        lo, hi = int(h["PTROBS_MIN"]), int(h["PTROBS_MAX"])
        assert 1 <= lo <= hi <= len(phot) and hi - lo + 1 == int(h["NOBS"])
        q = phot[lo - 1 : hi]
        assert np.all(q["MJD"] > 0)
        delta = q["MJD"] - peaks[name]
        window = (delta > -30) & (delta < 100)
        q = q[window]
        flags = (q["PHOTFLAG"] & 1016) == 0
        q = q[flags]
        assert len(q) > 0
        for point in q:
            band = str(point["BAND"]).strip()
            assert band in aliases, band
            rows.append(
                {
                    "SNID": name,
                    "MJD": float(point["MJD"]),
                    "FLT": aliases[band],
                    "FLUXCAL": float(point["FLUXCAL"]),
                    "FLUXCALERR": float(point["FLUXCALERR"]),
                    "HOSTGAL_SPECZ": float(h["REDSHIFT_FINAL"]),
                    "HOSTGAL_SPECZ_ERR": float(h["REDSHIFT_FINAL_ERR"]),
                }
            )
        metadata.append(
            {
                "CID": name,
                "n_raw": hi - lo + 1,
                "n_window_rejected": int((~window).sum()),
                "n_flag_rejected": int((~flags).sum()),
                "n_retained": len(q),
                "peak_mjd_measured": peaks[name],
            }
        )
    assert {m["CID"] for m in metadata} == ids
    frame = pd.DataFrame(rows)
    assert (
        np.isfinite(
            frame[
                ["MJD", "FLUXCAL", "FLUXCALERR", "HOSTGAL_SPECZ", "HOSTGAL_SPECZ_ERR"]
            ]
        )
        .all()
        .all()
    )
    return frame, pd.DataFrame(metadata)


def infer(frame, model, chunk_size=128):
    ids = sorted(frame.SNID.unique())
    result = []
    for start in range(0, len(ids), chunk_size):
        selected = set(ids[start : start + chunk_size])
        names, prediction = classify_lcs(
            frame[frame.SNID.isin(selected)].copy(), str(model), "cpu"
        )
        assert prediction.shape == (len(selected), 1, 2)
        result.extend((str(n), float(v[0, 0])) for n, v in zip(names, prediction))
    output = pd.DataFrame(result, columns=["CID", "pIa"])
    assert (
        output.CID.is_unique
        and len(output) == len(ids)
        and output.pIa.between(0, 1).all()
    )
    return output


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--name", required=True)
    args = p.parse_args()
    folder = WORK / "fits" / args.name
    run = json.loads((folder / "run.json").read_text())
    assert run["fit_graceful"], "Wait for completed native fit."
    source = Path(sys.modules[classify_lcs.__module__].__file__)
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    model = DATA / "models/classifiers/DES-SN5YR/SNNTRAINV19_z_TRAINDES_V19/model.pt"
    fitpath = folder / "fit.FITRES.TEXT"
    peakpath = folder / "fit.SNANA.TEXT"
    fitted = fitres(fitpath)
    peaks = fitres(peakpath)
    assert fitted.index.is_unique and peaks.index.is_unique
    peak = peaks.PKMJDINI.loc[fitted.index].to_dict()
    assert all(np.isfinite(list(peak.values())))
    stem = WORK / "simulations" / args.name / args.name
    frame, meta = load_curves(stem, set(fitted.index), peak)
    start = time.monotonic()
    predicted = infer(frame, model)
    elapsed = time.monotonic() - start
    # Deterministic batch-versus-single check on a prespecified evenly spaced subset.
    checkids = list(sorted(fitted.index))[:: max(1, len(fitted) // 5)][:5]
    singles = infer(frame[frame.SNID.isin(checkids)], model, 1)
    comparison = (
        predicted.set_index("CID")
        .join(singles.set_index("CID"), lsuffix="_batch", rsuffix="_single")
        .dropna()
    )
    delta = float(np.max(abs(comparison.pIa_batch - comparison.pIa_single)))
    assert delta < 1e-6
    predicted.to_csv(folder / "classifier.csv", index=False)
    meta.to_csv(folder / "classifier-preprocessing.csv", index=False)
    result = {
        "name": args.name,
        "code_sha256": sha(__file__),
        "wrapper_sha256": sha(source),
        "exact_band_aliases": {"DES-" + b: b for b in "griz"},
        "model_sha256": sha(model),
        "native_fit_sha256": sha(fitpath),
        "native_peak_sha256": sha(peakpath),
        "classified": len(predicted),
        "pIa_gt05": int((predicted.pIa > 0.5).sum()),
        "pIa_gt0999": int((predicted.pIa > 0.999).sum()),
        "inference_seconds": elapsed,
        "device": "CPU",
        "torch_version": torch.__version__,
        "single_batch_max_probability_difference": delta,
        "truth_inputs_used": False,
        "original_classifier_equivalence_certified": False,
        "output_sha256": sha(folder / "classifier.csv"),
        "preprocessing_sha256": sha(folder / "classifier-preprocessing.csv"),
    }
    (folder / "classifier.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
