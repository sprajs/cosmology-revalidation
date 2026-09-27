#!/usr/bin/env python3
"""Validate identities, numerical equations and observational result closure."""

from pathlib import Path
import json, datetime
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from scipy.optimize import minimize
from acquire import ROOT, WORK, OUT, sha, check_source
from query_sdss import ID_DTYPES
from spectra import pixel_weights
import diagnostics as diag


def main():
    lock = json.loads(Path(__file__).with_name("source-lock.json").read_text())["files"]
    for rel in lock:
        check_source(WORK / rel)
    d = pd.read_csv(WORK / "host-spectrum-age-ledger.csv", dtype=ID_DTYPES)
    assert len(d) == 810 and d.specObjID.nunique() == 794
    assert d.source_id.is_unique
    # Confirm ID identity against the original downloaded SQL response bytes.
    src = pd.concat(
        [
            pd.read_csv(p, comment="#", dtype=ID_DTYPES)
            for p in (WORK / "sdss-queries").glob("query-*.csv")
        ],
        ignore_index=True,
    ).drop_duplicates("specObjID")
    q = d[d.specObjID.notna()].merge(
        src[["specObjID", "bestObjID", "photometric_objID"]],
        on="specObjID",
        suffixes=("_saved", "_original"),
        validate="many_to_one",
    )
    for k in ["bestObjID", "photometric_objID"]:
        assert q[k + "_saved"].fillna("").equals(q[k + "_original"].fillna(""))
    # Verify original positive host match cuts and independently calculated ranks.
    assert d.position_z_accepted.all() and d.unique_photometric_counterpart.all()
    assert (d.separation_arcsec <= 3).all() and (abs(d.delta_z) < 0.001).all()
    r = json.loads((OUT / "diagnostics-summary.json").read_text())
    a = d[(d["sample"] == "ZTF") & d.index_quality & d.sed_age.notna()].drop_duplicates(
        "physical_host"
    )
    rho = float(spearmanr(a.sed_age, a.d4000_n).statistic)
    assert (
        abs(
            rho
            - r["age_index_relations"]["ZTF"]["associations"]["d4000_n"]["raw"][
                "spearman"
            ]
        )
        < 1e-12
    )
    # Flat Fnu and Flambda integral identities plus finite-noise variance check.
    w = np.linspace(3800, 4150, 35001)
    b = pixel_weights(w, 3850, 3950)
    c = pixel_weights(w, 4000, 4100)
    b /= b.sum()
    c /= c.sum()
    flatnu = float((c @ np.ones(len(w))) / (b @ np.ones(len(w))))
    flatlambda = float((c @ (w * w)) / (b @ (w * w)))
    exact = (4100**3 - 4000**3) / (3950**3 - 3850**3)
    assert abs(flatnu - 1) < 1e-12 and abs(flatlambda - exact) < 1e-8
    means = pd.read_csv(WORK / "remeasured-dn4000.csv")
    means = means[means.valid].copy()
    ratio_identity = float(
        np.max(abs(means.raw_red_Fnu_uJy / means.raw_blue_Fnu_uJy - means.raw_Dn4000))
    )
    assert ratio_identity < 1e-12
    # Direct likelihood reevaluation and independent optimizer on actual165SN fit.
    cohort = pd.read_csv(WORK / "brightness-cohort.csv", dtype=ID_DTYPES)
    cohort = diag.ztf.prepare(cohort)
    names = r["brightness"]["models"]["host_spectra_age"]
    fit = r["brightness"]["primary_fixed_indices_and_SED_medians"]["fits"][
        "host_spectra_age"
    ]
    x = diag.matrix(cohort, names)
    y = cohort.y.to_numpy()
    beta = np.array([fit["parameters"][k] for k in names])
    sigma = fit["intrinsic_sigma_mag"]

    def cost(v):
        var = diag.titan.variance(cohort, v[:-1], names, False) + np.exp(v[-1]) ** 2
        return float(np.sum((y - x @ v[:-1]) ** 2 / var + np.log(2 * np.pi * var)))

    point = np.r_[beta, np.log(sigma)]
    value = cost(point)
    assert abs(value - fit["minus2loglike"]) < 1e-9
    rng = np.random.default_rng(2026092801)
    other = minimize(
        cost,
        point + rng.normal(size=len(point)) * 0.02,
        method="SLSQP",
        options={"ftol": 1e-10, "maxiter": 1500},
    )
    difference = float(other.fun - value)
    assert abs(difference) < 1e-5
    # No host leaks into multiple predictive folds.
    folds = pd.read_csv(WORK / "brightness-heldout.csv")
    assert folds.groupby("physical_host").fold.nunique().max() == 1
    # Missing-value sentinels must not enter aperture-age differences.
    for sample in r["MaNGA_FIREFLY"]["samples"].values():
        for lib in ["miles", "mastar"]:
            q = sample[lib + "_central_minus_1Re_log10Gyr"]
            if q["n"]:
                assert max(abs(np.array(q["min_max"]))) < 5
    miles = pd.read_csv(WORK / "firefly-miles-selected-hosts.csv")
    mastar = pd.read_csv(WORK / "firefly-mastar-selected-hosts.csv")
    same = miles[["source_id", "ff_PLATEIFU"]].merge(
        mastar[["source_id", "ff_PLATEIFU"]],
        on="source_id",
        suffixes=("_a", "_b"),
        validate="one_to_one",
    )
    assert (same.ff_PLATEIFU_a == same.ff_PLATEIFU_b).all()
    m = json.loads((OUT / "manga-local-summary.json").read_text())
    ml = pd.read_csv(WORK / "manga-central-SN-indices.csv")
    good = ml.SN_site_valid
    assert (ml.loc[good, "SN_site_coordinate_error_arcsec"] < 0.36).all()
    assert (ml.loc[good, "SN_site_spaxel_SNR"] > 3).all()
    out = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "passed": True,
        "source_bytes_verified": len(lock),
        "exact_identifier_comparisons": int(d.specObjID.notna().sum()),
        "matched_host_rows": len(d),
        "distinct_spectra": d.specObjID.nunique(),
        "flat_Fnu_ratio": flatnu,
        "flat_Flambda_analytic_difference": flatlambda - exact,
        "observed_band_ratio_identity_max": ratio_identity,
        "independent_rank_correlation": rho,
        "actual_brightness_NLL_identity_difference": value - fit["minus2loglike"],
        "independent_SLSQP_minus_Minuit_NLL": difference,
        "SLSQP_status": str(other.message),
        "host_fold_leakage": False,
        "MaNGA_valid_SN_site_rows": int(good.sum()),
        "code_sha256": sha(Path(__file__)),
        "result_hashes": {
            p.name: sha(p)
            for p in [
                OUT / "diagnostics-summary.json",
                OUT / "manga-local-summary.json",
                OUT / "legacy-nnls-audit.json",
                OUT / "historical-age-summary.json",
            ]
        },
    }
    (OUT / "validation.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
