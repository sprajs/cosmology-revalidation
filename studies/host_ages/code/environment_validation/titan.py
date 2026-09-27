#!/usr/bin/env python3
"""Common-sample ZTF brightness comparison with public TITAN host summaries.

Age and delay are model-inferred global host quantities. All primary results
condition on posterior medians. The diagonal-error sensitivity is not a joint
age/dust likelihood and cannot identify a causal cosmological correction.
"""
import argparse, datetime, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm
import ztf
from analyse import ROOT, HERE, DEFAULT_ARCHIVE, DESIGN, sha

TD = json.loads((HERE / "titan-design.json").read_text())
BASE = ztf.MODELS["nonlinear"]
CONTROL = BASE + ["titan_mass", "titan_AV", "titan_metallicity"]
MODELS = {
    "standard": BASE,
    "host_controls": CONTROL,
    "host_controls_age": CONTROL + ["host_age"],
    "host_controls_delay": CONTROL + ["delay"],
}
FIELDS = {
    "titan_mass": "stellar_mass",
    "titan_AV": "dust:Av",
    "titan_metallicity": "iyer2019:metallicity",
    "host_age": "mass_weighted_age",
    "delay": "prog_age",
}
CENTRE = {
    "titan_mass": 10,
    "titan_AV": 0.3,
    "titan_metallicity": 1,
    "host_age": 5,
    "delay": 3,
}


def matrix(d, names):
    return np.column_stack(
        [
            (
                d[FIELDS[k] + "_50"].to_numpy() - CENTRE[k]
                if k in FIELDS
                else ztf.make_matrix(d, [k])[:, 0]
            )
            for k in names
        ]
    )


def variance(d, b, names, errors=False):
    v = ztf.variance(d, b, names)
    if errors:
        for k, value in zip(names, b):
            if k in FIELDS:
                f = FIELDS[k]
                v += (
                    value**2
                    * ((d[f + "_84"].to_numpy() - d[f + "_16"].to_numpy()) / 2) ** 2
                )
    return v


def run(d, cv=True, errors=False):
    varfun = lambda dd, bb, nn: variance(dd, bb, nn, errors)
    fits = {}
    scores = {}
    ledger = d[["ztfname", "iau_name", "host_group", "fold"]].copy()
    for name, cols in MODELS.items():
        f, b, bc, ss = ztf.fit(d, cols, matrix, varfun)
        fits[name] = f
        if cv:
            score = np.zeros(len(d))
            pred = np.zeros(len(d))
            vpred = np.zeros(len(d))
            for fold in range(5):
                tr = d.fold.to_numpy() != fold
                te = ~tr
                ff, bb, cc, sc = ztf.fit(d[tr], cols, matrix, varfun)
                xt = matrix(d[te], cols)
                pred[te] = xt @ bb
                vpred[te] = (
                    varfun(d[te], bb, cols)
                    + sc**2
                    + np.einsum("ni,ij,nj->n", xt, cc, xt)
                )
                delta = d.y.to_numpy()[te] - pred[te]
                score[te] = -0.5 * (
                    np.log(2 * np.pi * vpred[te]) + delta**2 / vpred[te]
                )
            scores[name] = score
            ledger[name + "_score"] = score
            f["heldout_log_predictive_score"] = float(score.sum())
            f["heldout_RMSE_mag"] = float(
                np.sqrt(np.mean((d.y.to_numpy() - pred) ** 2))
            )
    pairs = {}
    rng = np.random.default_rng(DESIGN["seed"])
    groups = d.host_group.to_numpy()
    for added, base in [
        ("host_controls", "standard"),
        ("host_controls_age", "host_controls"),
        ("host_controls_delay", "host_controls"),
    ]:
        if cv:
            delta = scores[added] - scores[base]
            units = np.array([delta[groups == g].sum() for g in sorted(set(groups))])
            boots = rng.choice(units, (2000, len(units)), replace=True).sum(axis=1)
            pairs[added + "_minus_" + base] = {
                "delta_log_predictive_score": float(delta.sum()),
                "host_bootstrap_95_conditional_on_fitted_folds": np.quantile(
                    boots, [0.025, 0.975]
                ).tolist(),
            }
    for key in ["host_age", "delay"]:
        name = "host_controls_age" if key == "host_age" else "host_controls_delay"
        f = fits[name]
        b = f["parameters"][key]
        se = f["standard_errors"][key]
        delta = 0.03 / se
        f["age_or_delay_slope_95"] = [b - 1.96 * se, b + 1.96 * se]
        f["normal_approximation_power_for_abs_003"] = float(
            norm.cdf(-1.96 - delta) + norm.sf(1.96 - delta)
        )
    return {
        "n": len(d),
        "hosts": int(d.host_group.nunique()),
        "age_error_treatment": (
            "independent marginal Gaussian propagation"
            if errors
            else "fixed posterior medians"
        ),
        "fits": fits,
        "predictive_comparisons": pairs,
    }, ledger


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    ap.add_argument(
        "--titan",
        type=Path,
        default=ROOT
        / ".work/population-transport/titan-author/host_props_with_SN_age_good_Mar18.csv",
    )
    ap.add_argument("--work", type=Path, default=ROOT / ".work/environment-validation")
    ap.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / "studies/host_ages/results/environment_validation/titan-ztf-summary.json",
    )
    args = ap.parse_args()
    assert sha(args.titan) == TD["input_sha256"]
    t = pd.read_csv(args.titan)
    assert t.transient.is_unique
    base = args.archive / "sources/updates/2026-09-20-ztf"
    paths = [
        base / "extracted/ztfsniadr2_lite/tables" / k
        for k in ["snia_data.csv", "globalhost_data.csv", "localhost_data.csv"]
    ] + [base / "originals/Ginolin25ab_masterlist--903d65d.csv"]
    sn, glo, loc, master = [pd.read_csv(p) for p in paths]
    d, _, _ = ztf.sample(sn, glo, loc, master)
    d = ztf.prepare(d)
    # Exact IAU IDs link physical SNe. TITAN has no sky coordinates here, so its
    # host identity cannot be independently coordinate-verified in this file.
    q = d.merge(
        t,
        left_on="iau_name",
        right_on="transient",
        suffixes=("", "_titan"),
        validate="one_to_one",
    )
    matches = len(q)
    cut = (q.d_dlr < 4) & (q.d_dlr_titan < 4)
    cols = [f + "_" + s for f in FIELDS.values() for s in ["16", "50", "84"]]
    cut &= np.isfinite(q[cols]).all(axis=1)
    for f in FIELDS.values():
        cut &= (q[f + "_84"] >= q[f + "_50"]) & (q[f + "_50"] >= q[f + "_16"])
    q["selected_titan"] = cut
    q.drop(columns=["lc_cov"]).to_csv(args.work / "titan-ztf-join.csv", index=False)
    q = q[cut].reset_index(drop=True)
    assert len(q) > len(CONTROL) + 10
    primary, ledger = run(q)
    ledger.to_csv(args.work / "titan-ztf-heldout.csv", index=False)
    sensitivity, _ = run(q, errors=True)
    x = matrix(q, CONTROL)
    age = q.mass_weighted_age_50.to_numpy()
    delay = q.prog_age_50.to_numpy()
    age_res = age - x @ np.linalg.lstsq(x, age, rcond=None)[0]
    delay_res = delay - x @ np.linalg.lstsq(x, delay, rcond=None)[0]
    result = {
        "schema": 1,
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "code_sha256": sha(__file__),
        "dependencies_sha256": {
            p.name: sha(p) for p in [HERE / "analyse.py", HERE / "ztf.py"]
        },
        "design_sha256": sha(HERE / "titan-design.json"),
        "inputs": [{"path": str(p), "sha256": sha(p)} for p in [args.titan] + paths],
        "titan_snapshot_rows": len(t),
        "ztf_primary_sample": len(d),
        "exact_name_matches": matches,
        "selected_rows": len(q),
        "excluded_after_host_distance_and_finite_summary_cuts": matches - len(q),
        "host_identity_limit": "exact unique IAU alias matches; TITAN file has no coordinates, so host assignment is conditioned on author association and d_DLR<4 in both releases",
        "primary": primary,
        "independent_summary_error_sensitivity": sensitivity,
        "z_at_least_002_sensitivity": run(
            q[q.redshift >= 0.02].reset_index(drop=True), cv=False
        )[0],
        "colour_below_01_sensitivity": run(
            q[q.c < 0.1].reset_index(drop=True), cv=False
        )[0],
        "support": {
            "host_age_range_Gyr": [float(age.min()), float(age.max())],
            "host_age_sd_Gyr": float(np.std(age)),
            "host_age_sd_after_controls_Gyr": float(np.std(age_res)),
            "host_age_variance_remaining_fraction": float(
                np.var(age_res) / np.var(age)
            ),
            "delay_sd_after_controls_Gyr": float(np.std(delay_res)),
            "delay_variance_remaining_fraction": float(
                np.var(delay_res) / np.var(delay)
            ),
            "median_host_age_half_68_width_Gyr": float(
                np.median((q.mass_weighted_age_84 - q.mass_weighted_age_16) / 2)
            ),
            "median_delay_half_68_width_Gyr": float(
                np.median((q.prog_age_84 - q.prog_age_16) / 2)
            ),
        },
        "limitations": [
            "The age-of-titans CSV is a public prepaper 8610-row snapshot; the exact 6983 final-paper sample is not reproduced.",
            "Age, dust, mass, metallicity and inferred delay share SED data/model/prior dependence; their joint posterior covariance and full SFHs are not supplied here.",
            "Host age is global mass-weighted stellar age; prog_age is a posterior summary of an SFH-specific mean delay, not a direct age of each progenitor.",
            "Host attenuation is not line-of-sight SN extinction.",
            "The age coefficient is incremental conditional prediction after potential mediators, not total causal age dependence.",
            "Existing selected ZTF photometry is not cosmology-ready; no survey calibration, selection retraining or transport validation is supplied by this comparison.",
            "Cross-validation tests new SNe conditional on shared training/calibration; bootstrap intervals condition on fitted overlapping folds.",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(args.output)
    print(
        json.dumps(
            {"sample": len(q), "primary": primary, "support": result["support"]},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
