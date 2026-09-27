#!/usr/bin/env python3
"""Independent training/evaluation correction response for declared physical mocks."""
import io
import argparse
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from extinction import fitzpatrick99
from common import ROOT, HERE, WORK, RESULTS, sha

sys.path.insert(0, str(ROOT))
from lib.records import fitres

EDGES = np.array([0.05, 0.3, 0.5, 0.7, 0.9, 1.2])
NAMES = [
    "SPH_TRAIN_NOMINAL",
    "SPH_EVAL_NOMINAL",
    "SPH_EVAL_AGE",
    "SPH_EVAL_DUST",
    "SPH_TRAIN_AGE",
]


def dump(name):
    path = WORK / "simulations" / name / (name + ".DUMP")
    text = path.read_text().splitlines()
    header = next(x.split()[1:] for x in text if x.startswith("VARNAMES:"))
    frame = pd.read_csv(
        io.StringIO("\n".join(x[3:] for x in text if x.startswith("SN:"))),
        sep=r"\s+",
        names=header,
    )
    assert not frame.duplicated(["CID", "LIBID"]).any()
    return frame


def load(name):
    path = WORK / "fits" / name / "fit.FITRES.TEXT"
    f = fitres(path)
    assert f.index.is_unique
    prob = pd.read_csv(path.parent / "classifier.csv", dtype={"CID": str}).set_index(
        "CID"
    )
    f = f.join(prob, validate="one_to_one")
    assert f.pIa.notna().all()
    d = dump(name)
    campaign = (
        "positive-dust-campaign.json" if name.startswith("SPP_") else "campaign.json"
    )
    job = next(
        j
        for j in json.loads((RESULTS / campaign).read_text())["jobs"]
        if j["name"] == name
    )
    assert len(d) == job["attempts"]
    x = f[["x0ERR", "x1ERR", "cERR"]].to_numpy()
    cov = np.zeros((len(f), 3, 3))
    cov[:, range(3), range(3)] = x * x
    for i, j, k in [(0, 1, "COV_x1_x0"), (0, 2, "COV_c_x0"), (1, 2, "COV_x1_c")]:
        cov[:, i, j] = cov[:, j, i] = f[k]
    # Correlation scaling prevents declaring an x0 covariance PSD merely by tolerance.
    corr = cov / x[:, :, None] / x[:, None, :]
    finite = np.isfinite(corr).all(axis=(1, 2))
    eig = np.full(len(f), -np.inf)
    eig[finite] = np.linalg.eigvalsh(corr[finite])[:, 0]
    jac = np.stack(
        [-2.5 / np.log(10) / f.x0, np.full(len(f), 0.15), np.full(len(f), -3.14)],
        axis=1,
    )
    var = np.einsum("ni,nij,nj->n", jac, cov, jac)
    f["var_mu"] = var + 0.1**2
    f["covariance_positive"] = eig > 0
    f["quality"] = (
        f.ERRFLAG_FIT.eq(0)
        & f.x1.between(-3, 3)
        & f.c.between(-0.3, 0.3)
        & f.x1ERR.between(0, 1, inclusive="neither")
        & f.cERR.between(0, 1.5, inclusive="neither")
        & f.PKMJDERR.between(0, 2, inclusive="neither")
        & f.FITPROB.gt(0.001)
        & f.zHD.between(0.05, 1.2)
        & f.covariance_positive
        & np.isfinite(var)
        & (var > 0)
        & f.x0.gt(0)
    )
    wave = np.linspace(3500, 8000, 91).astype(float)
    amin = np.array(
        [
            np.min(fitzpatrick99(wave, float(a), float(r)))
            for a, r in zip(f.SIM_AV, f.SIM_RV)
        ]
    )
    f["physical_screen"] = amin >= -1e-10
    f["y"] = f.mB - f.SIM_DLMAG
    f["mass_step"] = np.where(f.HOST_LOGMASS >= 10, 0.5, -0.5)
    f["age"] = f.SIM_HOSTLIB_SN_age
    f["block"] = f.SIM_LIBID.astype(int)
    f["bin"] = np.clip(np.digitize(f.zHD, EDGES) - 1, 0, 4)
    quality = f.quality
    counts = {
        "generated_attempts": len(d),
        "reused_CID_in_attempt_dump": int(d.CID.duplicated().sum()),
        "accepted_detector_host_and_minimumcuts": int(d.FLAG_ACCEPT.sum()),
        "fitres": len(f),
        "quality": int(quality.sum()),
        "nonpositive_covariance_all_fitres": int((eig <= 0).sum()),
        "quality_pIa_gt05": int((quality & f.pIa.gt(0.5)).sum()),
        "quality_pIa_gt0999": int((quality & f.pIa.gt(0.999)).sum()),
        "negative_optical_extinction_all_fitres": int((~f.physical_screen).sum()),
        "negative_optical_extinction_quality_pIa_gt05": int(
            (quality & f.pIa.gt(0.5) & ~f.physical_screen).sum()
        ),
        "files_sha256": {
            str(p): sha(p)
            for p in [
                path,
                path.parent / "classifier.csv",
                WORK / "simulations" / name / (name + ".DUMP"),
            ]
        },
    }
    return f, counts


def matrix(f):
    return np.column_stack([np.ones(len(f)), f.x1, f.c, f.mass_step])


def coordinates(f):
    return np.column_stack([f.zHD / 0.1, f.x1, f.c / 0.1, f.HOST_LOGMASS])


def fit_predict(train, test):
    x = matrix(train)
    w = 1 / train.var_mu.to_numpy()
    y = train.y.to_numpy()
    coef = np.linalg.solve(x.T @ (w[:, None] * x), x.T @ (w * y))
    residual = y - x @ coef
    distance, idx = cKDTree(coordinates(train)).query(coordinates(test), k=40)
    nearw = w[idx]
    correction = (nearw * residual[idx]).sum(axis=1) / nearw.sum(axis=1)
    valid = distance[:, -1] <= 3
    # Predetermined redshift occupancy gate.
    nt = np.bincount(train.bin, minlength=5)
    ne = np.bincount(test.bin, minlength=5)
    valid &= (nt[test.bin] >= 20) & (ne[test.bin] >= 10)
    return (
        test.y.to_numpy() - matrix(test) @ coef - correction,
        valid,
        coef,
        distance[:, -1],
    )


def summarize(f, res, valid, weights=None):
    w = 1 / f.var_mu.to_numpy() if weights is None else weights / f.var_mu.to_numpy()
    means = []
    for b in range(5):
        use = valid & (f.bin.to_numpy() == b) & (w > 0)
        means.append(np.average(res[use], weights=w[use]) if use.any() else np.nan)
    # Conditional linear age association with free redshift-bin intercepts.
    use = valid & (w > 0)
    x = np.column_stack([np.eye(5)[f.bin.to_numpy()], f.age.to_numpy()])[use]
    coef = np.linalg.lstsq(
        x * np.sqrt(w[use, None]), res[use] * np.sqrt(w[use]), rcond=None
    )[0]
    return np.r_[means, coef[-1]]


def json_safe(x):
    if isinstance(x, np.ndarray):
        return json_safe(x.tolist())
    if isinstance(x, (list, tuple)):
        return [json_safe(v) for v in x]
    if isinstance(x, dict):
        return {k: json_safe(v) for k, v in x.items()}
    if isinstance(x, (float, np.floating)):
        return float(x) if np.isfinite(x) else None
    if isinstance(x, np.integer):
        return int(x)
    return x


def selected(f, variant):
    keep = f.quality.copy()
    if variant != "no_classifier":
        keep &= f.pIa.gt(0.999 if variant == "strict_classifier" else 0.5)
    if variant == "positive_extinction":
        keep &= f.physical_screen
    return f[keep].copy()


def covariance_summary(a):
    # Report covariance on jointly supported bootstrap draws, and expose conditioning.
    active = np.isfinite(a[:, :5]).any(axis=0)
    complete = np.isfinite(a[:, :5][:, active]).all(axis=1)
    c = np.full((5, 5), np.nan)
    if complete.sum() > 1:
        c[np.ix_(active, active)] = np.atleast_2d(
            np.cov(a[complete, :5][:, active], rowvar=False)
        )
    return {
        "zbin_covariance": c,
        "zbin_supported_bootstrap_replicates": np.isfinite(a[:, :5]).sum(axis=0),
        "covariance_complete_support_replicates": int(complete.sum()),
        "covariance_support_note": "Covariance conditions on all active redshift-bin occupancy and neighbor-support gates passing. Unsupported bins are null; this is not independent coverage calibration.",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--positive-dust", action="store_true")
    args = parser.parse_args()
    names = (
        ["SPP_TRAIN_NOMINAL", "SPP_EVAL_NOMINAL", "SPP_EVAL_AGE", "SPP_TRAIN_AGE"]
        if args.positive_dust
        else NAMES
    )
    frames = {}
    counts = {}
    for name in names:
        frames[name], counts[name] = load(name)
    design = json.loads((HERE / "analysis-design.json").read_text())
    output = {
        "interpretation": "Simulation-conditional response, not an empirical correction. Released W22 ages are mocked. Fixed cadence, reconstructed classifier, no core-collapse contaminants and no SALT retraining.",
        "code_sha256": sha(__file__),
        "analysis_design_sha256": sha(HERE / "analysis-design.json"),
        "counts": counts,
        "redshift_edges": EDGES.tolist(),
        "variance": "J C(x0,x1,c) J^T with fixed alpha=.15,beta=3.14 plus(.1mag)^2; used as a working weight, not a fitted full survey covariance.",
        "diagnostic_age_slope": "A secondary conditional association after five redshift-bin intercepts, declared in code before contrasts; simulated true age is used for this diagnostic only, never as a predictor.",
        "variants": {},
    }
    for variant in (
        ["primary"]
        if args.positive_dust
        else ["primary", "no_classifier", "strict_classifier", "positive_extinction"]
    ):
        fs = {n: selected(f, variant) for n, f in frames.items()}
        specs = [
            ("nominal_closure", "SPH_TRAIN_NOMINAL", "SPH_EVAL_NOMINAL"),
            ("age_nominal_correction", "SPH_TRAIN_NOMINAL", "SPH_EVAL_AGE"),
            ("dust_nominal_correction", "SPH_TRAIN_NOMINAL", "SPH_EVAL_DUST"),
            ("age_matched_correction", "SPH_TRAIN_AGE", "SPH_EVAL_AGE"),
        ]
        if args.positive_dust:
            specs = [
                (a, t.replace("SPH_", "SPP_"), e.replace("SPH_", "SPP_"))
                for a, t, e in [specs[0], specs[1], specs[3]]
            ]
        result = {}
        point_predictions = {}
        for label, t, e in specs:
            res, valid, coef, radius = fit_predict(fs[t], fs[e])
            point_predictions[label] = (res, valid)
            result[label] = {
                "train_n": len(fs[t]),
                "eval_n": len(fs[e]),
                "supported_n": int(valid.sum()),
                "supported_per_zbin": np.bincount(fs[e].bin[valid], minlength=5),
                "training_coefficients_intercept_x1_c_highmass": coef,
                "statistics_zbin_means_and_age_slope": summarize(fs[e], res, valid),
                "median_40th_neighbor_distance": np.median(radius),
                "max_40th_neighbor_distance": radius.max(),
            }
            fs[e].assign(residual=res, supported=valid).to_csv(
                WORK
                / (
                    ("positive-dust-" if args.positive_dust else "")
                    + variant
                    + "-"
                    + label
                    + ".csv"
                )
            )
        # Conditional fixed-cadence block bootstrap; each pool resampled independently.
        # Blocks shared across arms preserve their simulation coupling when present.
        blocks = {
            pool: np.unique(
                np.concatenate([f.block for name, f in fs.items() if pool in name])
            )
            for pool in ["TRAIN", "EVAL"]
        }
        rng = np.random.default_rng(2026092702)
        draws = []
        common_support_draws = []
        repetitions = 500 if variant == "primary" else 200
        for _ in range(repetitions):
            multiplicity = {}
            for pool, keys in blocks.items():
                sampled = rng.choice(keys, len(keys), replace=True)
                ids, nums = np.unique(sampled, return_counts=True)
                multiplicity[pool] = dict(zip(ids, nums))
            boots = {}
            for name, f in fs.items():
                m = (
                    f.block.map(multiplicity["TRAIN" if "TRAIN" in name else "EVAL"])
                    .fillna(0)
                    .to_numpy()
                    .astype(int)
                )
                boots[name] = f.iloc[np.repeat(np.arange(len(f)), m)]
            values = []
            boot_predictions = {}
            for label, t, e in specs:
                r, v, _, _ = fit_predict(boots[t], boots[e])
                values.append(summarize(boots[e], r, v))
                boot_predictions[label] = (r, v)
            if "age_matched_correction" in point_predictions:
                rn, vn = boot_predictions["age_nominal_correction"]
                ra, va = boot_predictions["age_matched_correction"]
                ff = boots["SPP_EVAL_AGE" if args.positive_dust else "SPH_EVAL_AGE"]
                common_support_draws.append(
                    summarize(ff, ra, va & vn) - summarize(ff, rn, va & vn)
                )
            draws.append(values)
        draws = np.array(draws)
        for i, (label, _, _) in enumerate(specs):
            a = draws[:, i, :]
            result[label]["bootstrap_standard_error"] = np.nanstd(a, axis=0, ddof=1)
            result[label]["bootstrap_percentile95"] = np.nanpercentile(
                a, [2.5, 97.5], axis=0
            ).T
            result[label].update(covariance_summary(a))
        differences = {}
        for label, i, j in (
            [
                ("age_minus_nominal", 1, 0),
                ("matched_age_minus_nominal_closure", 2, 0),
                ("matched_minus_nominal_correction_on_age", 2, 1),
            ]
            if args.positive_dust
            else [
                ("age_minus_nominal", 1, 0),
                ("dust_minus_nominal", 2, 0),
                ("matched_minus_nominal_correction_on_age", 3, 1),
                ("matched_age_minus_nominal_closure", 3, 0),
            ]
        ):
            a = draws[:, i, :] - draws[:, j, :]
            obs = np.array(
                result[specs[i][0]]["statistics_zbin_means_and_age_slope"]
            ) - np.array(result[specs[j][0]]["statistics_zbin_means_and_age_slope"])
            differences[label] = {
                "estimate": obs,
                "bootstrap_standard_error": np.nanstd(a, axis=0, ddof=1),
                "bootstrap_percentile95": np.nanpercentile(a, [2.5, 97.5], axis=0).T,
                **covariance_summary(a),
                "complete_bootstrap_replicates": int(np.isfinite(a).all(axis=1).sum()),
            }
        if "age_matched_correction" in point_predictions:
            rn, vn = point_predictions["age_nominal_correction"]
            ra, va = point_predictions["age_matched_correction"]
            ff = fs["SPP_EVAL_AGE" if args.positive_dust else "SPH_EVAL_AGE"]
            a = np.array(common_support_draws)
            differences[
                "matched_minus_nominal_correction_common_evaluation_support"
            ] = {
                "n": int((vn & va).sum()),
                "estimate": summarize(ff, ra, va & vn) - summarize(ff, rn, va & vn),
                "bootstrap_standard_error": np.nanstd(a, axis=0, ddof=1),
                "bootstrap_percentile95": np.nanpercentile(a, [2.5, 97.5], axis=0).T,
                **covariance_summary(a),
            }
        output["variants"][variant] = {
            "arms": result,
            "differences": differences,
            "bootstrap_replicates": repetitions,
        }
        np.savez_compressed(
            WORK
            / (
                ("positive-dust-" if args.positive_dust else "")
                + variant
                + "-bootstrap.npz"
            ),
            draws=draws,
        )
    (
        RESULTS
        / (
            "positive-dust-response.json"
            if args.positive_dust
            else "correction-response.json"
        )
    ).write_text(json.dumps(json_safe(output), indent=2) + "\n")
    print(json.dumps(json_safe(output["variants"]["primary"]), indent=2))


if __name__ == "__main__":
    main()
