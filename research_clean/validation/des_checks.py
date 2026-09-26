#!/usr/bin/env python3
"""Independent algebraic checks of frozen DES replay outputs.

Run from any directory with the clean environment. This does not fit cosmology,
reconstruct selection, or replace the frozen native flux covariance.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
os.environ.setdefault("JAX_PLATFORMS", "cpu")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import multivariate_normal, norm, rankdata, t as student_t


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def split_rhat(x):
    """Classical split R-hat, independently recomputed; not rank-normalized."""
    x = np.asarray(x)
    n = x.shape[1] // 2
    x = np.concatenate([x[:, :n], x[:, -n:]], axis=0)
    within = np.mean(np.var(x, axis=1, ddof=1), axis=0)
    between = n * np.var(np.mean(x, axis=1), axis=0, ddof=1)
    return np.sqrt(((n - 1) / n * within + between / n) / within)


def rank_rhat(x):
    """Rank-normalized and median-folded split diagnostics (maximum)."""
    values = np.asarray(x).reshape(x.shape[0], x.shape[1], -1)
    output = []
    for i in range(values.shape[-1]):
        a = values[..., i]
        diagnostic = []
        for transformed in (a, np.abs(a-np.median(a))):
            ranks = rankdata(transformed.ravel(), method="average")
            normal = norm.ppf((ranks-.375)/(len(ranks)+.25)).reshape(a.shape)
            diagnostic.append(float(split_rhat(normal)))
        output.append(max(diagnostic))
    return np.asarray(output)


def numpy_predict(params, data, name, noise="gaussian"):
    """NumPy implementation independent of the JAX feature/projection code."""
    names = {
        "none": [], "stretch": [0], "colour": [1], "tripp": [0, 1],
        "host": [0, 1, 2], "host_colour": [0, 1, 2, 3],
        "broken_colour": [0, 1, 2, 4], "evolution": [0, 1, 2, 5, 6],
        "flexible": list(range(8)),
    }[name]
    x, c = data["y"][:, 1], data["y"][:, 2]
    g = data["z"] / (1 + data["z"]) - 0.2
    baseline = data["distance_reference"] + params["M"] + np.interp(
        data["z"], [0.01, .10, .20, .35, .50, .70, .90, 1.20],
        np.r_[0., params["offsets"]],
    )
    co = np.zeros(8)
    co[names] = params.get("coefficients", [])
    means, variances = [], []
    for h in (0, 1):
        hc = h - .5
        means.append(baseline + co[0]*x + co[1]*c + co[2]*hc
                     + co[3]*hc*c + co[4]*np.maximum(c, 0)
                     + co[5]*x*g + co[6]*c*g + co[7]*x*x)
        dx = co[0] + co[5]*g + 2*co[7]*x
        dc = co[1] + co[3]*hc + co[4]*(c > 0) + co[6]*g
        cov = data["cov"]
        variance = (cov[:, 0, 0] + dx*dx*cov[:, 1, 1] + dc*dc*cov[:, 2, 2]
                    - 2*dx*cov[:, 0, 1] - 2*dc*cov[:, 0, 2]
                    + 2*dx*dc*cov[:, 1, 2] + params["scatter"][h]**2)
        variances.append(variance)
    means, variances = np.array(means).T, np.array(variances).T
    y = data["y"][:, :1]
    lp = (norm.logpdf(y, means, np.sqrt(variances)) if noise == "gaussian"
          else student_t.logpdf(y, 4, means, np.sqrt(variances/2)))
    p = np.clip(data["host_prob"], 1e-10, 1-1e-10)
    logp = np.logaddexp(np.log1p(-p)+lp[:, 0], np.log(p)+lp[:, 1])
    mean = (1-p)*means[:, 0] + p*means[:, 1]
    variance = ((1-p)*variances[:, 0] + p*variances[:, 1]
                + p*(1-p)*(means[:, 1]-means[:, 0])**2)
    return logp, mean, variance


def predictor_checks(base, colour, second_seed):
    import jax
    import jax.numpy as jnp
    from lib.predictors import predictive, MODELS

    with np.load(ROOT / "data/des/predictors/data.npz", allow_pickle=False) as f:
        arrays = dict(f)
    rows = pd.read_csv(ROOT / "data/des/predictors/rows.csv", dtype={"CID": str})
    keep = (arrays["survey"] == 10) & (arrays["pIa"] > .999)
    test = keep & (arrays["fold"] == 0)
    cohort = read(ROOT / "data/des/predictors/cohort.json")
    assert len(rows) == len(arrays["y"])
    assert np.array_equal(rows.fold, arrays["fold"])
    assert np.array_equal(rows.IDSURVEY, arrays["survey"])
    assert np.isfinite(arrays["y"][keep]).all()
    assert np.isfinite(arrays["cov"][keep]).all()
    assert rows.loc[keep, "CID"].is_unique
    assert dict(zip(rows.loc[keep, "CID"], map(int, arrays["fold"][keep]))) == cohort["cid_to_fold"]
    cov = arrays["cov"][keep]
    assert np.max(abs(cov - cov.transpose(0, 2, 1))) < 1e-12
    eig = np.linalg.eigvalsh(cov)
    assert eig.min() > 0
    data = {key: value[test] for key, value in arrays.items()}
    # Exact coasting luminosity distance, with heliocentric redshift factor.
    data["distance_reference"] = (5*np.log10((1+data["zhel"])*np.log1p(data["z"])
                                          * 299792.458/70) + 25)
    reports = []
    all_ll = []
    for directory, name in [(base, "tripp"), (colour, "colour"), (second_seed, "tripp")]:
        with np.load(directory / "chains.npz", allow_pickle=False) as f:
            chains = dict(f)
        flat = {key: value.reshape((-1, *value.shape[2:])) for key, value in chains.items()}
        draws = len(flat["M"])
        lp, means, variances = [], [], []
        for i in range(draws):
            a, b, c = numpy_predict({k: v[i] for k, v in flat.items()}, data, name)
            lp.append(a); means.append(b); variances.append(c)
        lp, means, variances = map(np.asarray, (lp, means, variances))
        score = logsumexp(lp, axis=0)-np.log(draws)
        mean = means.mean(0)
        sd = np.sqrt(variances.mean(0)+means.var(0))
        table = pd.read_csv(directory / "heldout.csv", dtype={"CID": str})
        assert table.CID.tolist() == rows.loc[test, "CID"].tolist()
        assert np.array_equal(table.observed_mB.to_numpy(), data["y"][:, 0]) or np.allclose(table.observed_mB, data["y"][:, 0], atol=1e-13, rtol=0)
        errors = {
            "log_density": float(np.max(abs(score-table.log_predictive_density))),
            "mean": float(np.max(abs(mean-table.predictive_mean))),
            "sd": float(np.max(abs(sd-table.predictive_sd))),
        }
        assert max(errors.values()) < 1e-10, errors
        rhats = np.concatenate([np.atleast_1d(split_rhat(v)).ravel() for v in chains.values()])
        assert rhats.max() < 1.01
        robust_rhats = np.concatenate([rank_rhat(v) for v in chains.values()])
        assert robust_rhats.max() < 1.01
        lag_one = []
        for value in chains.values():
            centered = value.reshape(value.shape[0], value.shape[1], -1).copy()
            centered -= centered.mean(axis=1, keepdims=True)
            lag_one.extend(((centered[:, :-1]*centered[:, 1:]).sum(1)/(centered*centered).sum(1)).ravel())
        # Independent within-chain batch estimate of integration Monte Carlo error.
        batch_scores = []
        for chain in lp.reshape(4, -1, lp.shape[-1]):
            for block in np.array_split(chain, 5):
                batch_scores.append((logsumexp(block, axis=0)-np.log(len(block))).sum())
        reports.append({"model": name, "result_directory": directory.name, "draws": draws,
            "score_sum": float(score.sum()), "rmse_mag": float(np.sqrt(np.mean((data["y"][:, 0]-mean)**2))),
            "max_absolute_numpy_reconstruction_errors": errors,
            "max_independent_split_rhat": float(rhats.max()),
            "max_rank_normalized_folded_split_rhat": float(robust_rhats.max()),
            "within_chain_lag_one_correlation_range": [float(min(lag_one)), float(max(lag_one))],
            "approximate_total_log_score_batch_mcse": float(np.std(batch_scores, ddof=1)/np.sqrt(len(batch_scores))),
            "coefficients_mean": flat["coefficients"].mean(0).tolist(),
            "coefficients_sd": flat["coefficients"].std(0, ddof=1).tolist(),
            "saved_convergence": read(directory / "convergence.json"),
        })
        all_ll.append(score)
    # Explicit paired physical-object comparison; field blocks preserve field dependence.
    delta = all_ll[0]-all_ll[1]
    rng = np.random.default_rng(20260926)
    obj = delta[rng.integers(0, len(delta), size=(10000, len(delta)))].mean(1)
    fields = rows.loc[test, "field"].to_numpy()
    labels = np.unique(fields)
    sums = np.array([delta[fields == field].sum() for field in labels])
    counts = np.array([(fields == field).sum() for field in labels])
    ix = rng.integers(0, len(labels), size=(10000, len(labels)))
    block = sums[ix].sum(1)/counts[ix].sum(1)
    # Compiled NumPy/JAX checks cover both residual laws and every supplied mean family.
    numerical = []
    for name, features in MODELS.items():
        p = {"M": -19.2, "offsets": np.linspace(-.1, .1, 7),
             "scatter": np.array([.11, .18]), "coefficients": np.linspace(-.15, .45, len(features))}
        d = {k: jnp.asarray(v[:17]) for k, v in data.items()}
        dn = {k: np.asarray(v) for k, v in d.items()}
        for noise in ("gaussian", "student4"):
            expected = numpy_predict(p, dn, name, noise)
            actual = jax.jit(lambda q: predictive(q, d, name, noise))(p)
            error = max(float(np.max(abs(a-b))) for a, b in zip(expected, actual))
            assert error < 1e-9, (name, noise, error)
            if features:
                f = jax.jit(jax.grad(lambda co: predictive({**p, "coefficients": co}, d, name, noise)[0].sum()))
                ad = np.asarray(f(p["coefficients"]))
                fd = []
                for i in range(len(features)):
                    step = np.zeros(len(features)); step[i] = 1e-5
                    fd.append((numpy_predict({**p, "coefficients": p["coefficients"]+step}, dn, name, noise)[0].sum()
                               -numpy_predict({**p, "coefficients": p["coefficients"]-step}, dn, name, noise)[0].sum())/2e-5)
                gradient_error = float(np.max(abs(ad-fd)))
                assert gradient_error < 2e-5, (name, noise, gradient_error)
            else:
                gradient_error = 0.
            numerical.append({"model": name, "noise": noise, "max_value_error": error, "max_gradient_error": gradient_error})
    seed_delta = float(all_ll[2].sum()-all_ll[0].sum())
    seed_mcse = float(np.hypot(reports[0]["approximate_total_log_score_batch_mcse"], reports[2]["approximate_total_log_score_batch_mcse"]))
    return {"cohort": {"total": int(keep.sum()), "training": int((keep & ~test).sum()), "test": int(test.sum()),
                       "fields": len(labels), "minimum_covariance_eigenvalue": float(eig.min())},
            "reconstructed_runs": reports, "compiled_likelihood_checks": numerical,
            "second_seed_check": {"score_difference_nats": seed_delta, "combined_approximate_batch_mcse": seed_mcse,
                                   "difference_over_mcse": seed_delta/seed_mcse},
            "tripp_minus_colour": {"sum_nats": float(delta.sum()), "mean_nats": float(delta.mean()),
                                   "paired_object_95_percent_mean_interval": np.quantile(obj, [.025, .975]).tolist(),
                                   "field_block_95_percent_mean_interval": np.quantile(block, [.025, .975]).tolist(),
                                   "leave_one_field_out_mean_range": [float(min((delta[fields != x].mean() for x in labels))), float(max((delta[fields != x].mean() for x in labels)))],
                                   "interval_scope": "Object and field resampling of fixed predictions, not calibration or model-selection uncertainty."}}


def calibration_checks(directory):
    arrays, reports = [], []
    for name in ("discovery", "validation"):
        with np.load(ROOT / f"data/calibration/{name}-projected-modes.npz", allow_pickle=False) as f:
            q = dict(f)
        a = np.column_stack([q["calibration_design"], q["observer_design"]])
        f, u = a.T@a, a.T@q["residual"]
        assert np.allclose(f, q["joint_F"].sum(0), atol=2e-9, rtol=1e-10)
        assert np.allclose(u, q["joint_u"].sum(0), atol=1e-9, rtol=1e-10)
        for i, (start, stop) in enumerate(zip(q["row_offsets"][:-1], q["row_offsets"][1:])):
            assert np.allclose(a[start:stop].T@a[start:stop], q["joint_F"][i], atol=1e-9, rtol=1e-9)
            assert np.allclose(a[start:stop].T@q["residual"][start:stop], q["joint_u"][i], atol=1e-9, rtol=1e-9)
        reports.append({"cohort": name, "objects": len(q["CID"]), "projected_rows": len(a),
                        "gram_absolute_error": float(abs(f-q["joint_F"].sum(0)).max()),
                        "residual_projection_absolute_error": float(abs(u-q["joint_u"].sum(0)).max())})
        arrays.append((a, q["residual"], q["CID"].astype(str)))
    assert not set(arrays[0][2]) & set(arrays[1][2])
    with np.load(ROOT / "data/calibration/distance-responses.npz", allow_pickle=False) as f:
        response = dict(f)
    rows = pd.read_csv(ROOT / "data/des/predictors/rows.csv", dtype={"CID": str})
    with np.load(ROOT / "data/des/predictors/data.npz") as f:
        zmap = dict(zip(rows.CID, f["zhel"]))
    ids = response["validation_CID"].astype(str)
    order = sorted(range(len(ids)), key=lambda i: (zmap[ids[i]], ids[i]))
    contrast = response["validation"][order[-255:]].mean(0)-response["validation"][order[:255]].mean(0)
    assert np.allclose(contrast, response["contrast"], atol=1e-14, rtol=0)
    h = contrast.copy(); h[12:] /= .02
    bandmap = np.array([[1., 0, 0], [-1, -1, -1], [0, 1, 0], [0, 0, 1]])
    isotropic = np.eye(15)
    isotropic[12:, 12:] = np.linalg.cholesky(2*.02**2*np.linalg.inv(bandmap.T@bandmap))
    transforms = {"systematics_only": np.eye(15)[:, :12], "inherited_observer": np.diag(np.r_[np.ones(12), [.02]*3]), "isotropic_observer": isotropic}
    expected = {r["family"]: r for r in read(directory / "summary.json")["records"]}
    outputs = []
    for name, transform in transforms.items():
        ad, av = (entry[0]@transform for entry in arrays)
        yd, yv = (entry[1] for entry in arrays)
        target = h@transform
        covariances, means = [], []
        for design, residual in ((ad, yd), (np.vstack([ad, av]), np.r_[yd, yv])):
            augmented = np.vstack([design, np.eye(len(target))])
            mean = np.linalg.lstsq(augmented, np.r_[residual, np.zeros(len(target))], rcond=None)[0]
            _, r = np.linalg.qr(augmented, mode="reduced")
            rinv = np.linalg.solve(r, np.eye(len(r)))
            covariance = rinv@rinv.T
            means.append(mean); covariances.append(covariance)
        sd = [np.sqrt(target@c@target) for c in covariances]
        assert np.allclose(sd, [expected[name]["discovery_sd_mag"], expected[name]["combined_sd_mag"]], atol=1e-12)
        assert np.allclose([target@m for m in means], [expected[name]["discovery_mean_response_mag"], expected[name]["combined_mean_response_mag"]], atol=1e-12)
        b = av@np.linalg.cholesky(covariances[0])
        residual = yv-av@means[0]
        system = np.eye(len(target))+b.T@b
        v = b.T@residual
        evidence = .5*(yv@yv-residual@residual+v@np.linalg.solve(system, v)-np.linalg.slogdet(system)[1])
        assert abs(evidence-expected[name]["validation_log_evidence_increment"]) < 1e-9
        outputs.append({"family": name, "discovery_sd_mag": float(sd[0]), "combined_sd_mag": float(sd[1]),
                        "combined_mean_response_mag": float(target@means[1]),
                        "independent_predictive_evidence_increment": float(evidence),
                        "remaining_prior_variance_fraction": float(sd[1]**2/(target@target))})
    return {"sufficient_statistic_checks": reports, "contrast_reconstruction_max_error": float(abs(contrast-response["contrast"]).max()), "independent_qr_and_predictive_results": outputs,
            "scope": "Exact conditional Gaussian calculation for supplied projected designs and priors; upstream pixel/projection construction and unknown modes are not validated here."}


def flux_checks(directory):
    from lib.flux_engine import Engine
    from lib.flux_fit import Case
    engine = Engine()
    full, held = read(directory / "fits.json"), read(directory / "heldout.json")
    cases = {cid: Case(engine, cid) for cid in sorted({r["CID"] for r in full})}
    objective_error, native_error, stationarity, integration, forward = [], [], [], [], []
    for cid, case in cases.items():
        assert np.linalg.eigvalsh(case.cov).min() > 0
        inv = case.q["inverse_frozen_flux_covariance"]
        assert np.allclose(inv@case.cov, np.eye(len(inv)), atol=1e-9, rtol=0)
        residual = case.y-case.q["model_flux"]
        native = residual@np.linalg.solve(case.cov, residual)+float(case.q["prior_chi2"].ravel()[0])
        native_error.append(abs(native-float(case.q["chi2"].ravel()[0])))
        p = case.p0.copy()
        fine = engine.prepare(case.z, case.ebv, 1.)
        coarse = engine.flux(p, case.b, case.t, case.z, case.grid, interpolation="sncosmo")
        refined = engine.flux(p, case.b, case.t, case.z, fine, interpolation="sncosmo")
        public = engine.sncosmo_flux(p, case.b, case.t, case.z, case.ebv)
        scale = max(abs(coarse).max(), 1e-20)
        integration.append(float(abs(coarse-refined).max()/scale))
        forward.append(float(abs(coarse-public).max()/scale))
    for fit in full:
        case = cases[fit["CID"]]
        x = np.array(fit["x"])
        f = lambda p: case.flux(p, fit["family"], fit["interpolation"])
        residual = case.y-f(x)
        chi2 = residual@np.linalg.solve(case.cov, residual)+((x[3]-case.priorcentre)/10)**2
        objective_error.append(abs(chi2-fit["chi2_total"]))
        step = np.array([1e-4, 1e-3, 1e-4, 3e-3])
        jac = np.column_stack([(f(x-2*np.eye(4)[i]*step[i])-8*f(x-np.eye(4)[i]*step[i])+8*f(x+np.eye(4)[i]*step[i])-f(x+2*np.eye(4)[i]*step[i]))/(12*step[i]) for i in range(4)])
        grad = -jac.T@np.linalg.solve(case.cov, residual)
        grad[3] += (x[3]-case.priorcentre)/100
        hc = np.array(fit["hessian_covariance"])
        stationarity.append(float(np.sqrt(max(0, grad@hc@grad))))
    predictive_errors = []
    for fit in held:
        case = cases[fit["CID"]]
        train, test = np.array(fit["train_rows"]), np.array(fit["test_rows"])
        assert len(set(train) & set(test)) == 0 and sorted(np.r_[train, test]) == list(range(len(case.y)))
        for i in test:
            assert int(hashlib.sha256(f"independent-flux-v1|{case.cid}|{int(np.floor(case.t[i]+.5))}".encode()).hexdigest(), 16) % 5 == 0
        # Saved conditional mean/diagonal are checked against a full independent covariance.
        x = np.array(fit["x"])
        f = lambda p: case.flux(p, fit["family"], "sncosmo")
        step = np.array([2e-5, 2e-4, 2e-5, 2e-3]+([2e-5] if len(x)==5 else []))
        jac = np.column_stack([(f(x+np.eye(len(x))[i]*step[i])-f(x-np.eye(len(x))[i]*step[i]))/(2*step[i]) for i in range(len(x))])
        cross = case.cov[np.ix_(test, train)]
        a = cross@np.linalg.inv(case.cov[np.ix_(train, train)])
        mean = f(x)[test]+a@(case.y[train]-f(x)[train])
        noise = case.cov[np.ix_(test, test)]-a@cross.T
        j = jac[test]-a@jac[train]
        covariance = noise+j@np.array(fit["hessian_covariance"])@j.T
        score = multivariate_normal.logpdf(case.y[test], mean=mean, cov=covariance)
        predictive_errors.append(abs(score-fit["prediction"]["logpdf_laplace"]))
        assert np.allclose(mean, fit["prediction"]["test_model"], atol=1e-10)
        assert np.allclose(np.diag(covariance), fit["prediction"]["test_variance"], atol=1e-10)
    assert max(objective_error) < 1e-9
    assert max(native_error) < 1e-9
    assert max(predictive_errors) < 1e-8
    assert max(stationarity) < 1e-3
    assert max(integration) < 1e-4
    return {"full_fits": len(full), "heldout_fits": len(held),
            "native_chi2_reconstruction_max_error": max(native_error),
            "fitted_objective_reconstruction_max_error": max(objective_error),
            "heldout_full_covariance_logpdf_max_error": max(predictive_errors),
            "max_local_covariance_scaled_gradient_norm": max(stationarity),
            "max_5A_to_1A_flux_difference_over_object_peak": max(integration),
            "max_public_sncosmo_bandflux_difference_over_object_peak": max(forward),
            "max_two_start_chi2_gap": max(f["two_start_chi2_gap"] for f in full+held),
            "max_half_step_hessian_relative_difference": max(f["hessian_stepsize_relative_difference"] for f in full+held),
            "native_reference_max_abs_delta_mB_x1_c_t0": read(directory / "summary.json")["max_abs_delta_mB_x1_c_t0"],
            "scope": "Checks frozen accepted epochs, mean integration, fixed covariance, prior and local fit/prediction algebra; not historical detector reduction or selection."}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", default="revalidation-20260926-v1")
    parser.add_argument("--colour", default="revalidation-20260926-des-colour")
    parser.add_argument("--second-seed", default="revalidation-20260926-des-second-seed")
    args = parser.parse_args()
    result = {"schema": "des-independent-checks-v1", "source_sha256": sha(Path(__file__)),
              "status": "passed", "replay_prefix": args.prefix,
              "calibration": calibration_checks(ROOT / "results" / (args.prefix+"-calibration")),
              "flux": flux_checks(ROOT / "results" / (args.prefix+"-des-flux")),
              "predictors": predictor_checks(ROOT / "results" / (args.prefix+"-des-predictors"), ROOT / "results" / args.colour, ROOT / "results" / args.second_seed)}
    paths = list((ROOT / "data/des").rglob("*"))+list((ROOT / "data/calibration").rglob("*"))
    for directory in (args.prefix+"-calibration", args.prefix+"-des-flux", args.prefix+"-des-predictors", args.colour, args.second_seed):
        paths.extend((ROOT / "results" / directory).glob("*"))
    paths.extend(ROOT / x for x in ["lib/predictors.py", "lib/flux_engine.py", "lib/flux_fit.py", "workflows/calibration.py", "workflows/des_predictors.py", "workflows/des_flux.py"])
    result["inputs_sha256"] = {str(p.relative_to(ROOT)): sha(p) for p in sorted(paths) if p.is_file()}
    result["literature_scope"] = {"frozen_release": "DES-SN5YR assets, not DES-Dovekie",
        "original_cosmology_paper": "https://arxiv.org/abs/2401.02929",
        "later_reanalysis": "https://arxiv.org/abs/2511.07517v3",
        "checked_utc_date": "2026-09-26",
        "comparison_limit": "This audit independently checks fixed flux objectives and selected-sample predictors, not the full published cosmology or later recalibrated release."}
    out = ROOT / "validation/reports/des.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"status": result["status"], "output": str(out), "flux": result["flux"], "paired_predictor_comparison": result["predictors"]["tripp_minus_colour"]}, indent=2))


if __name__ == "__main__":
    main()
