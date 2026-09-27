#!/usr/bin/env python3
"""Validate exact input identities, independent pools, covariance and native outputs."""
import datetime
import json
import sys
import numpy as np
from common import *
from analyze import dump, load, json_safe


def main():
    campaigns = [
        json.loads((RESULTS / n).read_text())
        for n in ["campaign.json", "positive-dust-campaign.json"]
    ]
    jobs = [j for c in campaigns for j in c["jobs"]]
    runs = []
    classifiers = []
    independence = {}
    injections = {}
    covariances = {}
    for job in jobs:
        folder = WORK / "fits" / job["name"]
        r = json.loads((folder / "run.json").read_text())
        assert (
            r["fit_graceful"] and r["generation_returncode"] == r["fit_returncode"] == 0
        )
        assert sha(job["input"]) == job["input_sha256"] == r["hashes"][job["input"]]
        for name in ["fit.FITRES.TEXT", "fit.SNANA.TEXT"]:
            assert sha(folder / name) == r[name]["sha256"]
        cl = json.loads((folder / "classifier.json").read_text())
        assert cl["native_fit_sha256"] == sha(folder / "fit.FITRES.TEXT")
        assert (
            cl["single_batch_max_probability_difference"] < 1e-6
            and cl["truth_inputs_used"] is False
        )
        runs.append(r)
        classifiers.append(cl)
        f, _ = load(job["name"])
        if job["arm"] == "age":
            error = np.max(abs(f.SIM_gammaDM + 0.03 * (f.age - 3)))
            assert error < 2e-7
            injections[job["name"]] = {
                "maximum_gammaDM_identity_error_mag": float(error),
                "grey_shift_formula": "-0.030*(mock_SN_age_Gyr-3)",
            }
    for prefix in ["SPH", "SPP"]:
        train = dump(prefix + "_TRAIN_NOMINAL")
        ev = dump(prefix + "_EVAL_NOMINAL")
        cols = ["LIBID", "GENZ", "GALID", "PEAKMJD", "SN_age", "SALT2x1", "SALT2c"]
        overlap = train.merge(ev, on=cols)
        assert len(overlap) == 0
        independence[prefix] = {
            "train_seed": 270911 if prefix == "SPH" else 270921,
            "evaluation_seed": 270912 if prefix == "SPH" else 270922,
            "identical_full_latent_tuples_across_pools": 0,
            "shared_cadence_LIBIDs": len(set(train.LIBID) & set(ev.LIBID)),
            "scope": "Independent generated realizations conditional on shared cadence, W22 mocked population and released surface/classifier; not independent surveys.",
        }
    for filename, bootprefix in [
        ("correction-response.json", ""),
        ("positive-dust-response.json", "positive-dust-"),
    ]:
        r = json.loads((RESULTS / filename).read_text())
        assert r["code_sha256"] == sha(HERE / "analyze.py")
        boot = np.load(WORK / (bootprefix + "primary-bootstrap.npz"))["draws"]
        slopes = boot[:, :, -1]
        ratios = np.std(slopes[:250], axis=0, ddof=1) / np.std(slopes, axis=0, ddof=1)
        eigen = {}
        for name, row in r["variants"]["primary"]["differences"].items():
            c = np.array(row["zbin_covariance"], dtype=float)
            active = np.isfinite(np.diag(c))
            if active.any():
                m = c[np.ix_(active, active)]
                assert np.isfinite(m).all()
                e = np.linalg.eigvalsh(m)
                assert e.min() > -1e-12
                eigen[name] = float(e.min())
        covariances[filename] = {
            "minimum_supported_covariance_eigenvalues": eigen,
            "first250_over_full500_age_slope_SE_by_arm": ratios,
            "coverage_status": "Bootstrap convergence/sensitivity only; frequentist coverage across new surveys has not been established.",
        }
    review = (
        ROOT / "studies/host_ages/results/galaxy_validation/survey-physics-review.json"
    )
    result = {
        "validated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "code_sha256": sha(__file__),
        "generated_attempts": sum(j["attempts"] for j in jobs),
        "all_native_runs_graceful": True,
        "native_run_count": len(runs),
        "native_generation_seconds": sum(r["generation_seconds"] for r in runs),
        "native_fit_seconds": sum(r["fit_seconds"] for r in runs),
        "classifier_inference_seconds": sum(
            c["inference_seconds"] for c in classifiers
        ),
        "classifier_all_single_batch_checks_pass": True,
        "independence": independence,
        "injected_sign": injections,
        "covariances": covariances,
        "independent_physics_review": {
            "path": str(review.relative_to(ROOT)),
            "sha256": sha(review),
        },
        "source_code_current_sha256": {
            str(p.relative_to(ROOT)): sha(p)
            for p in sorted(HERE.iterdir())
            if p.is_file()
        },
        "binary_sha256": {
            str(p.relative_to(ROOT)): sha(p)
            for p in [
                WORK / "SNANA/bin/snlc_sim.exe",
                WORK / "SNANA/bin/snlc_fit.exe",
                WORK / "SNANA/bin/SALT2mu.exe",
                WORK / "SNANA-legacy-eff/bin/snlc_sim.exe",
                WORK / "SNANA-legacy-eff/bin/SALT2mu.exe",
            ]
        },
        "model_and_cadence_sha256": {
            str(p.relative_to(WORK)): sha(p)
            for p in [
                DATA / "simlib/DES/DES-SN5YR_DES.SIMLIB",
                DATA / "simlib/DES/DES-SN5YR_DES_FLUXERRMODEL_SIM.DAT",
                DATA / "models/searcheff/SEARCHEFF_PIPELINE_DES.DAT",
                DATA / "models/searcheff/SEARCHEFF_PIPELINE_LOGIC.DAT",
                DATA / "kcor/Dovekie/calib_DES-SN5YR_DES.fits.gz",
                *sorted((DATA / "models/SALT3/SALT3.DOVEKIE").glob("*")),
                *[
                    DATA / "models/classifiers/DES-SN5YR/SNNTRAINV19_z_TRAINDES_V19" / n
                    for n in ["model.pt", "data_norm.json", "cli_args.json"]
                ],
            ]
            if p.is_file()
        },
        "limitations": [
            "Production 4D BBC support fails; native 1D and 40-neighbor estimates are separate benchmarks.",
            "Mock host mass has zero observational error.",
            "No SALT surface or classifier retraining; pure Ia signal acceptance does not calibrate contamination.",
            "Nominal age/redshift closure is not zero at the achieved precision; injection differences must be distinguished from absolute unbiasedness.",
            "No measurement of real progenitor ages or cosmological correction.",
        ],
    }
    (RESULTS / "native-execution.json").write_text(json.dumps(runs, indent=2) + "\n")
    (RESULTS / "classifier-execution.json").write_text(
        json.dumps(classifiers, indent=2) + "\n"
    )
    (RESULTS / "validation.json").write_text(
        json.dumps(json_safe(result), indent=2) + "\n"
    )
    print(
        json.dumps(
            {
                k: result[k]
                for k in [
                    "generated_attempts",
                    "native_run_count",
                    "native_generation_seconds",
                    "native_fit_seconds",
                    "classifier_inference_seconds",
                ]
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
