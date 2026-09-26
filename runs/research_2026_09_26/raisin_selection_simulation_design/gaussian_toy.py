#!/usr/bin/env python3
"""Fixed-truth selection checks. No observed flux is read or fitted.

Run `freeze` before `run`; protocol/source hashes must remain unchanged.
The only input table read at freeze is cohort metadata to name a native pilot.
"""
import argparse
import csv
import hashlib
import json
import platform
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.stats import binom, norm

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def freeze():
    paths = [
        "docs/research-2026-09-26/raisin-flux-sign-audit.md",
        "docs/research-2026-09-26/raisin-simulation-sign-accounting.md",
        "docs/research-2026-09-26/raisin-signed-refit-design.md",
        "runs/research_2026_09_26/raisin_sign_source/cadence-result.json",
        "runs/research_2026_09_26/raisin_sign_source/cadence-protocol.json",
        "runs/research_2026_09_26/raisin_sign_source/sim/inputs/DES/sim_DES_SNOOPY.input",
        "runs/research_2026_09_26/raisin_sign_source/sim/inputs/DES/SIMGEN_MASTER_DESSPEC.INPUT",
        "runs/research_2026_09_26/raisin_sign_source/fit/sim/REFAC_DES_RAISIN_opt.nml",
        "runs/research_2026_09_26/raisin_sign_source/raisin_cosmo/genSimlib.py",
        "runs/research_2026_09_26/raisin_flux_sign_audit/lineage-handoff.json",
        "runs/research_2026_09_26/astra_design/raisin_signed_refit/cohort.csv",
        "runs/research_2026_09_26/astra_design/raisin_signed_refit/source-row-ledger.csv",
        "runs/research_2026_09_26/astra_design/raisin_signed_refit/input-manifest.json",
    ]
    cohort = list(csv.DictReader((ROOT / paths[-3]).open()))
    cohort.sort(key=lambda r: (float(r["zHEL"]), r["CID"]))
    pilot = [{k: r[k] for k in ["CID", "zHEL", "peak_header", "raw_path"]}
             for r in [cohort[0], cohort[-1]]]
    protocol = {
        "status": "Frozen before synthetic noise outcomes; no observed fits or flux-as-truth",
        "source_sha256": sha(__file__),
        "inputs_sha256": {p: sha(ROOT / p) for p in paths},
        "seeds": {"scalar": 26092601, "archive": 26092602, "paired": 26092603},
        "scalar": {"snr": [0, 0.5, 1, 2, 5], "sigma": 1.0, "draws_per_snr": 300000},
        "paired": {"truth_flux": 0.5, "sigma": 1.0, "epochs": 64, "realizations": 50000},
        "arms": {
            "F": "All 64 epochs, new signed noise",
            "T": "All 64 epochs, retain newly generated y>0",
            "M": "One independent synthetic archived positive mask, new signed noise",
            "MT": "Same archived mask, additionally retain newly generated y>0",
        },
        "estimator": "Unconstrained arithmetic flux mean of retained values; record every empty arm. No magnitude conversion.",
        "shared_noise": "One complete new Gaussian vector per realization, subset without redrawing; archive uses independent seed.",
        "gates": {"maximum_mc_z": 6.0, "quadrature_absolute": 1e-10},
        "native_pilot_proposal_not_executed": {
            "metadata_rule": "minimum and maximum zHEL in frozen ten-DES16 cohort; CID resolves ties",
            "objects": pilot,
            "truths_stretch_AV": [[1.0, 0.0], [0.85, 0.3], [1.15, 0.3]],
            "distance_reference": "Flat LCDM H0=70 km/s/Mpc, OmegaM=.3 at fixed source redshift; a signal-generation convention, not a cosmology result",
            "timing": "Fixed common header peak first; physical rest phases [-7,45] only; native support/accepted-row gates before noise",
            "noise": "64 shared draws per truth and cadence, fixed raw-author quoted-error diagonal as conditional engineering noise; source error-model simulation is a later distinct gate",
            "fit_count": 1536,
            "model": "Released native SNooPy.B18, fixed host RV1.518 and existing signed-refit native configuration",
            "resource": "At most 2 CPU processes, one thread each; time first 32 fits; stop after pilot and report failures/Monte Carlo precision",
            "scope": "No intrinsic-scatter or survey-selection claim; never use observed fitted parameters or flux as truth",
        },
    }
    path = OUT / "protocol.json"
    assert not path.exists(), "Preserve existing frozen protocol"
    save(path, protocol)
    print(json.dumps({"protocol_sha256": sha(path), "native_pilot": pilot}))


def stats(values):
    valid = np.isfinite(values)
    x = values[valid]
    return {"n": len(x), "failed": int((~valid).sum()), "mean": float(x.mean()),
            "sd": float(x.std(ddof=1)), "mcse": float(x.std(ddof=1) / np.sqrt(len(x)))}


def run():
    start = time.monotonic()
    protocol = json.loads((OUT / "protocol.json").read_text())
    assert sha(__file__) == protocol["source_sha256"]
    # Source provenance was frozen before the toy. Some other agents may append
    # reports later; toy execution depends only on the immutable protocol/code.
    gate_z = protocol["gates"]["maximum_mc_z"]
    rng = np.random.default_rng(protocol["seeds"]["scalar"])
    scalar, integrals = [], []
    for s in protocol["scalar"]["snr"]:
        y = s + rng.standard_normal(protocol["scalar"]["draws_per_snr"])
        keep = y > 0
        yp = y[keep]
        p = norm.cdf(s)
        lam = norm.pdf(s) / p
        mean = s + lam
        var = 1 - s * lam - lam**2
        scalar.append({"snr": s, "p_keep_exact": p, "p_keep_mc": float(keep.mean()),
                       "p_keep_mc_z": float((keep.mean()-p)/np.sqrt(p*(1-p)/len(y))),
                       "positive_mean_exact": mean, "positive_mean_mc": float(yp.mean()),
                       "positive_mean_mc_z": float((yp.mean()-mean)/np.sqrt(var/len(yp))),
                       "positive_variance_exact": var, "positive_variance_mc": float(yp.var(ddof=1)),
                       "all_signed_mean_mc": float(y.mean()),
                       "all_signed_mean_mc_z": float((y.mean()-s)*np.sqrt(len(y)))})
        # Truncated score d log p(y|y>0,s)/ds = y-s-lambda(s).
        # Censored joint score: y-s when retained, -phi(s)/Phi(-s) when missing.
        mass = quad(lambda x: norm.pdf(x-s), 0, np.inf, epsabs=1e-12)[0]
        trunc_mass = quad(lambda x: norm.pdf(x-s)/p, 0, np.inf, epsabs=1e-12)[0]
        naive_score = quad(lambda x: (x-s)*norm.pdf(x-s)/p, 0, np.inf, epsabs=1e-12)[0]
        trunc_score = quad(lambda x: (x-s-lam)*norm.pdf(x-s)/p, 0, np.inf, epsabs=1e-12)[0]
        censored_score = quad(lambda x: (x-s)*norm.pdf(x-s), 0, np.inf, epsabs=1e-12)[0]-norm.pdf(s)
        integrals.append({"snr": s, "unconditional_positive_mass": mass,
                          "conditional_positive_mass": trunc_mass,
                          "naive_positive_score_expectation": naive_score,
                          "truncated_score_expectation": trunc_score,
                          "censored_joint_score_expectation": censored_score,
                          "max_closure_error": max(abs(mass-p), abs(trunc_mass-1), abs(naive_score-lam), abs(trunc_score), abs(censored_score))})

    config = protocol["paired"]
    f, sigma, n, reps = config["truth_flux"], config["sigma"], config["epochs"], config["realizations"]
    archived_flux = f + sigma * np.random.default_rng(protocol["seeds"]["archive"]).standard_normal(n)
    archive = archived_flux > 0
    rng = np.random.default_rng(protocol["seeds"]["paired"])
    y = f + sigma * rng.standard_normal((reps, n))
    masks = {"F": np.ones((reps, n), dtype=bool), "T": y > 0,
             "M": np.broadcast_to(archive, y.shape), "MT": archive & (y > 0)}
    estimates, counts, summaries = {}, {}, {}
    p = norm.cdf(f/sigma)
    lam = norm.pdf(f/sigma)/p
    tv = sigma**2 * (1 - (f/sigma)*lam - lam**2)
    for name, mask in masks.items():
        count = mask.sum(axis=1)
        estimates[name] = np.divide(np.where(mask, y, 0).sum(axis=1), count,
                                    out=np.full(reps, np.nan), where=count>0)
        counts[name] = count
        result = stats(estimates[name])
        ne = n if name in ["F", "T"] else int(archive.sum())
        target = f + sigma*lam if name in ["T", "MT"] else f
        k = np.arange(1, ne+1)
        theoretical_variance = (tv * np.sum(binom.pmf(k, ne, p)/k) / (1-(1-p)**ne)
                                if name in ["T", "MT"] else sigma**2/ne)
        result.update(expected_mean=target, mean_mc_z=(result["mean"]-target)/np.sqrt(theoretical_variance/result["n"]),
                      expected_sd=float(np.sqrt(theoretical_variance)),
                      expected_empty_probability=float((1-p)**ne) if name in ["T", "MT"] else 0.0,
                      mean_epochs=float(count.mean()), min_epochs=int(count.min()), max_epochs=int(count.max()))
        summaries[name] = result
    differences = {label: stats(estimates[a]-estimates[b]) for label,a,b in [
        ("T-F", "T", "F"), ("M-F", "M", "F"), ("MT-M", "MT", "M"), ("MT-T", "MT", "T")]}
    differences["interaction_(MT-M)-(T-F)"] = stats((estimates["MT"]-estimates["M"])-(estimates["T"]-estimates["F"]))
    gates = {
        "scalar_moments": all(abs(row[key]) < gate_z for row in scalar for key in ["p_keep_mc_z", "positive_mean_mc_z", "all_signed_mean_mc_z"]),
        "paired_means": all(abs(row["mean_mc_z"]) < gate_z for row in summaries.values()),
        "likelihood_normalization_and_score": max(r["max_closure_error"] for r in integrals) < protocol["gates"]["quadrature_absolute"],
        "all_failure_counts_reported": True,
    }
    np.savez_compressed(OUT / "paired-draws.npz", archived_flux=archived_flux, archived_mask=archive,
                        **{f"estimate_{k}": v for k,v in estimates.items()},
                        **{f"count_{k}": v for k,v in counts.items()})
    result = {"protocol_sha256": sha(OUT/"protocol.json"), "source_sha256": sha(__file__),
              "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "seconds": time.monotonic()-start},
              "scope": "Synthetic constant-signal flux means only; no RAISIN distance or dust estimate",
              "scalar": scalar, "likelihood_integrals": integrals,
              "paired": {"archived_retained": int(archive.sum()), "total_epochs": n, "truth_flux": f,
                         "summaries": summaries, "paired_differences": differences},
              "gates": gates, "all_gates_pass": all(gates.values())}
    save(OUT / "result.json", result)
    print(json.dumps(result, indent=2))
    assert result["all_gates_pass"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["freeze", "run"])
    args = parser.parse_args()
    freeze() if args.action == "freeze" else run()
