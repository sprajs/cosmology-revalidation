#!/usr/bin/env python3
"""Common-support native BBC response; never use bin-centered MURES as drift."""

import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import expit
from astropy.cosmology import FlatLambdaCDM
from common import RESULTS, HERE, sha
from scaled_bbc import SCALE, EDGES, fitres

COSMO = FlatLambdaCDM(H0=70, Om0=0.315, Tcmb0=0)
CASES = ["nominal", "age_nominal", "age_literal", "age_retained"]


def reconstruct(frame, pars):
    """Native no-offset mB normalization; omit only the global constant."""
    alpha, beta, gamma = [pars[k]["value"] for k in ["alpha0", "beta0", "gamma0"]]
    host = gamma * (0.5 - expit((frame.HOST_LOGMASS - 10) / 0.001))
    return (
        -2.5 * np.log10(frame.x0)
        + alpha * frame.x1
        - beta * frame.c
        - host
        - frame.biasCor_mu
    )


def summarize(frame, residual, weights):
    age = frame.SIM_HOSTLIB_SN_age.to_numpy()
    z = frame.zHD.to_numpy()
    bins = np.digitize(z, EDGES) - 1
    residual = np.asarray(residual)
    weights = np.asarray(weights)
    centered = residual - np.average(residual, weights=weights)
    means = []
    for j in range(len(EDGES) - 1):
        keep = bins == j
        means.append(
            {
                "z_min": float(EDGES[j]),
                "z_max": float(EDGES[j + 1]),
                "n": int(keep.sum()),
                "mean_mag": (
                    float(np.average(centered[keep], weights=weights[keep]))
                    if keep.sum() >= 10
                    else None
                ),
            }
        )
    allowed = [j for j in range(len(EDGES) - 1) if (bins == j).sum() >= 10]
    keep = np.isin(bins, allowed)
    slope = None
    if keep.sum() >= 30:
        x = np.column_stack(
            [age[keep]] + [(bins[keep] == j).astype(float) for j in allowed]
        )
        w = weights[keep]
        if np.linalg.matrix_rank(x) == x.shape[1]:
            b = np.linalg.lstsq(
                x * np.sqrt(w[:, None]), residual[keep] * np.sqrt(w), rcond=None
            )[0]
            slope = float(b[0])
    return {
        "n": len(frame),
        "weighted_effective_n": float(weights.sum() ** 2 / np.sum(weights**2)),
        "slope_n": int(keep.sum()),
        "slope_weighted_effective_n": (
            float(weights[keep].sum() ** 2 / np.sum(weights[keep] ** 2))
            if keep.any()
            else 0.0
        ),
        "age_range_Gyr": [float(np.min(age)), float(np.max(age))],
        "slope_mag_per_Gyr": slope,
        "raw_weighted_global_offset_mag": float(np.average(residual, weights=weights)),
        "centered_redshift_means": means,
    }


def analyze(cases, output=None):
    frames = {}
    audits = {}
    for name in CASES:
        r = cases[name]
        if not r["graceful"]:
            continue
        path = Path(r["folder"]) / "bbc.FITRES"
        assert sha(path) == r["outputs"]["bbc.FITRES"]
        f = fitres(path)
        assert f.index.is_unique and (f.MUERR > 0).all()
        frames[name] = f
        raw = reconstruct(f, r["parameters"])
        offset = np.median(raw - f.MU)
        maxerr = float(np.max(abs(raw - offset - f.MU)))
        mures = float(np.max(abs(f.MU - f.MUMODEL - f.MURES - f.M0DIF)))
        # Header nuisance values and final distance columns are rounded natively.
        assert maxerr < 3e-4 and mures < 1.5e-4, (name, maxerr, mures)
        audits[name] = {
            "native_Tripp_shape_max_error_mag": maxerr,
            "fitted_global_constant_mag": float(offset),
            "MU_MUMODEL_MURES_M0DIF_identity_max_error_mag": mures,
            "n": len(f),
        }
    out = {
        "scope": "Conditional mock intervention; fitted ages are simulated W22 delay truths, not observations.",
        "native_identity_audit": audits,
        "contrasts": {},
        "closure": {},
        "native_paired_support": False,
        "analyzable_slope_support": False,
        "estimand_definitions": {
            "frozen_nominal": "Nominal training and nominal fitted nuisance values. Native map geometry is evaluation-cohort dependent unless the separate common-input geometry audit passes.",
            "age_nominal": "Nominal training, native nuisance mode of this case; same map/support as the externally frozen case.",
            "age_literal": "Age-generated training, native target subtracts true injected SIM_gammaDM.",
            "age_retained": "Age-generated training with SIM_gammaDM zeroed only in training target; also changes global gamma centering and potentially MUCOV.",
        },
    }
    if "nominal" not in frames:
        out["status"] = "unsupported_nominal_map"
        if output:
            out["code_sha256"] = sha(__file__)
            Path(output).write_text(json.dumps(out, indent=2) + "\n")
        return out
    nom = frames["nominal"]
    out["closure"]["nominal"] = summarize(nom, nom.MU - nom.SIM_DLMAG, 1 / nom.MUERR**2)
    if not cases["nominal"].get("native_fit_identified", False):
        out["closure"]["nominal"]["slope_mag_per_Gyr"] = None
        out["closure"]["nominal"][
            "reason"
        ] = "Native nuisance or covariance validity gate failed."
    complete = nom.index
    for f in frames.values():
        complete = complete.intersection(f.index)
    if len(frames) != len(CASES):
        complete = complete[:0]
    for name, frozen in [
        ("age_nominal", True),
        ("age_nominal", False),
        ("age_literal", False),
        ("age_retained", False),
    ]:
        if name not in frames:
            continue
        label = "frozen_nominal" if frozen else name
        a = frames[name]
        for support, ids in [
            ("pair", nom.index.intersection(a.index)),
            ("all_cases", complete),
        ]:
            if not len(ids):
                continue
            n = nom.loc[ids]
            f = a.loc[ids]
            for col in ["SIM_HOSTLIB_SN_age", "SIM_DLMAG", "SIM_LIBID", "SIM_ZCMB"]:
                assert np.max(abs(n[col] - f[col])) < 2e-6, (name, col)
            mu = reconstruct(f, cases["nominal"]["parameters"]) if frozen else f.MU
            # Different global offsets are unidentifiable and removed exactly once.
            delta = (mu - f.SIM_DLMAG) - (n.MU - n.SIM_DLMAG)
            weights = 1 / n.MUERR**2
            r = summarize(n, delta, weights)
            r["native_nuisance_identified"] = bool(
                cases["nominal"].get("native_fit_identified", False)
                and (frozen or cases[name].get("native_fit_identified", False))
            )
            if not r["native_nuisance_identified"]:
                r["slope_mag_per_Gyr"] = None
                r["reason"] = (
                    "Native nuisance design rank or covariance validity gate failed."
                )
            truez_n = n.SIM_DLMAG + 5 * np.log10(
                COSMO.luminosity_distance(n.zHD).value
                / COSMO.luminosity_distance(n.SIM_ZCMB).value
            )
            truez_a = f.SIM_DLMAG + 5 * np.log10(
                COSMO.luminosity_distance(f.zHD).value
                / COSMO.luminosity_distance(f.SIM_ZCMB).value
            )
            r["observed_z_truth_sensitivity"] = summarize(
                n, (mu - truez_a) - (n.MU - truez_n), weights
            )
            if not r["native_nuisance_identified"]:
                r["observed_z_truth_sensitivity"]["slope_mag_per_Gyr"] = None
            r["fraction_of_injected_slope"] = (
                r["slope_mag_per_Gyr"] / -0.03
                if r["slope_mag_per_Gyr"] is not None
                else None
            )
            r["intervention_shift_age_identity_max_mag"] = float(
                np.max(abs(f.SIM_gammaDM + 0.03 * (f.SIM_HOSTLIB_SN_age - 3)))
            )
            r["common_occurrence_ids"] = ids.to_list()
            out["contrasts"][label + "_" + support] = r
    out["native_paired_support"] = bool(out["contrasts"])
    out["analyzable_slope_support"] = any(
        v["slope_mag_per_Gyr"] is not None for v in out["contrasts"].values()
    )
    out["status"] = (
        "supported" if out["analyzable_slope_support"] else "no_analyzable_paired_slope"
    )
    if output:
        out["code_sha256"] = sha(__file__)
        out["analysis_design_sha256"] = sha(HERE / "scaled-analysis-design.json")
        Path(output).write_text(json.dumps(out, indent=2) + "\n")
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--shards", type=int, required=True)
    p.add_argument("--common", action="store_true")
    p.add_argument("--fixed", action="store_true")
    p.add_argument("--highmass-gauge", action="store_true")
    a = p.parse_args()
    suffix = (
        ("-common" if a.common or a.highmass_gauge else "")
        + ("-fixed" if a.fixed else "")
        + ("-highmass-gauge" if a.highmass_gauge else "")
    )
    r = json.loads(
        (RESULTS / (f"scaled-bbc-k{a.shards:03d}" + suffix + ".json")).read_text()
    )
    out = analyze(
        r["cases"], RESULTS / (f"scaled-response-k{a.shards:03d}" + suffix + ".json")
    )
    print(
        json.dumps(
            {k: v["slope_mag_per_Gyr"] for k, v in out["contrasts"].items()}, indent=2
        )
    )


if __name__ == "__main__":
    main()
