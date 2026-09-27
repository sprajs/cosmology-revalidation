#!/usr/bin/env python3
"""New SALT/environment predictive diagnostic from released ZTF fitted observables.

This is deliberately not a BayeSN-result reproduction, dust-mechanism fit, BBC
simulation or cosmology-ready Hubble diagram. Table cuts and every failed row
are saved; host groups remain together in validation folds.
"""
from __future__ import annotations
import argparse, datetime, hashlib, json, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from scipy.stats import norm
from iminuit import Minuit
from analyse import ROOT, HERE, DEFAULT_ARCHIVE, DESIGN, sha

C_LIGHT = 299792.458
COSMO = FlatLambdaCDM(H0=70, Om0=0.3)
MODELS = {
    "width_colour": ["intercept", "x1", "c", "z"],
    "mass": ["intercept", "x1", "c", "z", "mass_step"],
    "mass_localcolour": ["intercept", "x1", "c", "z", "mass_step", "localcolour"],
    "nonlinear": ["intercept", "x1", "c", "z", "mass_step", "localcolour", "c2", "x12"],
}


def make_matrix(d, names):
    columns = {
        "intercept": np.ones(len(d)),
        "x1": d.x1.to_numpy(),
        "c": d.c.to_numpy(),
        "z": d.redshift.to_numpy() - 0.035,
        "mass_step": norm.cdf(
            (d.mass_global.to_numpy() - 10)
            / np.maximum(d.mass_err_global.to_numpy(), 1e-6)
        ),
        "localcolour": d.restframe_gz_local.to_numpy() - 1,
        "c2": d.c.to_numpy() ** 2,
        "x12": d.x1.to_numpy() ** 2,
    }
    return np.column_stack([columns[k] for k in names])


def variance(d, b, names):
    param = dict(zip(names, b))
    dx = np.full(len(d), param["x1"])
    dc = np.full(len(d), param["c"])
    if "x12" in param:
        dx += 2 * param["x12"] * d.x1.to_numpy()
    if "c2" in param:
        dc += 2 * param["c2"] * d.c.to_numpy()
    g = np.column_stack([np.ones(len(d)), -dx, -dc])
    cov = np.stack(d.lc_cov.to_numpy())
    v = np.einsum("ni,nij,nj->n", g, cov, g) + d.mu_z_variance.to_numpy()
    if "localcolour" in param:
        v += param["localcolour"] ** 2 * d.restframe_gz_err_local.to_numpy() ** 2
    if "mass_step" in param:
        p = norm.cdf(
            (d.mass_global.to_numpy() - 10)
            / np.maximum(d.mass_err_global.to_numpy(), 1e-6)
        )
        v += param["mass_step"] ** 2 * p * (1 - p)
    return v


def fit(d, names, matrix_function=make_matrix, variance_function=variance):
    x = matrix_function(d, names)
    y = d.y.to_numpy()
    start = np.linalg.lstsq(x, y, rcond=None)[0]

    def cost(*par):
        b = np.array(par[:-1])
        sig = par[-1]
        v = variance_function(d, b, names) + sig**2
        if np.any(v <= 0):
            return 1e100
        return float(np.sum((y - x @ b) ** 2 / v + np.log(2 * np.pi * v)))

    m = Minuit(cost, *start, 0.12, name=names + ["intrinsic_sigma"])
    m.errordef = 1
    m.limits["intrinsic_sigma"] = (0.00001, 1)
    m.tol = 1e-5
    m.strategy = 2
    m.migrad(ncall=30000)
    m.hesse()
    if not m.valid:
        m.migrad(ncall=30000)
        m.hesse()
    if not m.valid or not m.fmin.has_posdef_covar:
        raise RuntimeError(str(m.fmin))
    b = np.array([m.values[k] for k in names])
    bc = np.array(m.covariance)[:-1, :-1]
    sig = float(m.values["intrinsic_sigma"])
    v = variance_function(d, b, names) + sig**2
    return (
        {
            "parameters": dict(zip(names, b.tolist())),
            "standard_errors": {n: float(m.errors[n]) for n in names},
            "intrinsic_sigma_mag": sig,
            "minus2loglike": float(m.fval),
            "coefficient_covariance": bc.tolist(),
            "valid_minimum": m.valid,
            "edm": float(m.fmin.edm),
            "residual_std_mag": float(np.std(y - x @ b)),
            "residual_nmad_mag": float(
                1.4826 * np.median(abs((y - x @ b) - np.median(y - x @ b)))
            ),
        },
        b,
        bc,
        sig,
    )


def groups(d):
    co = SkyCoord(d.ra_host.to_numpy() * u.deg, d.dec_host.to_numpy() * u.deg)
    i, j, _, _ = search_around_sky(co, co, 3 * u.arcsec)
    parent = np.arange(len(d))

    def root(x):
        while parent[x] != x:
            x = parent[x]
        return x

    for a, b in zip(i, j):
        parent[root(b)] = root(a)
    group = [
        str(min(d.ztfname.iloc[k] for k in range(len(d)) if root(k) == root(t)))
        for t in range(len(d))
    ]
    return group


def sample(sn, glo, loc, master):
    assert (
        sn.ztfname.is_unique
        and glo.ztfname.is_unique
        and loc.ztfname.is_unique
        and master.ztfname.is_unique
    )
    d = sn.merge(glo, on="ztfname", validate="one_to_one").merge(
        loc, on="ztfname", suffixes=("_global", "_local"), validate="one_to_one"
    )
    mask = pd.DataFrame(index=d.index)
    mask["masterlist"] = d.ztfname.isin(master.ztfname)
    mask["host_redshift"] = d.source.eq("z_gal")
    mask["redshift_range"] = d.redshift.between(0.01, 0.06)
    mask["coverage"] = d.lccoverage_flag.eq(1)
    mask["fitquality"] = d.fitquality_flag.eq(1)
    mask["SALT_cuts"] = (
        (d.fitprob > 1e-7)
        & d.x1.between(-3, 3)
        & (d.x1_err < 1)
        & d.c.between(-0.2, 0.8)
        & (d.c_err < 0.1)
        & (d.t0_err < 1)
        & (d.x0 > 0)
        & (d.x0_err > 0)
    )
    cols = [
        "mass_global",
        "mass_err_global",
        "restframe_gz_global",
        "restframe_gz_err_global",
        "mass_local",
        "mass_err_local",
        "restframe_gz_local",
        "restframe_gz_err_local",
        "ra_host",
        "dec_host",
    ]
    mask["finite_hosts"] = np.isfinite(d[cols]).all(axis=1)
    mask["normal_or_91T"] = d.sub_type.isin(["norm", "91T"])
    d["selected"] = mask.all(axis=1)
    audit = pd.concat(
        [d[["ztfname", "source", "sub_type"]], mask, d[["selected"]]], axis=1
    )
    cur = np.ones(len(d), bool)
    counts = {}
    for k in mask:
        cur &= mask[k].to_numpy()
        counts[k] = int(cur.sum())
    return d[d.selected].copy().reset_index(drop=True), audit, counts


def prepare(d):
    d = d.copy()
    coord = SkyCoord(d.ra.to_numpy() * u.deg, d.dec.to_numpy() * u.deg)
    dip = SkyCoord(l=264.021 * u.deg, b=48.253 * u.deg, frame="galactic")
    projection = np.cos(coord.separation(dip).rad)
    beta = 369.82 / C_LIGHT
    d["zCMB"] = (1 + d.redshift.to_numpy()) * np.sqrt(1 - beta**2) / (
        1 - beta * projection
    ) - 1
    mu = (
        5
        * np.log10(
            COSMO.comoving_distance(d.zCMB.to_numpy()).value
            * (1 + d.redshift.to_numpy())
        )
        + 25
    )
    d["y"] = -2.5 * np.log10(d.x0.to_numpy()) - mu
    z = d.zCMB.to_numpy()
    dz = 1e-6
    deriv = (
        5
        * np.log10(
            COSMO.comoving_distance(z + dz).value * (1 + d.redshift.to_numpy() + dz)
        )
        - 5
        * np.log10(
            COSMO.comoving_distance(z - dz).value * (1 + d.redshift.to_numpy() - dz)
        )
    ) / (2 * dz)
    d["mu_z_variance"] = deriv**2 * (
        d.redshift_err.to_numpy() ** 2 + (250 / C_LIGHT) ** 2
    )
    n = len(d)
    cov = np.zeros((n, 3, 3))
    dm = -2.5 / np.log(10) / d.x0.to_numpy()
    cov[:, 0, 0] = (dm * d.x0_err.to_numpy()) ** 2
    cov[:, 1, 1] = d.x1_err.to_numpy() ** 2
    cov[:, 2, 2] = d.c_err.to_numpy() ** 2
    cov[:, 0, 1] = cov[:, 1, 0] = dm * d.cov_x0_x1.to_numpy()
    cov[:, 0, 2] = cov[:, 2, 0] = dm * d.cov_x0_c.to_numpy()
    cov[:, 1, 2] = cov[:, 2, 1] = d.cov_x1_c.to_numpy()
    assert np.linalg.eigvalsh(cov).min() > 0
    d["lc_cov"] = list(cov)
    d["host_group"] = groups(d)
    d["fold"] = [
        int(hashlib.sha256((str(DESIGN["seed"]) + v).encode()).hexdigest()[:8], 16) % 5
        for v in d.host_group
    ]
    return d


def analyse(d, do_cv=True):
    fits = {}
    scores = {}
    predictions = {}
    for label, names in MODELS.items():
        f, b, bc, sig = fit(d, names)
        fits[label] = f
        if do_cv:
            score = np.zeros(len(d))
            pred = np.zeros(len(d))
            var = np.zeros(len(d))
            foldfit = []
            for fold in range(5):
                tr = d.fold.to_numpy() != fold
                te = ~tr
                dt = d[te]
                ff, bb, cc, ss = fit(d[tr], names)
                xx = make_matrix(dt, names)
                pred[te] = xx @ bb
                var[te] = (
                    variance(dt, bb, names)
                    + ss**2
                    + np.einsum("ni,ij,nj->n", xx, cc, xx)
                )
                r = dt.y.to_numpy() - pred[te]
                score[te] = -0.5 * (np.log(2 * np.pi * var[te]) + r * r / var[te])
                foldfit.append(
                    {
                        "fold": fold,
                        "training_count": int(tr.sum()),
                        "test_count": int(te.sum()),
                        "fit_valid": ff["valid_minimum"],
                    }
                )
            scores[label] = score
            predictions[label] = pred
            f["heldout"] = {
                "log_predictive_score": float(score.sum()),
                "rmse_mag": float(np.sqrt(np.mean((d.y.to_numpy() - pred) ** 2))),
                "folds": foldfit,
            }
    if do_cv:
        rng = np.random.default_rng(DESIGN["seed"])
        pairs = {}
        for label, base in [
            ("mass", "width_colour"),
            ("mass_localcolour", "mass"),
            ("nonlinear", "mass_localcolour"),
        ]:
            delta = scores[label] - scores[base]
            group = np.array(d.host_group)
            unique = np.unique(group)
            unit = np.array([delta[group == g].sum() for g in unique])
            boot = rng.choice(unit, (2000, len(unit)), replace=True).sum(axis=1)
            pairs[label + "_minus_" + base] = {
                "delta_log_predictive_score": float(delta.sum()),
                "host_bootstrap_95_conditional_on_fitted_folds": np.quantile(
                    boot, [0.025, 0.975]
                ).tolist(),
                "host_count": len(unique),
                "limit": "bootstrap treats held-out scores conditional on the fitted overlapping training folds; not complete training or calibration uncertainty",
            }
    else:
        pairs = {}
    return (
        {
            "n": len(d),
            "host_groups": int(d.host_group.nunique()),
            "fold_counts": d.fold.value_counts().sort_index().to_dict(),
            "fits": fits,
            "predictive_comparisons": pairs,
        },
        scores,
        predictions,
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    p.add_argument("--work", type=Path, default=ROOT / ".work/environment-validation")
    p.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / "studies/host_ages/results/environment_validation/ztf-summary.json",
    )
    args = p.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    base = args.archive / "sources/updates/2026-09-20-ztf"
    folder = base / "extracted/ztfsniadr2_lite/tables"
    paths = [
        folder / k
        for k in ["snia_data.csv", "globalhost_data.csv", "localhost_data.csv"]
    ] + [base / "originals/Ginolin25ab_masterlist--903d65d.csv"]
    sn, glo, loc, master = [pd.read_csv(f) for f in paths]
    d, audit, counts = sample(sn, glo, loc, master)
    audit.to_csv(args.work / "ztf-selection.csv", index=False)
    d = prepare(d)
    mainresult, scores, predictions = analyse(d)
    out = {
        "schema": 1,
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "code_sha256": sha(__file__),
        "shared_code_sha256": sha(HERE / "analyse.py"),
        "design_sha256": sha(HERE / "design.json"),
        "inputs": [{"path": str(f), "sha256": sha(f)} for f in paths],
        "selection_cumulative": counts,
        "selected_names_sha256": hashlib.sha256(
            "\n".join(sorted(d.ztfname)).encode()
        ).hexdigest(),
        "models": MODELS,
        "primary": mainresult,
        "sensitivity_z_at_least_002": analyse(
            d[d.redshift >= 0.02].reset_index(drop=True), False
        )[0],
        "sensitivity_colour_below_01": analyse(
            d[d.c < 0.1].reset_index(drop=True), False
        )[0],
        "assumptions": {
            "redshift": "release heliocentric frame verified in Rigault et al. 2409.04346v2 Sec 4.5/Table 3; CMB dipole 369.82 km/s at Galactic (264.021,48.253) degrees; exact Lorentz Doppler factor; no local flow reconstruction",
            "reference_cosmology": "flat LCDM Om=.3 H0=70; intercept and linear redshift nuisance are fitted",
            "variance": "first-order SALT covariance, 250 km/s velocity floor, fitted scatter, local-colour Gaussian summary error, soft mass mixture variance approximation",
            "environment": "local and global host mass/rest-frame g-z; none are measured age or SN extinction",
        },
        "limitations": [
            "Not the exact unpublished 932/929 author mask; smaller explicit host-z sample.",
            "Not a BayeSN or BBC fit; released SALT parameter fits are starting observables, not detector pixels.",
            "Neither mass nor colour identifies age, dust or explosion channel.",
            "Calibration and cross-SN velocity correlations are missing; no precision cosmological claim.",
            "Training and host follow-up selection are not modelled.",
            "C<0.1 restricts SALT colour, not actual low dust extinction.",
        ],
    }
    export = d.drop(columns=["lc_cov"]).copy()
    for k in scores:
        export[k + "_heldout_logscore"] = scores[k]
        export[k + "_heldout_prediction"] = predictions[k]
    export.to_csv(args.work / "ztf-analysis-ledger.csv", index=False)
    args.output.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    print(args.output)
    print(
        json.dumps(
            {
                "counts": counts,
                "fits": mainresult["fits"],
                "comparisons": mainresult["predictive_comparisons"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
