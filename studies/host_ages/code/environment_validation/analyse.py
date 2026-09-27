#!/usr/bin/env python3
"""Observed host-aperture checks and an independent-age calibrator comparison.

No data are distributed here. The CIGALE Bayesian summaries are measurements
conditional on that paper's SED model, not direct progenitor ages. Shared age
posterior/calibration covariance is unavailable, so errors-in-variables results
are explicitly sensitivity analyses. Downloaded and expanded tables stay local.
"""
from __future__ import annotations
import argparse, datetime, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize_scalar
from scipy.stats import norm, spearmanr, wilcoxon

ROOT = Path(__file__).resolve().parents[4]
DEFAULT_ARCHIVE = (
    Path.home()
    / ".local/share/cosmology-revalidation/archive/17487bf659fcbdeeea072221492bac14b04a0a85"
)
HERE = Path(__file__).resolve().parent
DESIGN = json.loads((HERE / "design.json").read_text())


def sha(p):
    h = hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda: f.read(2**20), b""):
            h.update(b)
    return h.hexdigest()


def gls(x, y, c):
    f = cho_factor(c, lower=True)
    ci_x = cho_solve(f, x)
    cov = np.linalg.inv(x.T @ ci_x)
    b = cov @ (x.T @ cho_solve(f, y))
    r = y - x @ b
    return b, cov, float(r @ cho_solve(f, r))


def aperture(a, b, label, rng):
    m = a.merge(
        b,
        on="SN" if "SN" in b else "Host",
        suffixes=("_local", "_other"),
        validate="one_to_one" if "SN" in b else "many_to_one",
    )
    host = m["Host_local"] if "Host_local" in m else m.Host
    groups = [np.flatnonzero(host.to_numpy() == h) for h in sorted(host.unique())]
    samples = [
        np.concatenate([groups[i] for i in rng.integers(0, len(groups), len(groups))])
        for _ in range(DESIGN["bootstrap_draws"])
    ]
    out = {"rows": len(m), "hosts": len(groups), "direction": label, "properties": {}}
    for p in ["AgeMW", "AV", "u-r", "logsSFR", "logM"]:
        scale = 1000 if p == "AgeMW" else 1
        d = (m[p + "_other"] - m[p + "_local"]).to_numpy() / scale
        med = float(np.median(d))
        draws = np.array([np.median(d[s]) for s in samples])
        err = np.hypot(m[p + "_err_other"], m[p + "_err_local"]).to_numpy() / scale
        out["properties"][p] = {
            "median": med,
            "cluster_bootstrap_95": np.quantile(draws, [0.025, 0.975]).tolist(),
            "nmad": float(1.4826 * np.median(abs(d - med))),
            "wilcoxon_p_ignoring_host_dependence": float(wilcoxon(d).pvalue),
            "fraction_abs_difference_exceeding_quadrature_error": float(
                np.mean(abs(d) > err)
            ),
            "error_limit": "quadrature assumes independent local/global errors, which is not justified by overlapping photometry",
        }
    return out, m


def cov_diagnostic(d):
    cols = ["AgeMW", "AV", "logM", "u-r", "logsSFR"]
    v = d[cols].to_numpy().copy()
    v[:, 0] /= 1000
    x = np.column_stack([np.ones(len(v)), v[:, 1:3]])
    r = v[:, 0] - x @ np.linalg.lstsq(x, v[:, 0], rcond=None)[0]
    return {
        "columns": ["age_Gyr"] + cols[1:],
        "pearson": np.corrcoef(v.T).tolist(),
        "spearman": spearmanr(v).statistic.tolist(),
        "age_sd_Gyr": float(np.std(v[:, 0], ddof=1)),
        "age_residual_sd_after_AV_and_mass_Gyr": float(np.std(r, ddof=1)),
        "age_fraction_variance_left_after_AV_mass": float(np.var(r) / np.var(v[:, 0])),
        "median_quoted_age_error_Gyr": float(np.median(d.AgeMW_err) / 1000),
        "scope": "between-object covariance of posterior summaries, not within-object age/dust posterior covariance",
    }


ALIASES = {"2008fv_comb": "2008fv", "1994drichmond": "1994d"}


def normalize(v):
    name = str(v).strip().lower().removeprefix("sn")
    return ALIASES.get(name, name)


def calibrators(local, phot, pantheon, cov, work, rng):
    loc = local.merge(
        phot[["SN", "RA", "Dec", "Dist"]], on="SN", validate="one_to_one"
    ).rename(columns={"RA": "host_table_SN_RA", "Dec": "host_table_SN_Dec"})
    loc["key"] = loc.SN.map(normalize)
    pan = pantheon.copy()
    pan["source_row"] = np.arange(len(pan))
    pan["key"] = pan.CID.map(normalize)
    assert loc.key.is_unique
    m = loc.merge(pan, on="key", validate="one_to_many")
    pcoord = SkyCoord(pan.RA.to_numpy() * u.deg, pan.DEC.to_numpy() * u.deg)
    lcoord = SkyCoord(
        loc.host_table_SN_RA.to_numpy() * u.deg,
        loc.host_table_SN_Dec.to_numpy() * u.deg,
    )
    idx, sep, _ = lcoord.match_to_catalog_sky(pcoord)
    new = []
    for i in np.flatnonzero(sep.arcsec <= 3):
        if loc.iloc[i]["key"] != pan.iloc[idx[i]]["key"]:
            new.append(
                {
                    "SN": str(loc.iloc[i].SN),
                    "CID": str(pan.iloc[idx[i]].CID),
                    "separation_arcsec": float(sep.arcsec[i]),
                }
            )
    assert not new, "New coordinate aliases need manual epoch validation before joining"
    m["coordinate_separation_arcsec"] = (
        SkyCoord(
            m.host_table_SN_RA.to_numpy() * u.deg,
            m.host_table_SN_Dec.to_numpy() * u.deg,
        )
        .separation(SkyCoord(m.RA.to_numpy() * u.deg, m.DEC.to_numpy() * u.deg))
        .arcsec
    )
    m.to_csv(work / "dustpedia-pantheon-crosswalk.csv", index=False)
    rejected = m.loc[
        m.coordinate_separation_arcsec > 3,
        ["SN", "CID", "coordinate_separation_arcsec"],
    ].to_dict("records")
    m = m[m.coordinate_separation_arcsec <= 3].copy()
    cal = m[m.IS_CALIBRATOR == 1].copy().reset_index(drop=True)
    ids = cal.source_row.to_numpy()
    c = cov[np.ix_(ids, ids)]
    y = (cal.m_b_corr - cal.CEPH_DIST).to_numpy()
    age = cal.AgeMW.to_numpy() / 1000
    av = cal.AV.to_numpy()
    mass = cal.global_logM.to_numpy()
    xbase = np.column_stack([np.ones(len(cal))])
    xa = np.column_stack([xbase, age - age.mean()])
    xfull = np.column_stack([xa, av - av.mean(), mass - mass.mean()])
    models = {"intercept": xbase, "age": xa, "age_AV_globalmass": xfull}
    fits = {}
    for name, x in models.items():
        b, bc, chi = gls(x, y, c)
        pred = np.zeros(len(y))
        logscore = 0.0
        square = []
        for host in sorted(cal.Host.unique()):
            te = np.flatnonzero(cal.Host.to_numpy() == host)
            tr = np.flatnonzero(cal.Host.to_numpy() != host)
            xt = x[tr]
            xv = x[te]
            ct = c[np.ix_(tr, tr)]
            cv = c[np.ix_(te, te)]
            cvt = c[np.ix_(te, tr)]
            bt, btcov, _ = gls(xt, y[tr], ct)
            fac = cho_factor(ct, lower=True)
            op = cho_solve(fac, cvt.T).T
            pred[te] = xv @ bt + op @ (y[tr] - xt @ bt)
            design = xv - op @ xt
            pc = cv - op @ cvt.T + design @ btcov @ design.T
            diff = y[te] - pred[te]
            pf = cho_factor(pc, lower=True)
            logscore += float(
                -0.5
                * (
                    len(te) * np.log(2 * np.pi)
                    + 2 * np.log(np.diag(pf[0])).sum()
                    + diff @ cho_solve(pf, diff)
                )
            )
            square.append(float(np.mean(diff**2)))
        fit = {
            "coefficients": b.tolist(),
            "coefficient_se": np.sqrt(np.diag(bc)).tolist(),
            "coefficient_covariance": bc.tolist(),
            "chi2": chi,
            "degrees_freedom": len(y) - x.shape[1],
            "leave_host_out_log_score": logscore,
            "leave_host_out_host_weighted_rmse_mag": float(np.sqrt(np.mean(square))),
        }
        if x.shape[1] > 1:
            se = float(np.sqrt(bc[1, 1]))
            delta = 0.03 / se
            fit.update(
                age_slope_mag_per_Gyr=float(b[1]),
                age_slope_fixed_age_95=[
                    float(b[1] - 1.96 * se),
                    float(b[1] + 1.96 * se),
                ],
                fixed_age_power_for_abs_003=float(
                    norm.cdf(-1.96 - delta) + norm.sf(1.96 - delta)
                ),
                fixed_age_80pct_detection_slope=float(
                    (norm.ppf(0.975) + norm.ppf(0.8)) * se
                ),
            )
            # Repeat synthetic correlated errors on this fixed observed design.
            estimates = (
                np.linalg.inv(x.T @ np.linalg.solve(c, x))
                @ x.T
                @ np.linalg.solve(
                    c, np.linalg.cholesky(c) @ rng.normal(size=(len(y), 10000))
                )
            )[1]
            fit["null_95_interval_coverage_10000"] = float(
                np.mean(abs(estimates / se) < 1.96)
            )
        fits[name] = fit
    # The same SN's age error is perfectly shared across its repeated light-curve rows.
    keys = sorted(cal.SN.unique())
    a = np.array([[sn == k for k in keys] for sn in cal.SN], float)
    ageerr = np.array([cal.loc[cal.SN == k, "AgeMW_err"].iloc[0] / 1000 for k in keys])
    agecov = (a * ageerr**2) @ a.T
    eiv = {}
    for name, x in [("age", xa), ("age_AV_globalmass", xfull)]:
        nuisance = np.delete(x, 1, axis=1)

        def nll(slope):
            cc = c + slope * slope * agecov
            yy = y - slope * x[:, 1]
            b, _, q = gls(nuisance, yy, cc)
            return 0.5 * (q + np.linalg.slogdet(cc)[1])

        opt = minimize_scalar(
            nll, bounds=(-0.3, 0.3), method="bounded", options={"xatol": 1e-12}
        )
        grid = np.linspace(-0.3, 0.3, 6001)
        delta = np.array([2 * (nll(s) - opt.fun) for s in grid])
        good = grid[delta <= 3.841459]
        eiv[name] = {
            "slope_mle": float(opt.x),
            "nominal_profile_95": [float(good.min()), float(good.max())],
            "grid_step": 0.0001,
            "scope": "Gaussian local-age-summary working likelihood with no population prior; no within-host or age/dust posterior covariance; interval coverage not certified",
        }
    # Shared-host pairs cancel the Cepheid distance and most global-host terms.
    siblings = []
    for host, g in cal.groupby("Host"):
        sns = g.SN.unique()
        if len(sns) < 2:
            continue
        for i, sn1 in enumerate(sns):
            for sn2 in sns[i + 1 :]:
                w = np.zeros(len(cal))
                for sign, sn in [(1, sn1), (-1, sn2)]:
                    rows = np.flatnonzero(cal.SN.to_numpy() == sn)
                    cs = c[np.ix_(rows, rows)]
                    iw = np.linalg.solve(cs, np.ones(len(rows)))
                    iw /= iw.sum()
                    w[rows] = sign * iw
                da = float(w @ age)
                dy = float(w @ y)
                e1 = cal.loc[cal.SN == sn1, "AgeMW_err"].iloc[0] / 1000
                e2 = cal.loc[cal.SN == sn2, "AgeMW_err"].iloc[0] / 1000
                siblings.append(
                    {
                        "host": host,
                        "first_SN": sn1,
                        "second_SN": sn2,
                        "first_minus_second_age_Gyr": da,
                        "age_difference_error_independent_summaries_Gyr": float(
                            np.hypot(e1, e2)
                        ),
                        "corrected_magnitude_difference": dy,
                        "magnitude_difference_se": float(np.sqrt(w @ c @ w)),
                        "predicted_difference_if_slope_minus003": -0.03 * da,
                    }
                )
    cal.to_csv(work / "dustpedia-calibrators.csv", index=False)
    return {
        "identifier_match_rows": len(m) + len(rejected),
        "coordinate_verified_match_rows": len(m),
        "unique_SNe": int(m.SN.nunique()),
        "verified_aliases": ALIASES,
        "coordinate_only_alias_candidates": new,
        "coordinate_conflicts": rejected,
        "maximum_verified_separation_arcsec": float(
            m.coordinate_separation_arcsec.max()
        ),
        "matches": sorted(m.SN.unique().tolist()),
        "noncalibrator_zHD_range": [
            float(m.loc[m.IS_CALIBRATOR == 0, "zHD"].min()),
            float(m.loc[m.IS_CALIBRATOR == 0, "zHD"].max()),
        ],
        "calibrator_rows": len(cal),
        "calibrator_SNe": int(cal.SN.nunique()),
        "calibrator_hosts": int(cal.Host.nunique()),
        "covariance_minimum_eigenvalue": float(np.linalg.eigvalsh(c).min()),
        "covariance_includes_Cepheid_errors": True,
        "statistical_scope": "fixed released distance covariance and CIGALE summary estimates; repeated light curves and Cepheid hosts retain full covariance; observational host follow-up selection is not modelled",
        "fits": fits,
        "gaussian_age_error_sensitivity": eiv,
        "calibrator_siblings": siblings,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    ap.add_argument("--work", type=Path, default=ROOT / ".work/environment-validation")
    ap.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / "studies/host_ages/results/environment_validation/dustpedia-summary.json",
    )
    args = ap.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    folder = (
        args.archive
        / "sources/updates/2026-09-20-ztf/extracted/kelsey2026-stag1765-supplement"
    )
    files = [
        folder / f
        for f in [
            "cigale_local_3kpc.txt",
            "cigale_local_1kpc.txt",
            "cigale_global.txt",
            "photometry_local_3kpc.txt",
        ]
    ]
    loc, one, glo, phot = [pd.read_csv(f, comment="#") for f in files]
    assert (
        len(loc) == 90
        and loc.SN.is_unique
        and len(glo) == 78
        and glo.Host.is_unique
        and len(one) == 41
        and one.SN.is_unique
    )
    assert set(loc.Host) == set(glo.Host) and set(one.SN) <= set(loc.SN)
    assert not loc[["AgeMW", "AgeMW_err", "AV", "AV_err", "logM"]].isna().any().any()
    rng = np.random.default_rng(DESIGN["seed"])
    global_ap, _ = aperture(loc, glo, "global minus local 3 kpc", rng)
    one_ap, _ = aperture(
        loc[loc.SN.isin(one.SN)], one, "local 1 kpc minus local 3 kpc", rng
    )
    loc = loc.merge(
        glo[["Host", "logM"]].rename(columns={"logM": "global_logM"}),
        on="Host",
        validate="many_to_one",
    )
    pf = ROOT / "data/distances/Pantheon+SH0ES.dat"
    cf = ROOT / "data/distances/Pantheon+SH0ES_STAT+SYS.cov"
    pan = pd.read_csv(pf, sep=r"\s+")
    c = np.loadtxt(cf, skiprows=1).reshape(len(pan), len(pan))
    asym = float(np.max(abs(c - c.T)))
    assert asym < 3.1e-8
    c = (c + c.T) / 2
    readme = (
        args.archive
        / "sources/repos/PantheonPlusSH0ES__DataRelease@7fc6805/Pantheon+_Data/4_DISTANCES_AND_COVAR/README"
    )
    result = {
        "schema": 1,
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "study": "independent UV-IR host measurement and corrected-distance comparison",
        "design_sha256": sha(HERE / "design.json"),
        "code_sha256": sha(__file__),
        "inputs": [
            {"path": str(f), "sha256": sha(f)} for f in files + [pf, cf, readme]
        ],
        "input_covariance_max_asymmetry_mag2": asym,
        "input_covariance_symmetrized": True,
        "local_global": global_ap,
        "one_vs_three_kpc": one_ap,
        "local_covariation": cov_diagnostic(loc),
        "calibrator_comparison": calibrators(loc, phot, pan, c, args.work, rng),
        "limitations": [
            "CIGALE age is a model-dependent mass-weighted stellar age, not progenitor delay.",
            "Host aperture attenuation is not SN line-of-sight extinction.",
            "No joint age/dust/mass posterior or cross-aperture covariance is released.",
            "DustPedia follow-up and overlap are selected, nearby and not independent of Pantheon distances.",
            "This does not measure redshift evolution or validate a cosmological age correction.",
        ],
    }
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(args.output)
    print(
        json.dumps(
            {
                "aperture": result["local_global"],
                "calibrators": result["calibrator_comparison"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
