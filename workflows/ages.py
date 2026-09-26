"""Extract host ages and fit explicit Gaussian-summary regression models."""

import zipfile
import numpy as np
import pandas as pd
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import minimize
from lib.paths import DATA
from lib.ages import read_age, mu, wls

DEFAULTS = {"overlap": "G11_first", "zmin": 0.06, "zmax": 0.42, "latent_age_fit": True}


def run(out, cfg):
    if (
        cfg["overlap"] not in {"G11_first", "R19_first"}
        or not 0 < cfg["zmin"] < cfg["zmax"]
    ):
        raise ValueError("Invalid overlap policy or redshift interval")
    # Re-extract exact bytes from the archived supplement; do not parse PDFs.
    extracted = out / "extracted"
    extracted.mkdir()
    with zipfile.ZipFile(DATA / "ages/supplement.zip") as archive:
        for i in [1, 2]:
            matches = [
                n for n in archive.namelist() if n.split("/")[-1] == f"table{i}.dat"
            ]
            if len(matches) != 1:
                raise ValueError("Ambiguous supplementary table")
            raw = archive.read(matches[0])
            if raw != (DATA / f"ages/table{i}.dat").read_bytes():
                raise ValueError("Supplementary extraction differs from frozen table")
            (extracted / f"table{i}.dat").write_bytes(raw)
    frames = [
        read_age(extracted / f"table{i}.dat", label)
        for i, label in [(1, "G11"), (2, "R19")]
    ]
    ages = pd.concat(frames, ignore_index=True)
    pp = pd.read_csv(
        DATA / "distances/Pantheon+SH0ES.dat", sep=r"\s+", dtype={"CID": str}
    )
    pp["pp_row"] = np.arange(len(pp))
    joined = ages.merge(
        pp[pp.IDSURVEY == 1],
        on="CID",
        how="left",
        indicator=True,
        validate="many_to_one",
    )
    joined.to_csv(out / "all_age_crosswalk.csv", index=False)
    matched = joined[joined._merge == "both"].drop(columns="_merge").copy()
    matched["mu_model"] = mu(matched.zHD, matched.zHEL)
    matched["hr_corrected"] = matched.MU_SH0ES - matched.mu_model
    matched["hr_no_bias"] = matched.hr_corrected + matched.biasCor_m_b
    matched["hr_tripp"] = (
        matched.mB + 0.148 * matched.x1 - 3.112 * matched.c + 19.253 - matched.mu_model
    )
    selected = matched.drop_duplicates(
        "CID", keep="first" if cfg["overlap"] == "G11_first" else "last"
    )
    selected = selected[
        (selected.zHD > cfg["zmin"]) & (selected.zHD < cfg["zmax"])
    ].copy()
    selected.to_csv(out / "selected.csv", index=False)
    if len(selected) < 10:
        raise ValueError("Too few matched objects")
    wls_results = []
    for column in ["hr_corrected", "hr_no_bias", "hr_tripp"]:
        b, c, chi = wls(
            selected.age.to_numpy(),
            selected[column].to_numpy(),
            selected.MU_SH0ES_ERR_DIAG.to_numpy(),
        )
        wls_results.append(
            {
                "outcome": column,
                "slope": b[1],
                "slope_se": np.sqrt(c[1, 1]),
                "chi2": chi,
                "dof": len(selected) - 2,
            }
        )
    result = {
        "source_rows": {"G11": len(frames[0]), "R19": len(frames[1])},
        "unmatched_age_rows": int((joined._merge != "both").sum()),
        "matched_rows_before_deduplication": len(matched),
        "selected_unique_objects": len(selected),
        "wls": wls_results,
        "scope": "Host-age summary statistics joined to released SDSS distances. WLS ignores age errors. The latent Gaussian fit is a distributional sensitivity, not original age PDFs, LINMIX, or selection-normalized inference.",
    }
    if not cfg["latent_age_fit"]:
        return result
    raw = np.loadtxt(DATA / "distances/Pantheon+SH0ES_STAT+SYS.cov")
    n = int(raw[0])
    if n != len(pp) or len(raw) != 1 + n * n:
        raise ValueError("Covariance row order mismatch")
    ids = selected.pp_row.to_numpy(dtype=int)
    cov = raw[1:].reshape(n, n)[np.ix_(ids, ids)]
    cov = (cov + cov.T) / 2
    age, error = selected.age.to_numpy(), selected.age_err.to_numpy()
    bounds = [(0, 14), (-4, 2.6), (-3, 3), (-0.2, 0.15), (-9, 0)]
    latent = []
    for column in ["hr_corrected", "hr_no_bias"]:
        y = selected[column].to_numpy()

        def nll(theta):
            mean, logtau, intercept, slope, logscatter = theta
            tau2 = np.exp(2 * logtau)
            av = tau2 + error**2
            posterior_mean = mean + tau2 / av * (age - mean)
            posterior_variance = tau2 * error**2 / av
            matrix = cov + np.diag(
                np.exp(2 * logscatter) + slope * slope * posterior_variance
            )
            chol = cholesky(matrix, lower=True)
            white = solve_triangular(
                chol, y - intercept - slope * posterior_mean, lower=True
            )
            return (
                0.5 * np.sum(np.log(av) + (age - mean) ** 2 / av)
                + np.log(np.diag(chol)).sum()
                + 0.5 * (white @ white)
            )

        fits = [
            minimize(
                nll,
                [age.mean(), np.log(t), y.mean() - b * age.mean(), b, np.log(0.07)],
                method="L-BFGS-B",
                bounds=bounds,
                options={"maxiter": 2000, "ftol": 1e-12, "gtol": 1e-6},
            )
            for t, b in [(1, -0.03), (2, -0.01), (0.3, -0.06)]
        ]
        good = [x for x in fits if x.success and np.isfinite(x.fun)]
        if not good:
            raise RuntimeError("No latent-age optimizer converged")
        best = min(good, key=lambda x: x.fun)
        latent.append(
            {
                "outcome": column,
                "theta": best.x.tolist(),
                "nll": float(best.fun),
                "slope": float(best.x[3]),
                "successful_starts": len(good),
                "boundary": any(
                    min(abs(best.x[i] - b[0]), abs(best.x[i] - b[1])) < 1e-4
                    for i, b in enumerate(bounds)
                ),
            }
        )
    result["gaussian_latent_age_mle"] = latent
    return result
