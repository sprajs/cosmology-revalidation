#!/usr/bin/env python3
"""Small low-z Dovekie/TITAN overlap: conditional slopes and achievable precision."""
import argparse, datetime, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
from scipy.linalg import cho_factor, cho_solve, solve_triangular
from scipy.stats import norm
from analyse import ROOT, HERE, DEFAULT_ARCHIVE, sha, normalize, gls

sys.path.insert(0, str(ROOT))
from lib.records import fitres


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument(
        "--titan",
        type=Path,
        default=ROOT
        / ".work/population-transport/titan-author/host_props_with_SN_age_good_Mar18.csv",
    )
    args = parser.parse_args()
    folder = args.archive / "sources/repos/des-science__DES-SN5YR/4_DISTANCES_COVMAT"
    hp = folder / "DES-Dovekie_HD.csv"
    cp = folder / "STAT+SYS.npz"
    tp = args.titan
    hd = fitres(hp).reset_index()
    t = pd.read_csv(tp)
    hd["source_row"] = np.arange(len(hd))
    hd["key"] = hd.CID.map(normalize)
    t["key"] = t.transient.map(normalize)
    assert hd.key.is_unique and t.key.is_unique
    m = hd.merge(t, on="key", validate="one_to_one")
    nmatch = len(m)
    valid = (m.d_dlr < 4) & np.isfinite(
        m[
            [
                "mass_weighted_age_50",
                "mass_weighted_age_16",
                "mass_weighted_age_84",
                "dust:Av_50",
                "stellar_mass_50",
            ]
        ]
    ).all(axis=1)
    m["selected"] = valid
    m.to_csv(
        ROOT / ".work/environment-validation/titan-dovekie-crosswalk.csv", index=False
    )
    m = m[valid].reset_index(drop=True)
    assert (m.mass_weighted_age_16 <= m.mass_weighted_age_50).all()
    assert (m.mass_weighted_age_50 <= m.mass_weighted_age_84).all()
    with np.load(cp) as pack:
        n = int(pack[pack.files[0]][0])
        precision = np.zeros((n, n))
        precision[np.triu_indices(n)] = pack[pack.files[1]]
        low = np.tril_indices(n, -1)
        precision[low] = precision.T[low]
    assert n == len(hd)
    ids = m.source_row.to_numpy()
    selector = np.eye(n)[:, ids]
    sol = cho_solve(cho_factor(precision, lower=True), selector)
    c = sol[ids, :]
    closure = float(np.max(abs(precision @ sol - selector)))
    assert closure < 1e-10
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    mu = (
        5
        * np.log10(
            cosmo.comoving_distance(m.zHD.to_numpy()).value * (1 + m.zHEL.to_numpy())
        )
        + 25
    )
    y = m.MU.to_numpy() - mu
    age = m.mass_weighted_age_50.to_numpy()
    av = m["dust:Av_50"].to_numpy()
    mass = m.stellar_mass_50.to_numpy()
    designs = {
        "age": np.column_stack([np.ones(len(m)), age - age.mean()]),
        "age_AV_mass": np.column_stack(
            [np.ones(len(m)), age - age.mean(), av - av.mean(), mass - mass.mean()]
        ),
    }
    fits = {}
    l = np.linalg.cholesky(c)
    for name, x in designs.items():
        b, bc, chi = gls(x, y, c)
        wx = solve_triangular(l, x, lower=True)
        wy = solve_triangular(l, y, lower=True)
        independent = np.linalg.lstsq(wx, wy, rcond=None)[0]
        check = float(np.max(abs(independent - b)))
        assert check < 1e-10
        se = np.sqrt(bc[1, 1])
        delta = 0.03 / se
        fits[name] = {
            "coefficients": b.tolist(),
            "coefficient_se": np.sqrt(np.diag(bc)).tolist(),
            "age_slope_mag_per_Gyr": float(b[1]),
            "nominal_fixed_age_95": [float(b[1] - 1.96 * se), float(b[1] + 1.96 * se)],
            "chi2": chi,
            "df": len(m) - x.shape[1],
            "independent_whitened_lstsq_max_coefficient_difference": check,
            "fixed_age_detection_power_for_abs_003": float(
                norm.cdf(-1.96 - delta) + norm.sf(1.96 - delta)
            ),
            "fixed_age_80pct_detection_slope": float(
                (norm.ppf(0.975) + norm.ppf(0.8)) * se
            ),
        }
    result = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "code_sha256": sha(__file__),
        "helper_sha256": sha(HERE / "analyse.py"),
        "reader_sha256": sha(ROOT / "lib/records.py"),
        "scope": "New descriptive, non-preregistered low-z overlap check. No age effect selected using results; both declared simple and host-controlled fits reported.",
        "inputs": [
            {"path": str(p), "sha256": sha(p)}
            for p in [
                hp,
                cp,
                tp,
                folder / "README.md",
                folder / "DES-Dovekie-SN_Likelihood.py",
            ]
        ],
        "normalized_name_matches": nmatch,
        "selected_valid_and_host_dDLR_lt4": len(m),
        "identities": m[["CID", "transient", "IDSURVEY", "zHD", "d_dlr"]].to_dict(
            "records"
        ),
        "zHD_range": [float(m.zHD.min()), float(m.zHD.max())],
        "covariance": "Unpacked full released total precision in Hubble-diagram order, solved for covariance columns, then subset; metadata order not used; no statistical covariance added again.",
        "released_distance_semantics": "MU is released width/colour/mass-step and bias-corrected distance; MUERR renormalizes BEAMS probability. STAT+SYS includes the release's statistical and systematic covariance. Existing corrections are accepted rather than refitted or added again; their underlying physical accuracy is not established by this test.",
        "precision_inverse_subset_max_closure_error": closure,
        "covariance_minimum_eigenvalue": float(np.linalg.eigvalsh(c).min()),
        "fits": fits,
        "limitations": [
            "TITAN host assignment lacks coordinates in this snapshot; exact unique event names and d_DLR used, not independently verified host positions.",
            "Global posterior-median ages and dust/mass are treated as fixed; no joint host likelihood or covariance is available.",
            "Selected nearby Foundation overlap cannot measure high-redshift transport or represent DES main survey.",
            "Reference LCDM only supplies low-z distance shape; free intercept removes absolute scale.",
            "No claim that a weak slope proves standard correction adequate.",
        ],
    }
    out = ROOT / "studies/host_ages/results/environment_validation/dovekie-summary.json"
    out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
