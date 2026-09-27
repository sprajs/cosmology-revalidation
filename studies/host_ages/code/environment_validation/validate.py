#!/usr/bin/env python3
"""Numerical/source validation distinct from scientific identification gates."""
import datetime, json
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.time import Time
from scipy.linalg import solve_triangular
from analyse import ROOT, HERE, sha, gls


def main():
    outdir = ROOT / "studies/host_ages/results/environment_validation"
    work = ROOT / ".work/environment-validation"
    names = [
        "dustpedia-summary.json",
        "ztf-summary.json",
        "titan-ztf-summary.json",
        "titan-calibration.json",
        "robustness.json",
    ]
    records = {n: json.loads((outdir / n).read_text()) for n in names}
    verified = {}
    for d in records.values():
        for item in d.get("inputs", []):
            path = Path(item["path"])
            assert path.exists() and sha(path) == item["sha256"], str(path)
            verified[str(path)] = item["sha256"]
    own = {
        "dustpedia-summary.json": "analyse.py",
        "ztf-summary.json": "ztf.py",
        "titan-ztf-summary.json": "titan.py",
        "titan-calibration.json": "calibrate.py",
        "robustness.json": "robustness.py",
    }
    for result, file in own.items():
        assert records[result]["code_sha256"] == sha(HERE / file), (
            result,
            "code changed since run",
        )
    for result in [
        "titan-ztf-summary.json",
        "titan-calibration.json",
        "robustness.json",
    ]:
        for file, digest in records[result]["dependencies_sha256"].items():
            assert sha(HERE / file) == digest, (result, file, "dependency changed")
    assert records["ztf-summary.json"]["shared_code_sha256"] == sha(HERE / "analyse.py")
    for result in ["dustpedia-summary.json", "ztf-summary.json"]:
        assert records[result]["design_sha256"] == sha(HERE / "design.json")
    assert records["titan-ztf-summary.json"]["design_sha256"] == sha(
        HERE / "titan-design.json"
    )
    assert records["titan-calibration.json"]["input_join_sha256"] == sha(
        work / "titan-ztf-join.csv"
    )
    p = pd.read_csv(ROOT / "data/distances/Pantheon+SH0ES.dat", sep=r"\s+")
    c = np.loadtxt(
        ROOT / "data/distances/Pantheon+SH0ES_STAT+SYS.cov", skiprows=1
    ).reshape(len(p), len(p))
    c = (c + c.T) / 2
    cal = pd.read_csv(work / "dustpedia-calibrators.csv")
    ids = cal.source_row.to_numpy()
    cc = c[np.ix_(ids, ids)]
    a = cal.AgeMW.to_numpy() / 1000
    x = np.column_stack([np.ones(len(cal)), a - a.mean()])
    y = (cal.m_b_corr - cal.CEPH_DIST).to_numpy()
    b, bc, _ = gls(x, y, cc)
    l = np.linalg.cholesky(cc)
    wx = solve_triangular(l, x, lower=True)
    wy = solve_triangular(l, y, lower=True)
    ind = np.linalg.lstsq(wx, wy, rcond=None)[0]
    err = float(np.max(abs(ind - b)))
    assert err < 1e-10
    cross = pd.read_csv(work / "dustpedia-pantheon-crosswalk.csv")
    alias = cross[cross.CID.isin(["2008fv_comb", "1994DRichmond"])]
    years = Time(alias.PKMJD.to_numpy(), format="mjd").to_datetime()
    assert all(int(str(cid)[:4]) == t.year for cid, t in zip(alias.CID, years))
    assert alias.coordinate_separation_arcsec.max() < 3
    r = records["robustness.json"]
    assert abs(r["uniform_magnitude_offset_test"]["age_slope_change"]) < 1e-7
    assert (
        abs(r["free_redshift_bin_intercepts"]["age_slope_change_under_bin_offsets"])
        < 1e-7
    )
    for model in records["titan-calibration.json"]["injections"].values():
        assert (
            model["converged"] == model["replicates_requested"]
            and not model["failures"]
        )
        assert abs(model["coverage95"] - 0.95) < 3 * np.sqrt(
            0.95 * 0.05 / model["converged"]
        )
    results = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "code_sha256": sha(__file__),
        "status": "pass for the declared numerical/source checks; physical age/dust identification remains unresolved",
        "input_files_verified": len(verified),
        "source_identities": verified,
        "result_sha256": {n: sha(outdir / n) for n in names},
        "whitened_lstsq_vs_full_covariance_GLS_max_coefficient_difference": err,
        "alias_epochs_and_coordinates_verified": alias[
            ["SN", "CID", "PKMJD", "coordinate_separation_arcsec"]
        ].to_dict("records"),
        "full_calibrator_covariance_minimum_eigenvalue": float(
            np.linalg.eigvalsh(cc).min()
        ),
        "checks": [
            "All recorded source, current code, shared dependency and design hashes match.",
            "Full covariance GLS agrees with an independent Cholesky-whitened least-squares solution.",
            "New catalogue aliases agree by sky and event year.",
            "Magnitude and redshift-bin nuisance-shift invariance pass.",
            "600 fitted-working-model simulations converge with nominal coverage consistent within Monte Carlo precision.",
        ],
        "not_validated": [
            "Raw photometric calibration and full peculiar-velocity covariance for ZTF.",
            "Exact public x0 blinding transformation code.",
            "Joint host age/dust/metallicity posterior, age-prior robustness and local-to-global mapping.",
            "Survey selection, population transport and residual cosmological bias.",
        ],
    }
    (outdir / "validation.json").write_text(
        json.dumps(results, indent=2, allow_nan=False) + "\n"
    )
    print(
        json.dumps(
            {k: v for k, v in results.items() if k not in ["source_identities"]},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
