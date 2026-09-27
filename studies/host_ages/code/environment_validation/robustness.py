#!/usr/bin/env python3
"""Adversarial checks prompted by source/method review, not blind holdout tests."""
import datetime, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import titan, ztf
from analyse import ROOT, HERE, sha, gls


def main():
    work = ROOT / ".work/environment-validation"
    p = work / "titan-ztf-join.csv"
    d = pd.read_csv(p)
    d = d[d.selected_titan].reset_index(drop=True)
    d = ztf.prepare(d)
    result = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "scope": "Post-outcome sensitivity controls; no retrospective preregistration.",
        "code_sha256": sha(__file__),
        "input_join_sha256": sha(p),
        "dependencies_sha256": {
            f.name: sha(f)
            for f in [HERE / "analyse.py", HERE / "ztf.py", HERE / "titan.py"]
        },
    }
    cols = titan.CONTROL + ["host_age"]
    base, b, bc, sig = ztf.fit(
        d, cols, titan.matrix, lambda dd, bb, nn: titan.variance(dd, bb, nn, False)
    )
    q = d.copy()
    q["y"] += 0.2
    shift, _, _, _ = ztf.fit(
        q, cols, titan.matrix, lambda dd, bb, nn: titan.variance(dd, bb, nn, False)
    )
    result["uniform_magnitude_offset_test"] = {
        "offset_mag": 0.2,
        "age_slope_change": shift["parameters"]["host_age"]
        - base["parameters"]["host_age"],
        "intercept_change": shift["parameters"]["intercept"]
        - base["parameters"]["intercept"],
    }
    bins = np.floor((d.redshift.to_numpy() - 0.01) / 0.005).astype(int)
    d["zbin"] = bins
    present = sorted(set(bins))
    bin_names = ["zbin_" + str(i) for i in present[1:]]
    names = [k for k in cols if k != "z"] + bin_names

    def mx(dd, nn):
        return np.column_stack(
            [
                (
                    (dd.zbin.to_numpy() == int(k.split("_")[1])).astype(float)
                    if k.startswith("zbin_")
                    else titan.matrix(dd, [k])[:, 0]
                )
                for k in nn
            ]
        )

    f, _, _, _ = ztf.fit(
        d, names, mx, lambda dd, bb, nn: titan.variance(dd, bb, nn, False)
    )
    q = d.copy()
    q["y"] += 0.2 * np.sin(q.zbin.to_numpy())
    fs, _, _, _ = ztf.fit(
        q, names, mx, lambda dd, bb, nn: titan.variance(dd, bb, nn, False)
    )
    result["free_redshift_bin_intercepts"] = {
        "width": 0.005,
        "counts": d.zbin.value_counts().sort_index().to_dict(),
        "fit": f,
        "injected_bin_offset_max_abs_mag": float(
            abs(0.2 * np.sin(q.zbin.to_numpy())).max()
        ),
        "age_slope_change_under_bin_offsets": fs["parameters"]["host_age"]
        - f["parameters"]["host_age"],
        "interpretation": "The age slope uses within-bin variation; arbitrary bin-constant magnitude blinding is absorbed, but arbitrary intra-bin/object-specific perturbations are not.",
    }
    dm = d.stellar_mass_50 - d.mass_global
    dc = d.gz_colour_50 - d.restframe_gz_global
    result["host_catalogue_agreement"] = {
        "stellar_mass_median_TITAN_minus_ZTF_dex": float(np.median(dm)),
        "stellar_mass_nmad_dex": float(1.4826 * np.median(abs(dm - np.median(dm)))),
        "stellar_mass_absolute_difference_over_05dex": int((abs(dm) > 0.5).sum()),
        "restframe_gz_median_TITAN_minus_ZTF": float(np.median(dc)),
        "restframe_gz_nmad": float(1.4826 * np.median(abs(dc - np.median(dc)))),
        "mass_spearman": float(spearmanr(d.stellar_mass_50, d.mass_global).statistic),
        "note": "Different SED models/priors and host associations may produce differences; agreement is not a sky-coordinate validation.",
    }
    result["brightness_blinding_gate"] = {
        "release_table_source": "https://arxiv.org/html/2409.04346v2#S6.T3",
        "release_table_statement": "blinded flux zeropoint near 30",
        "companion_source": "https://inspirehep.net/files/91e12bac894170dcfd70a703f09aac5c",
        "companion_location": "A&A 697 A125 (2025), page 12, discussion of M",
        "companion_statement": "a blinding factor has been added to M for the ZTF sample",
        "status": "Common magnitude-intercept blinding is supported by the companion description and is harmless to fitted contrasts with a free intercept. Exact public-release transformation code was not recovered, so this is not a fully certified unblinded luminosity test. Free-redshift-bin checks protect against bin-constant magnitude shifts only.",
    }
    result["calibration_limits"] = {
        "catalogue": "DR2 known photometric nonlinearity; no full cross-object calibration or peculiar-velocity covariance supplied in three tables",
        "age_prior": "TITAN posterior-median summaries have substantial prior/model dependence; working-model calibration does not test physical age identification",
    }
    out = ROOT / "studies/host_ages/results/environment_validation/robustness.json"
    out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
