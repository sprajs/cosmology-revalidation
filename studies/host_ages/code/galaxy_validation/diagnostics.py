#!/usr/bin/env python3
"""Independent observing-modality checks of SN-host SED age summaries.

No spectral index is treated as a model-free age; all brightness fits are
conditional predictive diagnostics, not causal or cosmological corrections.
"""

from pathlib import Path
import argparse, datetime, hashlib, json, sys
from functools import lru_cache
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr, ks_2samp
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import extinction
from acquire import ROOT, WORK, OUT, sha
from query_sdss import ID_DTYPES
from match import ARCHIVE

ENV = ROOT / "studies/host_ages/code/environment_validation"
sys.path.insert(0, str(ENV))
import ztf, titan

HERE = Path(__file__).parent
DESIGN = json.loads((HERE / "analysis-design.json").read_text())
SEED = DESIGN["seed"]
KEY = ["plateid", "mjd", "fiberid"]
COSMO = FlatLambdaCDM(H0=70, Om0=0.3)


def quantile(x):
    x = np.asarray(x)
    x = x[np.isfinite(x)]
    return (
        {
            "n": len(x),
            "median": float(np.median(x)),
            "q16_q84": np.quantile(x, [0.16, 0.84]).tolist(),
            "min_max": [float(x.min()), float(x.max())],
        }
        if len(x)
        else {"n": 0}
    )


def association(frame, field, controls=(), draws=2000):
    cols = ["sed_age", field] + list(controls)
    d = frame[cols].replace([np.inf, -np.inf], np.nan).dropna().to_numpy()
    if len(d) < 10:
        return {"n": len(d), "status": "too few observations"}

    def statistic(a):
        rr = np.column_stack([rankdata(a[:, i]) for i in range(a.shape[1])])
        if controls:
            X = np.column_stack([np.ones(len(rr)), rr[:, 2:]])
            rr[:, :2] -= X @ np.linalg.lstsq(X, rr[:, :2], rcond=None)[0]
        return float(np.corrcoef(rr[:, :2].T)[0, 1])

    point = statistic(d)
    rng = np.random.default_rng(SEED)
    boot = [statistic(d[rng.integers(len(d), size=len(d))]) for _ in range(draws)]
    return {
        "n": len(d),
        "spearman" if not controls else "partial_spearman": point,
        "physical_host_bootstrap_95": np.quantile(boot, [0.025, 0.975]).tolist(),
        "rank_controls": list(controls),
    }


@lru_cache(maxsize=1)
def load_brightness():
    base = ARCHIVE / "sources/updates/2026-09-20-ztf"
    paths = [
        base / "extracted/ztfsniadr2_lite/tables" / k
        for k in ["snia_data.csv", "globalhost_data.csv", "localhost_data.csv"]
    ] + [base / "originals/Ginolin25ab_masterlist--903d65d.csv"]
    sn, glo, loc, master = [pd.read_csv(p) for p in paths]
    d, _, _ = ztf.sample(sn, glo, loc, master)
    d = ztf.prepare(d)
    t = pd.read_csv(WORK / "titan-host-properties.csv")
    q = d.merge(
        t,
        left_on="iau_name",
        right_on="transient",
        suffixes=("", "_titan"),
        validate="one_to_one",
    )
    cut = (q.d_dlr < 4) & (q.d_dlr_titan < 4)
    cols = [
        f + "_" + suffix for f in titan.FIELDS.values() for suffix in ["16", "50", "84"]
    ]
    cut &= np.isfinite(q[cols]).all(axis=1)
    for f in titan.FIELDS.values():
        cut &= (q[f + "_84"] >= q[f + "_50"]) & (q[f + "_50"] >= q[f + "_16"])
    q = q[cut].copy().reset_index(drop=True)
    return q, paths


def ledger():
    h = pd.read_csv(WORK / "mpajhu-selected-hosts.csv").rename(
        columns={"mpa_PLATEID": "plateid", "mpa_MJD": "mjd", "mpa_FIBERID": "fiberid"}
    )
    s = pd.read_csv(WORK / "sdss-spectral-measurements.csv", dtype=ID_DTYPES)
    assert s.specObjID.is_unique and not s.duplicated(KEY).any()
    q = h.merge(s, on=KEY, how="left", validate="many_to_one", indicator="sdss_match")
    r = pd.read_csv(WORK / "remeasured-dn4000.csv")
    q = q.merge(r, on=KEY, how="left", validate="many_to_one")
    t, _ = load_brightness()
    q["in_401_TITAN_ZTF_brightness_sample"] = q.source_id.isin(t.ztfname)
    q["index_quality"] = (
        (q.reliable == 1)
        & (q.d4000_n > 0)
        & (q.d4000_n_err > 0)
        & (q.lick_hd_a_sub_err > 0)
        & np.isfinite(
            q[["d4000_n", "d4000_n_err", "lick_hd_a_sub", "lick_hd_a_sub_err"]]
        ).all(axis=1)
    )
    q["physical_host"] = q.mpa_PHOTOID.astype(str)
    # Legacy SDSS fixed aperture and actual SN separation are separate quantities.
    valid = np.isfinite(q[["sn_ra", "sn_dec", "spectral_ra", "spectral_dec"]]).all(
        axis=1
    )
    q["SN_to_fibre_arcsec"] = np.nan
    q.loc[valid, "SN_to_fibre_arcsec"] = (
        SkyCoord(
            q.loc[valid, "sn_ra"].to_numpy() * u.deg,
            q.loc[valid, "sn_dec"].to_numpy() * u.deg,
        )
        .separation(
            SkyCoord(
                q.loc[valid, "spectral_ra"].to_numpy() * u.deg,
                q.loc[valid, "spectral_dec"].to_numpy() * u.deg,
            )
        )
        .arcsec
    )
    q["fibre_diameter_kpc"] = 3 * COSMO.kpc_proper_per_arcmin(q.z.to_numpy()).value / 60
    q["Halpha_Hbeta"] = q.h_alpha_flux / q.h_beta_flux
    q["logNII_Halpha"] = np.log10(
        q.nii_6584_flux.where(q.nii_6584_flux > 0)
        / q.h_alpha_flux.where(q.h_alpha_flux > 0)
    )
    q["logOIII_Hbeta"] = np.log10(
        q.oiii_5007_flux.where(q.oiii_5007_flux > 0)
        / q.h_beta_flux.where(q.h_beta_flux > 0)
    )
    q["BPT_SF"] = True
    for name in ["h_alpha", "h_beta", "oiii_5007", "nii_6584"]:
        q["BPT_SF"] &= (q[name + "_flux_err"] > 0) & (
            q[name + "_flux"] / q[name + "_flux_err"] > 3
        )
    q["BPT_SF"] &= (q.logNII_Halpha < 0.05) & (
        q.logOIII_Hbeta < 0.61 / (q.logNII_Halpha - 0.05) + 1.3
    )
    k = extinction.fitzpatrick99(np.array([4861.0, 6563.0]), 3.1, 3.1)
    q["Balmer_EBV_F99Rv31"] = (
        2.5 / (k[0] - k[1]) * np.log10(q.Halpha_Hbeta.where(q.Halpha_Hbeta > 0) / 2.86)
    )
    q["Balmer_AV_F99Rv31"] = 3.1 * q.Balmer_EBV_F99Rv31
    q.to_csv(WORK / "host-spectrum-age-ledger.csv", index=False)
    # Finite catalogue photometric IDs must survive CSV round trip as exact digits.
    saved = pd.read_csv(WORK / "host-spectrum-age-ledger.csv", dtype=ID_DTYPES)
    for k in ID_DTYPES:
        assert saved[k].fillna("").equals(q[k].astype("string").fillna(""))
    return q


def spectral_checks(q):
    u = q[q.sdss_match == "both"].drop_duplicates("specObjID")
    v = u[u.valid.eq(True) & (u.d4000_n > 0) & (u.d4000_n_err > 0)]
    result = {
        "unique_spectra": len(u),
        "valid_reintegrated_catalogue_pairs": len(v),
        "raw_Dn_minus_MPA_raw": quantile(v.raw_Dn4000 - v.d4000_n),
        "MW_Dn_minus_MPA_raw": quantile(v.MW_corrected_Dn4000 - v.d4000_n),
        "MW_Dn_minus_MPA_emission_subtracted": quantile(
            v.MW_corrected_Dn4000 - v.d4000_n_sub
        ),
        "redshift_abs_difference_raw_vs_SQL": quantile(abs(v.spec_z - v.spectral_z)),
        "MW_delta_Dn": quantile(v.MW_corrected_Dn4000 - v.raw_Dn4000),
    }
    v = v.copy()
    v["delta"] = v.MW_corrected_Dn4000 - v.d4000_n
    result["largest_absolute_raw_index_discrepancies"] = [
        {
            "specObjID": r.specObjID,
            "difference": r.delta,
            "catalogue": r.d4000_n,
            "remeasured": r.MW_corrected_Dn4000,
            "SNR": r.sn_median,
        }
        for r in v.sort_values("delta", key=abs).tail(8).itertuples()
    ]
    return result


def observed_relations(q):
    out = {}
    for sample in ["ZTF", "G11", "R19"]:
        d = (
            q[(q["sample"] == sample) & q.index_quality & np.isfinite(q.sed_age)]
            .drop_duplicates("physical_host")
            .copy()
        )
        # Invalid absent catalogue masses are not physical control measurements.
        d.loc[d.sed_mass <= 0, "sed_mass"] = np.nan
        out[sample] = {
            "n_distinct_hosts": len(d),
            "ages": quantile(d.sed_age),
            "central_fibre_diameter_kpc": quantile(d.fibre_diameter_kpc),
            "SN_outside_fibre": int((d.SN_to_fibre_arcsec > 1.5).sum()),
            "SN_position_available": int(d.SN_to_fibre_arcsec.notna().sum()),
            "associations": {},
        }
        for idx in ["d4000_n", "lick_hd_a_sub", "MW_corrected_Dn4000"]:
            out[sample]["associations"][idx] = {
                "raw": association(d, idx),
                "mass_redshift_conditioned": association(d, idx, ["sed_mass", "z"]),
            }
        sf = d[d.BPT_SF]
        if len(sf) > 5:
            x = sf[["sed_Av", "Balmer_AV_F99Rv31"]].dropna()
            out[sample]["starforming_balmer"] = {
                "n": len(sf),
                "decrement": quantile(sf.Halpha_Hbeta),
                "conditional_AV": quantile(sf.Balmer_AV_F99Rv31),
                "SED_AV_vs_Balmer_AV_spearman": (
                    float(spearmanr(x.sed_Av, x.Balmer_AV_F99Rv31).statistic)
                    if len(x) > 5
                    else None
                ),
                "paired_n": len(x),
                "Balmer_minus_SED_AV": quantile(x.Balmer_AV_F99Rv31 - x.sed_Av),
                "negative_unclipped_Balmer_AV_count": int(
                    (sf.Balmer_AV_F99Rv31 < 0).sum()
                ),
            }
    return out


def manga():
    a = pd.read_csv(WORK / "firefly-miles-selected-hosts.csv")
    b = pd.read_csv(WORK / "firefly-mastar-selected-hosts.csv")
    d = a.merge(
        b[
            ["source_id", "ff_MANGAID", "ff_PLATEIFU"]
            + [k for k in b if k.startswith("ff_MW_AGE") or k.startswith("ff_LW_AGE")]
        ],
        on=["source_id", "ff_MANGAID", "ff_PLATEIFU"],
        suffixes=("_miles", "_mastar"),
        validate="one_to_one",
    )
    d.to_csv(WORK / "firefly-paired-models.csv", index=False)
    out = {
        "matched_host_rows": len(d),
        "unique_MaNGA_galaxies": d.ff_MANGAID.nunique(),
        "samples": {},
    }
    for sm in ["ZTF", "G11", "R19", "all"]:
        q = (
            (d if sm == "all" else d[d["sample"] == sm])
            .drop_duplicates("ff_MANGAID")
            .copy()
        )
        r = {"n": len(q)}
        for weighting in ["MW", "LW"]:
            for loc in ["3ARCSEC", "1Re"]:
                key = f"ff_{weighting}_AGE_{loc}"
                good = (
                    np.isfinite(q[key + "_miles"])
                    & np.isfinite(q[key + "_mastar"])
                    & (abs(q[key + "_miles"]) < 5)
                    & (abs(q[key + "_mastar"]) < 5)
                )
                v = q[good]
                r[f"{weighting}_{loc}_MaStar_minus_MILES_log10Gyr"] = quantile(
                    v[key + "_mastar"] - v[key + "_miles"]
                )
        for lib in ["miles", "mastar"]:
            col = "ff_MW_AGE_3ARCSEC_" + lib
            q["spectral_mass_weighted_age_Gyr"] = 10 ** q[col].where(abs(q[col]) < 5)
            r[lib + "_central_mass_weighted_age_Gyr"] = quantile(
                q.spectral_mass_weighted_age_Gyr
            )
            r[lib + "_central_vs_SED_age"] = association(
                q, "spectral_mass_weighted_age_Gyr"
            )
            r[lib + "_central_minus_SED_Gyr"] = quantile(
                q.spectral_mass_weighted_age_Gyr - q.sed_age
            )
            r[lib + "_central_minus_1Re_log10Gyr"] = quantile(
                q[col].where(abs(q[col]) < 5)
                - q["ff_MW_AGE_1Re_" + lib].where(abs(q["ff_MW_AGE_1Re_" + lib]) < 5)
            )
            r[lib + "_central_older_than_fiducial_universe"] = int(
                (
                    q.spectral_mass_weighted_age_Gyr > COSMO.age(q.z.to_numpy()).value
                ).sum()
            )
        out["samples"][sm] = r
    return out


def matrix(d, names):
    return np.column_stack(
        [
            (
                d["d4000_n"].to_numpy() - 1.5
                if k == "Dn4000"
                else (
                    d.lick_hd_a_sub.to_numpy() - 2
                    if k == "Hdelta"
                    else titan.matrix(d, [k])[:, 0]
                )
            )
            for k in names
        ]
    )


def brightness(q):
    d = load_brightness()[0].copy()
    s = q[q["sample"].eq("ZTF") & q.index_quality]
    d = d.merge(
        s[
            [
                "source_id",
                "physical_host",
                "specObjID",
                "d4000_n",
                "d4000_n_err",
                "lick_hd_a_sub",
                "lick_hd_a_sub_err",
            ]
        ],
        left_on="ztfname",
        right_on="source_id",
        validate="one_to_one",
    )
    d["fold"] = [
        int(hashlib.sha256((str(SEED) + v).encode()).hexdigest()[:8], 16) % 5
        for v in d.physical_host
    ]
    d.drop(columns=["lc_cov"]).to_csv(WORK / "brightness-cohort.csv", index=False)
    models = {
        "SALT": titan.BASE,
        "host": titan.CONTROL,
        "host_age": titan.CONTROL + ["host_age"],
        "host_spectra": titan.CONTROL + ["Dn4000", "Hdelta"],
        "host_spectra_age": titan.CONTROL + ["Dn4000", "Hdelta", "host_age"],
    }

    def var(dd, bb, nn, errors=False):
        v = titan.variance(dd, bb, nn, False)
        if errors:
            for name, key in [
                ("Dn4000", "d4000_n_err"),
                ("Hdelta", "lick_hd_a_sub_err"),
            ]:
                if name in nn:
                    v += bb[nn.index(name)] ** 2 * dd[key].to_numpy() ** 2
        return v

    def run(errors=False, cv=True):
        res = {}
        scores = {}
        preds = {}
        for label, names in models.items():
            vf = lambda dd, bb, nn: var(dd, bb, nn, errors)
            f, b, c, sig = ztf.fit(d, names, matrix, vf)
            res[label] = f
            if cv:
                pred = np.zeros(len(d))
                pv = np.zeros(len(d))
                for fold in range(5):
                    tr = d.fold != fold
                    te = ~tr
                    ff, bb, cc, ss = ztf.fit(d[tr], names, matrix, vf)
                    xx = matrix(d[te], names)
                    pred[te] = xx @ bb
                    pv[te] = (
                        vf(d[te], bb, names)
                        + ss**2
                        + np.einsum("ni,ij,nj->n", xx, cc, xx)
                    )
                score = -0.5 * (
                    np.log(2 * np.pi * pv) + (d.y.to_numpy() - pred) ** 2 / pv
                )
                scores[label] = score
                preds[label] = pred
                f["heldout_log_predictive_score"] = float(score.sum())
                f["heldout_RMSE_mag"] = float(np.sqrt(np.mean((d.y - pred) ** 2)))
        pairs = {}
        rng = np.random.default_rng(SEED)
        if cv:
            for added, base in [
                ("host", "SALT"),
                ("host_age", "host"),
                ("host_spectra", "host"),
                ("host_spectra_age", "host_spectra"),
            ]:
                delta = scores[added] - scores[base]
                units = np.array(
                    [
                        delta[d.physical_host == g].sum()
                        for g in sorted(set(d.physical_host))
                    ]
                )
                boot = rng.choice(units, (2000, len(units)), replace=True).sum(axis=1)
                pairs[added + "_minus_" + base] = {
                    "delta_log_predictive_score": float(delta.sum()),
                    "host_bootstrap_95_conditional_on_fitted_folds": np.quantile(
                        boot, [0.025, 0.975]
                    ).tolist(),
                }
            pd.DataFrame(
                dict(
                    source_id=d.source_id,
                    physical_host=d.physical_host,
                    fold=d.fold,
                    **{k + "_score": v for k, v in scores.items()},
                    **{k + "_prediction": v for k, v in preds.items()},
                )
            ).to_csv(WORK / "brightness-heldout.csv", index=False)
        return {"fits": res, "predictive_comparisons": pairs}

    return {
        "n": len(d),
        "distinct_hosts": d.physical_host.nunique(),
        "fold_counts": d.groupby("fold").size().to_dict(),
        "models": models,
        "primary_fixed_indices_and_SED_medians": run(),
        "independent_index_error_propagation_sensitivity": run(True, False),
    }


def selection(q):
    allhosts = pd.read_csv(WORK / "host-coordinates.csv")
    allhosts = allhosts[
        (allhosts["sample"] == "ZTF")
        & allhosts.host_accepted
        & allhosts.valid_position_redshift
        & np.isfinite(allhosts.sed_age)
    ]
    matched = allhosts[allhosts.source_id.isin(q.loc[q.index_quality, "source_id"])]
    missing = allhosts[~allhosts.source_id.isin(matched.source_id)]
    answer = {
        "status": "Descriptive selection audit after primary association outcomes, not a new discovery test",
        "all_accepted_TITAN_ZTF_SN_rows": len(allhosts),
        "with_good_spectra_SN_rows": len(matched),
        "without_good_spectra_SN_rows": len(missing),
        "comparison": {},
    }
    for key in ["sed_age", "sed_mass", "sed_Av", "z"]:
        x = matched[key].dropna()
        y = missing[key].dropna()
        answer["comparison"][key] = {
            "with_spectra": quantile(x),
            "without_spectra": quantile(y),
            "KS_distance": float(ks_2samp(x, y).statistic),
        }
    return answer


def main():
    global ARCHIVE
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, default=ARCHIVE)
    args = parser.parse_args()
    ARCHIVE = args.archive
    q = ledger()
    print("Ledger ready", len(q), flush=True)
    out = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "crossmatch_closure": {
            "accepted_host_rows": len(q),
            "missing_SQL_rows": q.loc[
                q.sdss_match != "both", ["source_id"] + KEY
            ].to_dict("records"),
            "returned_host_rows": int((q.sdss_match == "both").sum()),
            "unique_returned_spectra": int(q.specObjID.nunique()),
            "TITAN_ZTF_401_overlap": int(q.in_401_TITAN_ZTF_brightness_sample.sum()),
            "all_identifier_fields_read_as_exact_strings": list(ID_DTYPES),
        },
        "spectroscopic_selection": selection(q),
        "spectral_remeasurement": spectral_checks(q),
        "age_index_relations": observed_relations(q),
        "MaNGA_FIREFLY": manga(),
        "brightness": brightness(q),
    }
    out["provenance"] = {
        "code_sha256": {p.name: sha(p) for p in HERE.glob("*.py")},
        "design_sha256": {p.name: sha(p) for p in HERE.glob("*design.json")},
        "dependencies_sha256": {
            p.name: sha(p)
            for p in [ENV / "ztf.py", ENV / "titan.py", ENV / "analyse.py"]
        },
        "input_sha256": {
            str(p.relative_to(ROOT)): sha(p)
            for p in [
                WORK / "mpajhu-selected-hosts.csv",
                WORK / "sdss-spectral-measurements.csv",
                WORK / "remeasured-dn4000.csv",
                WORK / "firefly-miles-selected-hosts.csv",
                WORK / "firefly-mastar-selected-hosts.csv",
                WORK / "titan-host-properties.csv",
            ]
        },
        "inherited_input_sha256": {
            str(p.relative_to(ARCHIVE)): sha(p) for p in load_brightness()[1]
        },
        "output_sha256": {
            str(p.relative_to(ROOT)): sha(p)
            for p in [
                WORK / "host-spectrum-age-ledger.csv",
                WORK / "brightness-cohort.csv",
                WORK / "brightness-heldout.csv",
                WORK / "firefly-paired-models.csv",
            ]
        },
    }
    out["limits"] = [
        "Host positions and redshifts verify the SDSS counterpart, conditional on the ZTF/TITAN host association; the TITAN file itself has no independent host coordinates.",
        "Spectroscopy is an independent observing modality but shares stellar population, dust, IMF and library assumptions; MPA line deblending/emission subtraction is model dependent.",
        "Central SDSS fibres and FIREFLY apertures differ from global SED ages and local SN environments; a difference is not automatically an age calibration error.",
        "SED posterior medians are not age likelihoods, and nuisance controls can mediate a physical age effect. No null conditional coefficient establishes absence of age physics.",
        "Index error propagation ignores shared spectral-error covariance; raw spectrum errors ignore supplied resampling correlations.",
        "Brightness sample is low-redshift selected ZTF photometry with no survey-selection retraining, global calibration covariance or cosmology recovery.",
        "Five-fold prediction intervals condition on overlapping fitted training sets and do not include complete training or survey uncertainty.",
    ]
    (OUT / "diagnostics-summary.json").write_text(
        json.dumps(out, indent=2, allow_nan=False) + "\n"
    )
    print(
        json.dumps(
            {
                "output": "diagnostics-summary.json",
                "brightness_n": out["brightness"]["n"],
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
