#!/usr/bin/env python3
"""Independent read-only audit of signed photon-level survey age interventions."""

from pathlib import Path
import datetime, hashlib, io, json, sys
import numpy as np, pandas as pd
from astropy.io import fits
from extinction import fitzpatrick99
from acquire import ROOT, WORK, OUT, sha

sys.path.insert(0, str(ROOT))
from lib.records import fitres

SURVEY = ROOT / ".work/survey-physics"


def dump(name):
    p = SURVEY / "simulations" / name / (name + ".DUMP")
    lines = p.read_text().splitlines()
    header = next(x.split()[1:] for x in lines if x.startswith("VARNAMES:"))
    return (
        pd.read_csv(
            io.StringIO("\n".join(x[3:] for x in lines if x.startswith("SN:"))),
            sep=r"\s+",
            names=header,
        ),
        p,
    )


def raw(name):
    p = SURVEY / "simulations" / name / name
    with fits.open(str(p) + "_HEAD.FITS") as h:
        head = h[1].data.copy()
    with fits.open(str(p) + "_PHOT.FITS") as h:
        phot = h[1].data.copy()
    mapping = {(str(row["SNID"]).strip(), int(row["SIM_LIBID"])): row for row in head}
    assert len(mapping) == len(head)
    return (
        mapping,
        phot,
        [Path(str(p) + "_" + suffix + ".FITS") for suffix in ["HEAD", "PHOT"]],
    )


def check(prefix):
    nname = prefix + "_EVAL_NOMINAL"
    aname = prefix + "_EVAL_AGE"
    n, npth = dump(nname)
    a, apth = dump(aname)
    assert (
        not n.duplicated(["CID", "LIBID"]).any()
        and not a.duplicated(["CID", "LIBID"]).any()
    )
    matched = n.merge(
        a, on=["CID", "LIBID"], suffixes=("_nominal", "_age"), validate="one_to_one"
    )
    assert len(matched) == len(n) == len(a)
    latents = [
        "GENZ",
        "GALID",
        "SN_age",
        "SALT2x1",
        "SALT2c",
        "AV",
        "RV",
        "MU",
        "PEAKMJD",
    ]
    latent_valid_n = n.SN_age.between(0, 14) & (n.RV > 0) & (n.AV >= 0)
    latent_valid_a = a.SN_age.between(0, 14) & (a.RV > 0) & (a.AV >= 0)
    assert not ((~latent_valid_n) & n.FLAG_ACCEPT.astype(bool)).any()
    assert not ((~latent_valid_a) & a.FLAG_ACCEPT.astype(bool)).any()
    equality = {
        k: bool((matched[k + "_nominal"] == matched[k + "_age"]).all()) for k in latents
    }
    assert all(equality.values())
    nh, nphot, npaths = raw(nname)
    ah, aphot, apaths = raw(aname)
    shared = sorted(nh.keys() & ah.keys())
    epochs = 0
    magerrors = []
    raterrors = []
    gamma = []
    mBdiff = []
    mu = []
    noise = []
    flags = 0
    for key in shared:
        x, y = nh[key], ah[key]
        age = float(x["SIM_HOSTLIB(SN_age)"])
        shift = -0.03 * (age - 3)
        assert float(y["SIM_HOSTLIB(SN_age)"]) == age
        px = nphot[int(x["PTROBS_MIN"]) - 1 : int(x["PTROBS_MAX"])]
        py = aphot[int(y["PTROBS_MIN"]) - 1 : int(y["PTROBS_MAX"])]
        assert len(px) == int(x["NOBS"]) and len(py) == int(y["NOBS"])
        # Match actual exposures; do not assume array positions coincide after cuts.
        xx = pd.DataFrame(
            {
                k: np.array(px[k]).astype(str if k == "BAND" else float)
                for k in [
                    "MJD",
                    "BAND",
                    "SIM_MAGOBS",
                    "FLUXCAL",
                    "FLUXCALERR",
                    "PHOTFLAG",
                ]
            }
        )
        yy = pd.DataFrame(
            {
                k: np.array(py[k]).astype(str if k == "BAND" else float)
                for k in [
                    "MJD",
                    "BAND",
                    "SIM_MAGOBS",
                    "FLUXCAL",
                    "FLUXCALERR",
                    "PHOTFLAG",
                ]
            }
        )
        # Native float32 MJD can repeat within a band. Verify the complete
        # cadence sequence, then use occurrence number to retain every exposure.
        assert len(px) == len(py)
        for field in ["MJD", "BAND", "ZEROPT", "SKY_SIG"]:
            assert np.array_equal(px[field], py[field]), (key, field)
        xx["occurrence"] = xx.groupby(["MJD", "BAND"]).cumcount()
        yy["occurrence"] = yy.groupby(["MJD", "BAND"]).cumcount()
        paired = xx.merge(
            yy,
            on=["MJD", "BAND", "occurrence"],
            suffixes=("_nominal", "_age"),
            validate="one_to_one",
        )
        good = paired.SIM_MAGOBS_nominal.between(
            -30, 80
        ) & paired.SIM_MAGOBS_age.between(-30, 80)
        paired = paired[good]
        delta = paired.SIM_MAGOBS_age.to_numpy() - paired.SIM_MAGOBS_nominal.to_numpy()
        epochs += len(paired)
        magerrors.extend((delta - shift).tolist())
        raterrors.extend((10 ** (-0.4 * delta) / 10 ** (-0.4 * shift) - 1).tolist())
        gamma.append(float(y["SIM_SALT2gammaDM"] - x["SIM_SALT2gammaDM"]) - shift)
        mBdiff.append(float(y["SIM_SALT2mB"] - x["SIM_SALT2mB"]))
        mu.append(float(y["SIM_DLMU"] - x["SIM_DLMU"]))
        flags += int((paired.PHOTFLAG_age != paired.PHOTFLAG_nominal).sum())
    # Float32 SIM_MAGOBS controls this ~few micro-magnitude tolerance.
    assert max(abs(np.array(magerrors))) < 1e-5 and max(abs(np.array(gamma))) < 3e-8
    assert max(abs(np.array(mu))) == 0
    # Transform released native x0 covariance independently to magnitude covariance.
    fp = SURVEY / "fits" / nname / "fit.FITRES.TEXT"
    f = fitres(fp)
    dm = -2.5 / np.log(10) / f.x0.to_numpy()
    c = np.zeros((len(f), 3, 3))
    c[:, 0, 0] = (dm * f.x0ERR) ** 2
    c[:, 1, 1] = f.x1ERR**2
    c[:, 2, 2] = f.cERR**2
    c[:, 0, 1] = c[:, 1, 0] = dm * f.COV_x1_x0
    c[:, 0, 2] = c[:, 2, 0] = dm * f.COV_c_x0
    c[:, 1, 2] = c[:, 2, 1] = f.COV_x1_c
    j = np.array([1, 0.15, -3.14])
    v1 = np.einsum("i,nij,j->n", j, c, j)
    cv = c.copy()
    cv[:, 0, 0] = f.x0ERR**2
    cv[:, 0, 1] = cv[:, 1, 0] = f.COV_x1_x0
    cv[:, 0, 2] = cv[:, 2, 0] = f.COV_c_x0
    J = np.column_stack([dm, np.full(len(f), 0.15), np.full(len(f), -3.14)])
    v2 = np.einsum("ni,nij,nj->n", J, cv, J)
    result = {
        "generated_attempt_pairs": len(matched),
        "populated_SN_dust_latents_nominal": int(latent_valid_n.sum()),
        "unpopulated_preselection_latents_nominal": int((~latent_valid_n).sum()),
        "latent_equality": equality,
        "accepted_both_photometric_objects": len(shared),
        "compared_noiseless_epochs": epochs,
        "injection": "delta_mag=-0.030*(true_SN_age_Gyr-3), before detection; older simulated progenitors brighter",
        "maximum_noiseless_magnitude_identity_error": float(
            max(abs(np.array(magerrors)))
        ),
        "maximum_noiseless_fluxratio_relative_error": float(
            max(abs(np.array(raterrors)))
        ),
        "maximum_gammaDM_identity_error": float(max(abs(np.array(gamma)))),
        "maximum_SIM_mB_change": float(max(abs(np.array(mBdiff)))),
        "maximum_true_distance_modulus_change": float(max(abs(np.array(mu)))),
        "semantics": "SIM_mB/SALT2mB amplitude excludes separately applied gammaDM shift; its unchanged value does not indicate failed injection. SIM_MAGOBS confirms the photons carry the signed shift.",
        "accepted_nominal": int(n.FLAG_ACCEPT.sum()),
        "accepted_age": int(a.FLAG_ACCEPT.sum()),
        "attempt_selection_disagreements": int(
            (matched.FLAG_ACCEPT_nominal != matched.FLAG_ACCEPT_age).sum()
        ),
        "epoch_PHOTFLAG_disagreements_commonobjects": flags,
        "variance_transform_max_difference": float(np.max(abs(v1 - v2))),
        "source_sha256": {
            str(p.relative_to(ROOT)): sha(p) for p in [npth, apth, fp] + npaths + apaths
        },
    }
    if prefix == "SPP":
        wave = np.linspace(1000.0, 30000.0, 2901)
        minimum = []
        populated = pd.concat([n[latent_valid_n], a[latent_valid_a]], ignore_index=True)
        for _, event in populated.iterrows():
            assert abs(event.RV - 3.1) < 1e-5 and event.AV >= 0
            minimum.append(
                float(np.min(fitzpatrick99(wave, float(event.AV), float(event.RV))))
            )
        result["physical_extinction_support"] = {
            "all_generated_RV_3p1": True,
            "all_generated_AV_nonnegative": True,
            "generated_evaluation_occurrences_tested": len(populated),
            "rejected_before_SN_dust_parameters_populated": int(
                (~latent_valid_n).sum() + (~latent_valid_a).sum()
            ),
            "rest_wavelength_grid_A": [1000, 30000],
            "grid_points": len(wave),
            "minimum_extinction_mag": min(minimum),
            "limit": "Validates only the imposed F99-screen support on this declared wavelength grid; does not establish its empirical host/SN population distribution or other generator physics.",
        }
    return result


def main():
    result = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "code_sha256": sha(Path(__file__)),
        "main_W22_campaign": check("SPH"),
        "fixed_RV31_campaign": check("SPP"),
        "limits": [
            "Checks the implemented mock intervention, not the reality of the assumed progenitor-age population.",
            "Main W22 population contains nonphysical low-Rv F99 screens; a successful signed-photon check does not validate that physical support.",
            "Pairs use CID+LIBID occurrence identities. Failed attempts reuse CID, so CID-only matching would be wrong.",
            "Pure-Ia classifier filtering cannot establish contaminant rejection, probability calibration or a full BEAMS likelihood.",
            "No observed-pixel noise equality is asserted; truth magnitude and exposure identity checks are distinct from fitting/selection response.",
        ],
    }
    (OUT / "survey-physics-review.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
