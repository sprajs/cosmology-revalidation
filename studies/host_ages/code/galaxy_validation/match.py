#!/usr/bin/env python3
"""Build explicit SN-host/spectrum crosswalks before examining spectral outcomes."""

from pathlib import Path
import argparse, json, datetime, sys
import numpy as np, pandas as pd
from astropy.io import fits, ascii
from astropy.coordinates import SkyCoord, search_around_sky
from astropy import units as u
from acquire import ROOT, WORK, OUT, sha

ARCHIVE = (
    Path.home()
    / ".local/share/cosmology-revalidation/archive/17487bf659fcbdeeea072221492bac14b04a0a85"
)


def frame(data):
    return pd.DataFrame(
        {
            k: (
                np.array(data[k]).astype(str)
                if data[k].dtype.kind in "SU"
                else np.array(data[k]).astype(data[k].dtype.newbyteorder("="))
            )
            for k in data.names
            if data[k].ndim == 1
        }
    )


def age_table(p, lab):
    rows = []
    for line in p.read_text().splitlines()[1:]:
        c = [x.replace("\\", "").strip() for x in line.split("&")]
        if len(c) == 5:
            rows.append([str(int(c[0])), lab, float(c[1]), float(c[2])])
    return pd.DataFrame(rows, columns=["CID", "sample", "sed_age", "sed_age_error"])


def hosts(archive):
    base = archive / "sources/updates/2026-09-20-ztf/extracted/ztfsniadr2_lite/tables"
    paths = [
        base / "snia_data.csv",
        base / "globalhost_data.csv",
        WORK / "titan-host-properties.csv",
    ]
    sn, g, t = [pd.read_csv(p) for p in paths]
    d = sn.merge(g, on="ztfname", validate="one_to_one").merge(
        t,
        left_on="iau_name",
        right_on="transient",
        how="left",
        suffixes=("", "_titan"),
        validate="many_to_one",
    )
    z = pd.DataFrame(
        {
            "source_id": d.ztfname,
            "physical_SN": d.iau_name,
            "sample": "ZTF",
            "host_ra": d.ra_host,
            "host_dec": d.dec_host,
            "sn_ra": d.ra,
            "sn_dec": d.dec,
            "z": d.redshift,
            "sed_age": d.mass_weighted_age_50,
            "sed_age_error": (d.mass_weighted_age_84 - d.mass_weighted_age_16) / 2,
            "sed_mass": d.mass,
            "sed_Av": d["dust:Av_50"],
            "d_dlr": d.d_dlr,
            "d_dlr_titan": d.d_dlr_titan,
            "host_accepted": (d.d_dlr < 4)
            & ((d.d_dlr_titan < 4) | d.d_dlr_titan.isna()),
        }
    )
    allhosts = [z]
    gp = archive / "data/host_ages/gupta2011"
    paths += [
        gp / "table2.dat",
        gp / "ReadMe",
        archive / "data/host_ages/chung2025/table1.dat",
        archive / "data/host_ages/chung2025/table2.dat",
        ROOT / "data/distances/Pantheon+SH0ES.dat",
    ]
    old = ascii.read(
        gp / "table2.dat", format="cds", readme=str(gp / "ReadMe")
    ).to_pandas()
    old["CID"] = old.SNID.astype(str)
    aa = age_table(paths[-3], "G11")
    gg = aa.merge(old, on="CID", validate="one_to_one")
    pp = pd.read_csv(paths[-1], sep=r"\s+", dtype={"CID": str})
    gg = gg.merge(
        pp[pp.IDSURVEY == 1][["CID", "RA", "DEC"]],
        on="CID",
        how="left",
        validate="one_to_one",
    )
    allhosts.append(
        pd.DataFrame(
            {
                "source_id": "G11:" + gg.CID,
                "physical_SN": "SDSS:" + gg.CID,
                "sample": "G11",
                "host_ra": gg.RAdeg,
                "host_dec": gg.DEdeg,
                "sn_ra": gg.RA,
                "sn_dec": gg.DEC,
                "z": gg.z,
                "sed_age": gg.sed_age,
                "sed_age_error": gg.sed_age_error,
                "sed_mass": gg.M,
                "host_accepted": True,
            }
        )
    )
    pp = pd.read_csv(paths[-1], sep=r"\s+", dtype={"CID": str})
    r19 = age_table(paths[-2], "R19")
    p = WORK / "R19-campbell_global.tsv"
    paths.append(p)
    rg = pd.read_csv(p, sep="\t")
    rg["CID"] = rg.SNID.astype(str)
    rr = r19.merge(rg, on="CID", validate="one_to_one").merge(
        pp[pp.IDSURVEY == 1][["CID", "RA", "DEC", "HOST_LOGMASS"]],
        on="CID",
        how="left",
        validate="one_to_one",
    )
    allhosts.append(
        pd.DataFrame(
            {
                "source_id": "R19:" + rr.CID,
                "physical_SN": "SDSS:" + rr.CID,
                "sample": "R19",
                "host_ra": rr.ra,
                "host_dec": rr.dec,
                "sn_ra": rr.RA,
                "sn_dec": rr.DEC,
                "z": rr.redshift,
                "sed_age": rr.sed_age,
                "sed_age_error": rr.sed_age_error,
                "sed_mass": rr.HOST_LOGMASS,
                "host_accepted": True,
            }
        )
    )
    p = ROOT / ".work/environment-validation/dustpedia-pantheon-crosswalk.csv"
    if p.exists():
        paths.append(p)
        dd = pd.read_csv(p).drop_duplicates("CID")
    else:
        dd = pd.DataFrame(
            columns=[
                "CID",
                "HOST_RA",
                "HOST_DEC",
                "RA",
                "DEC",
                "zHEL",
                "AgeMW",
                "AgeMW_err",
                "global_logM",
                "AV",
            ]
        )
    allhosts.append(
        pd.DataFrame(
            {
                "source_id": "DustPedia:" + dd.CID,
                "physical_SN": dd.CID,
                "sample": "DustPedia",
                "host_ra": dd.HOST_RA,
                "host_dec": dd.HOST_DEC,
                "sn_ra": dd.RA,
                "sn_dec": dd.DEC,
                "z": dd.zHEL,
                "sed_age": dd.AgeMW,
                "sed_age_error": dd.AgeMW_err,
                "sed_mass": dd.global_logM,
                "sed_Av": dd.AV,
                "host_accepted": True,
            }
        )
    )
    h = pd.concat(allhosts, ignore_index=True)
    h["valid_position_redshift"] = (
        np.isfinite(h[["host_ra", "host_dec", "z"]]).all(axis=1)
        & h.host_ra.between(0, 360)
        & h.host_dec.between(-90, 90)
        & (h.z > 0)
    )
    return h, paths


def match(h, c, ra, dec, z, col):
    good = (
        np.isfinite(c[[ra, dec, z]]).all(axis=1)
        & c[ra].between(0, 360)
        & c[dec].between(-90, 90)
        & (c[z] > 0)
    )
    c = c[good].copy()
    c["catalogue_row"] = c.index
    hh = h[h.valid_position_redshift].copy()
    hh["host_row"] = hh.index
    hc = SkyCoord(hh.host_ra.to_numpy() * u.deg, hh.host_dec.to_numpy() * u.deg)
    cc = SkyCoord(c[ra].to_numpy() * u.deg, c[dec].to_numpy() * u.deg)
    i, j, sep, _ = search_around_sky(hc, cc, 3 * u.arcsec)
    q = (
        hh.iloc[i]
        .reset_index(drop=True)
        .join(c.iloc[j].reset_index(drop=True).add_prefix(col + "_"))
    )
    q["separation_arcsec"] = sep.arcsec
    q["delta_z"] = q[col + "_" + z] - q.z
    q["position_z_accepted"] = (abs(q.delta_z) < 0.001) & q.host_accepted
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", type=Path, default=ARCHIVE)
    args = ap.parse_args()
    h, paths = hosts(args.archive)
    h.to_csv(WORK / "host-coordinates.csv", index=False)
    summary = {
        "host_rows_by_sample": h.groupby("sample").size().to_dict(),
        "catalogues": {},
    }
    with fits.open(WORK / "gal_info_dr7_v5_2.fit.gz") as f:
        c = frame(f[1].data)
        c["PHOTOID"] = ["-".join(map(str, x)) for x in f[1].data["PHOTOID"]]
    q = match(h, c, "RA", "DEC", "Z", "mpa")
    q["quality_accepted"] = (
        q.position_z_accepted
        & (q.mpa_Z_WARNING == 0)
        & q.mpa_SPECTROTYPE.eq("GALAXY")
        & (q.mpa_SN_MEDIAN > 0)
    )
    # Repeat spectra of the same SDSS photometric object are one counterpart.
    q["unique_photometric_counterpart"] = (
        q.groupby("source_id").mpa_PHOTOID.transform("nunique") == 1
    )
    q["accepted"] = q.quality_accepted & q.unique_photometric_counterpart
    q.to_csv(WORK / "mpajhu-all-candidates.csv", index=False)
    selected = (
        q[q.accepted]
        .sort_values("mpa_SN_MEDIAN", ascending=False)
        .drop_duplicates("source_id")
        .sort_values("source_id")
    )
    selected.to_csv(WORK / "mpajhu-selected-hosts.csv", index=False)
    summary["catalogues"]["mpajhu"] = {
        "catalogue_rows": len(c),
        "candidate_spectra": len(q),
        "accepted_best_spectra": len(selected),
        "accepted_by_sample": selected.groupby("sample").size().to_dict(),
        "unique_SDSS_photometric_objects": selected.mpa_PHOTOID.nunique(),
        "unique_spectrum_keys": selected[["mpa_PLATEID", "mpa_MJD", "mpa_FIBERID"]]
        .drop_duplicates()
        .shape[0],
        "ambiguous_host_rows": q.loc[
            ~q.unique_photometric_counterpart, "source_id"
        ].nunique(),
    }
    for lib in ["miles", "mastar"]:
        p = WORK / f"manga-firefly-globalprop-v3_1_1-{lib}.fits"
        paths.append(p)
        with fits.open(p) as f:
            c = frame(f[1].data).join(frame(f[2].data)).join(frame(f[3].data))
        qq = match(h, c, "OBJRA", "OBJDEC", "REDSHIFT", "ff")
        qq["unique_galaxy"] = (
            qq.groupby("source_id").ff_MANGAID.transform("nunique") == 1
        )
        qq["accepted"] = qq.position_z_accepted & qq.unique_galaxy
        qq.to_csv(WORK / f"firefly-{lib}-all-candidates.csv", index=False)
        ss = (
            qq[qq.accepted]
            .sort_values(["source_id", "ff_PLATEIFU"])
            .drop_duplicates("source_id")
        )
        ss.to_csv(WORK / f"firefly-{lib}-selected-hosts.csv", index=False)
        summary["catalogues"][lib] = {
            "candidate_rows": len(qq),
            "accepted_host_rows": len(ss),
            "accepted_by_sample": ss.groupby("sample").size().to_dict(),
            "unique_MANGA_galaxies": ss.ff_MANGAID.nunique(),
        }
    paths.append(WORK / "gal_info_dr7_v5_2.fit.gz")
    summary["provenance"] = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "code_sha256": sha(Path(__file__)),
        "design_sha256": sha(Path(__file__).with_name("design.json")),
        "inputs": [{"name": p.name, "sha256": sha(p)} for p in paths],
        "output_hashes": {p.name: sha(p) for p in WORK.glob("*hosts.csv")},
    }
    (OUT / "crossmatch.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
