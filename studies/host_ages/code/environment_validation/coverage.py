#!/usr/bin/env python3
"""Bounded physical-identity audit for DES and RAISIN age-test coverage.

Numeric DES CIDs are only joined to DES HEAD SNIDs, never interpreted as IAU
names. Cross-survey sky candidates also require compatible peak epochs.
"""
import argparse, datetime, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.io import fits
import astropy.units as u
from analyse import ROOT, HERE, DEFAULT_ARCHIVE, sha, normalize

sys.path.insert(0, str(ROOT))
from workflows.raisin import distance
from lib.photometry import read_raisin
from lib.records import fitres


def coordinates(frame, ra="ra", dec="dec"):
    return SkyCoord(frame[ra].to_numpy() * u.deg, frame[dec].to_numpy() * u.deg)


def sky_candidates(a, b, label):
    i, j, sep, _ = search_around_sky(coordinates(a), coordinates(b), 3 * u.arcsec)
    out = []
    for ia, ib, s in zip(i, j, sep.arcsec):
        delta = abs(float(a.iloc[ia].peak) - float(b.iloc[ib].peak))
        out.append(
            {
                "left": str(a.iloc[ia].name_id),
                "right": str(b.iloc[ib].name_id),
                "separation_arcsec": float(s),
                "peak_difference_days": delta,
                "accepted_same_event": delta <= 30,
                "label": label,
            }
        )
    return out


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
    archive = args.archive
    work = ROOT / ".work/environment-validation"
    base = archive / "sources/updates/2026-09-20-ztf"
    tp = args.titan
    t = pd.read_csv(tp)
    tkeys = set(t.transient.map(normalize))
    zp = base / "extracted/ztfsniadr2_lite/tables/snia_data.csv"
    z = pd.read_csv(zp)
    z["key"] = z.iau_name.map(normalize)
    z = z[z.key.isin(tkeys)].copy()
    z = z.rename(columns={"ztfname": "name_id", "t0": "peak"})
    z = z[np.isfinite(z[["ra", "dec", "peak"]]).all(axis=1)].reset_index(drop=True)
    hp = archive / "data/des-diffimg/DES-SN5YR_DIFFIMG_HEAD.FITS.gz"
    with fits.open(hp) as h:
        a = h[1].data
        des = pd.DataFrame(
            {
                "name_id": a["SNID"].astype(str),
                "alias": a["IAUC"].astype(str),
                "ra": a["RA"],
                "dec": a["DEC"],
                "peak": a["PEAKMJD"],
            }
        )
    mp = (
        archive
        / "sources/repos/des-science__DES-SN5YR/4_DISTANCES_COVMAT/DES-Dovekie_Metadata.csv"
    )
    meta = fitres(mp)
    mids = set(meta.index.astype(str))
    ds = sky_candidates(z, des, "TITAN via ZTF to DES HEAD")
    desliteral = [
        {"SNID": str(r.name_id), "alias": str(r.alias)}
        for _, r in des.iterrows()
        if normalize(r.alias) in tkeys
    ]
    desdirect = [
        str(k) for k in meta.index if not str(k).isdigit() and normalize(k) in tkeys
    ]
    for row in ds:
        row["in_Dovekie_distance_table"] = str(row["right"]) in mids
    rp = ROOT / "data/raisin/distances/w/nir_dist/RAISIN_combined_FITOPT000.FITRES"
    nir = distance("nir")
    opt = distance("optical")
    both = distance("opticalnir")
    assert list(nir.index) == list(opt.index) == list(both.index) and len(nir) == 79
    photpaths = []
    ph = []
    for p in sorted((ROOT / "data/raisin/photometry/RAISIN").rglob("*")):
        if not p.is_file() or p.suffix.lower() != ".dat":
            continue
        hdr, _ = read_raisin(p)
        cid = hdr["SNID"].split()[0]
        if cid not in nir.index:
            continue
        ph.append(
            {
                "name_id": cid,
                "ra": float(hdr["RA"].split()[0]),
                "dec": float(hdr.get("DEC", hdr.get("DECL")).split()[0]),
                "peak": float(hdr["PEAKMJD"].split()[0]),
            }
        )
        photpaths.append(p)
    r = pd.DataFrame(ph)
    assert r.name_id.is_unique and len(r) == 79
    lpath = base / "extracted/kelsey2026-stag1765-supplement/photometry_local_3kpc.txt"
    l = pd.read_csv(lpath, comment="#")
    lkeys = set(l.SN.map(normalize))
    coord_l = coordinates(l, "RA", "Dec")
    ri, li, sep, _ = search_around_sky(coordinates(r), coord_l, 3 * u.arcsec)
    dustsky = [
        {
            "RAISIN": str(r.iloc[i].name_id),
            "DustPedia": str(l.iloc[j].SN),
            "separation_arcsec": float(s),
            "same_normalized_identifier": normalize(r.iloc[i].name_id)
            == normalize(l.iloc[j].SN),
        }
        for i, j, s in zip(ri, li, sep.arcsec)
    ]
    pp = ROOT / "data/distances/Pantheon+SH0ES.dat"
    pan = pd.read_csv(pp, sep=r"\s+")
    pan["key"] = pan.CID.map(normalize)
    pall = pan.rename(
        columns={"CID": "name_id", "RA": "ra", "DEC": "dec", "PKMJD": "peak"}
    )
    pal = pall[pall.key.isin(tkeys)].copy().reset_index(drop=True)
    rztf = sky_candidates(r, z, "RAISIN to TITAN via ZTF IAU")
    rpan = sky_candidates(r, pal, "RAISIN to TITAN via Pantheon alias")
    rall = sky_candidates(
        r, pall, "RAISIN to all Pantheon (age availability diagnostic)"
    )
    result = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "code_sha256": sha(__file__),
        "dependencies_sha256": {
            str(p.relative_to(ROOT)): sha(p)
            for p in [
                HERE / "analyse.py",
                ROOT / "lib/records.py",
                ROOT / "lib/photometry.py",
                ROOT / "workflows/raisin.py",
            ]
        },
        "method": "Exact named-event joins plus <=3 arcsec coordinate bridges with <=30-day peak agreement. A repeated host or sky transient at another epoch is not the same SN.",
        "inputs": [
            {"path": str(p), "sha256": sha(p)} for p in [tp, zp, hp, mp, rp, lpath, pp]
        ],
        "raisin_selected_header_sha256": {
            str(p.relative_to(ROOT)): sha(p) for p in photpaths
        },
        "TITAN_rows": len(t),
        "TITAN_unique_names": len(tkeys),
        "TITAN_to_ZTF_physical_name_bridge": len(z),
        "DES": {
            "HEAD_rows": len(des),
            "Dovekie_metadata_rows": len(meta),
            "TITAN_exact_IAUC_alias_matches": desliteral,
            "TITAN_non_numeric_Dovekie_CID_matches": desdirect,
            "ZTF_sky_candidates": ds,
            "clean_bridged_Dovekie_matches": [
                d
                for d in ds
                if d["accepted_same_event"] and d["in_Dovekie_distance_table"]
            ],
            "limits": "TITAN CSV has no sky coordinates. The ZTF bridge covers its matched named events only. No claim that every external alias in the universe was resolved.",
        },
        "RAISIN": {
            "common_distance_SNe": len(nir),
            "lowz_count": int((nir.zHD < 0.1).sum()),
            "highz_count": int((nir.zHD > 0.2).sum()),
            "direct_TITAN_names": [str(k) for k in nir.index if normalize(k) in tkeys],
            "direct_DustPedia_names": [
                str(k) for k in nir.index if normalize(k) in lkeys
            ],
            "DustPedia_sky_candidates": dustsky,
            "TITAN_ZTF_sky_candidates": rztf,
            "TITAN_Pantheon_sky_candidates": rpan,
            "all_Pantheon_same_event_match_count": len(
                {d["left"] for d in rall if d["accepted_same_event"]}
            ),
            "cross_band_covariance_limit": "The three released branch covariance matrices are not their joint optical-NIR cross covariance; the existing RAISIN workflow retains this missing cross-block as a bound, not zero.",
            "analysis_decision": "No new age-versus-optical-minus-NIR regression is justified unless a clean age-matched event cohort is found. Optical and NIR distances exist, but host age on the same physical events is the binding observation here.",
        },
        "validation": {
            "RAISIN_three_branch_identical_order": True,
            "RAISIN_unique_photometry_identity_count": len(r),
            "numeric_DES_CID_never_used_as_IAU": True,
        },
    }
    (
        ROOT / "studies/host_ages/results/environment_validation/coverage.json"
    ).write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    pd.DataFrame(ds + rztf + rpan + rall).to_csv(
        work / "coverage-sky-candidates.csv", index=False
    )
    print(
        json.dumps(
            {
                k: v
                for k, v in result.items()
                if k not in ["inputs", "raisin_selected_header_sha256"]
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
