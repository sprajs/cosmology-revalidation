#!/usr/bin/env python3
"""Acquire matched MaNGA maps and test the actual SN-site/central aperture gap.

Nearest-spaxel indices are spatially PSF-mixed, not independent local stellar
ages. No averaging of neighbouring spaxels as independent measurements.
"""

from pathlib import Path
import concurrent.futures, datetime, json, time, urllib.request
import numpy as np, pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy.stats import spearmanr
from acquire import ROOT, WORK, OUT, sha, check_source

TYPE = "HYB10-MILESHC-MASTARSSP"


def get(pi):
    plate, ifu = pi.split("-")
    p = WORK / "manga-maps" / f"manga-{pi}-MAPS.fits.gz"
    url = f"https://data.sdss.org/sas/dr17/manga/spectro/analysis/v3_1_1/3.1.0/{TYPE}/{plate}/{ifu}/manga-{pi}-MAPS-{TYPE}.fits.gz"
    record = {
        "plateifu": pi,
        "url": url,
        "path": str(p.relative_to(ROOT)),
        "access_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    }
    try:
        if not p.exists():
            for attempt in range(3):
                try:
                    b = urllib.request.urlopen(url, timeout=90).read()
                    p.with_suffix(".partial").write_bytes(b)
                    with fits.open(p.with_suffix(".partial")) as f:
                        assert "SPECINDEX" in f
                    p.with_suffix(".partial").replace(p)
                    break
                except Exception:
                    if attempt == 2:
                        raise
                    time.sleep(2)
            record["access"] = "downloaded"
        else:
            record["access"] = "reused"
        check_source(p)
        record.update(bytes=p.stat().st_size, sha256=sha(p))
    except Exception as e:
        record["error"] = str(e)
    print(json.dumps(record), flush=True)
    return record


def measure(pi, rows):
    p = WORK / "manga-maps" / f"manga-{pi}-MAPS.fits.gz"
    out = []
    with fits.open(p) as f:
        header = f["SPECINDEX"].header
        w = WCS(header).celestial
        channels = {
            v.strip(): int(k[1:]) - 1
            for k, v in header.items()
            if len(k) == 3 and k.startswith("C") and k[1:].isdigit()
        }
        for r in rows.itertuples():
            result = {
                "source_id": r.source_id,
                "sample": r.sample,
                "physical_SN": r.physical_SN,
                "MANGAID": r.ff_MANGAID,
                "plateifu": pi,
                "sed_age": r.sed_age,
                "z": r.z,
                "host_ra": r.host_ra,
                "host_dec": r.host_dec,
                "sn_ra": r.sn_ra,
                "sn_dec": r.sn_dec,
                "r_band_PSF_FWHM_arcsec": f[0].header.get("RFWHM", np.nan),
            }
            for lab, ra, dec in [
                ("central", r.host_ra, r.host_dec),
                ("SN_site", r.sn_ra, r.sn_dec),
            ]:
                result[lab + "_valid"] = False
                if not np.isfinite([ra, dec]).all():
                    result[lab + "_failure"] = "missing SN position"
                    continue
                x, y = w.world_to_pixel_values(ra, dec)
                ix, iy = int(np.round(x)), int(np.round(y))
                ny, nx = f["SPX_SNR"].data.shape
                result[lab + "_pixel_x"] = ix
                result[lab + "_pixel_y"] = iy
                if not (0 <= ix < nx and 0 <= iy < ny):
                    result[lab + "_failure"] = "outside map rectangle"
                    continue
                result[lab + "_spaxel_SNR"] = float(f["SPX_SNR"].data[iy, ix])
                result[lab + "_BINID"] = int(f["BINID"].data[-1, iy, ix])
                wr, wd = w.pixel_to_world_values(ix, iy)
                sep = (
                    SkyCoord(ra * u.deg, dec * u.deg)
                    .separation(SkyCoord(wr * u.deg, wd * u.deg))
                    .arcsec
                )
                result[lab + "_coordinate_error_arcsec"] = float(sep)
                valid = result[lab + "_spaxel_SNR"] > 3 and result[lab + "_BINID"] >= 0
                for idx in ["Dn4000", "HDeltaA"]:
                    c = channels[idx]
                    val = float(f["SPECINDEX"].data[c, iy, ix])
                    ivar = float(f["SPECINDEX_IVAR"].data[c, iy, ix])
                    mask = int(f["SPECINDEX_MASK"].data[c, iy, ix])
                    corr = float(f["SPECINDEX_CORR"].data[c, iy, ix])
                    good = (
                        np.isfinite([val, ivar, corr]).all()
                        and ivar > 0
                        and mask == 0
                        and corr > 0
                        and (val > 0 if idx == "Dn4000" else True)
                    )
                    result[lab + "_" + idx + "_raw"] = val
                    result[lab + "_" + idx + "_mask"] = mask
                    result[lab + "_" + idx + "_dispersion_correction"] = corr
                    if good:
                        result[lab + "_" + idx] = val * corr
                        result[lab + "_" + idx + "_error"] = abs(corr) / np.sqrt(ivar)
                    valid &= good
                result[lab + "_valid"] = bool(valid)
                if not valid:
                    result[lab + "_failure"] = (
                        "SNR<=3, uncovered, masked or nonpositive index inverse variance/correction"
                    )
            if np.isfinite([r.sn_ra, r.sn_dec]).all():
                result["SN_host_separation_arcsec"] = float(
                    SkyCoord(r.sn_ra * u.deg, r.sn_dec * u.deg)
                    .separation(SkyCoord(r.host_ra * u.deg, r.host_dec * u.deg))
                    .arcsec
                )
            out.append(result)
    return out


def summarise(d):
    from diagnostics import quantile, association

    out = {
        "host_rows": len(d),
        "distinct_galaxies": d.MANGAID.nunique(),
        "valid_central": int(d.central_valid.sum()),
        "valid_SN_site": int(d.SN_site_valid.sum()),
        "failed_SN_site": d.loc[~d.SN_site_valid, "SN_site_failure"]
        .value_counts()
        .to_dict(),
        "samples": {},
    }
    for sample in ["ZTF", "G11", "R19", "all"]:
        q = (d if sample == "all" else d[d["sample"] == sample]).drop_duplicates(
            "MANGAID"
        )
        paired = q[q.central_valid & q.SN_site_valid].copy()
        out["samples"][sample] = {
            "distinct_galaxies": len(q),
            "paired_central_SN": len(paired),
            "PSF_FWHM_arcsec": quantile(q.r_band_PSF_FWHM_arcsec),
            "SN_host_separation_arcsec": quantile(paired.SN_host_separation_arcsec),
        }
        for idx in ["Dn4000", "HDeltaA"]:
            out["samples"][sample][idx + "_SN_minus_central"] = quantile(
                paired["SN_site_" + idx] - paired["central_" + idx]
            )
            out["samples"][sample][idx + "_central_vs_SED_age"] = association(
                q[q.central_valid], "central_" + idx
            )
            out["samples"][sample][idx + "_SNsite_vs_SED_age"] = association(
                q[q.SN_site_valid], "SN_site_" + idx
            )
    return out


def main():
    (WORK / "manga-maps").mkdir(exist_ok=True)
    p = WORK / "firefly-miles-selected-hosts.csv"
    d = pd.read_csv(p)
    inputsha = sha(p)
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        records = list(pool.map(get, sorted(set(d.ff_PLATEIFU))))
    measurements = []
    for r in records:
        if "error" not in r:
            measurements += measure(r["plateifu"], d[d.ff_PLATEIFU.eq(r["plateifu"])])
    q = pd.DataFrame(measurements)
    q.to_csv(WORK / "manga-central-SN-indices.csv", index=False)
    out = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "code_sha256": sha(Path(__file__)),
        "summary_dependency_sha256": sha(Path(__file__).with_name("diagnostics.py")),
        "input_crosswalk_sha256": inputsha,
        "acquisition": records,
        "output_sha256": sha(WORK / "manga-central-SN-indices.csv"),
        "analysis": summarise(q),
        "limits": [
            "MaNGA DAP emission-line-subtracted indices, model-based velocity-dispersion correction; independent observing modality, not model-free age.",
            "Nearest0.5arcsec spaxel only, SNR>3, full index mask zero and positive ivar; no spatially independent-pixel claim. Seeing/PSF blends different local populations.",
            "Local means the observed SN coordinate, conditional on released astrometry; no inference of progenitor travel, birth site, line-of-sight extinction or exact age.",
            "Descriptive central-minus-SN distributions are not confidence intervals; age rank bootstrap intervals use distinct galaxies. Small overlapping heterogeneous samples are reported separately.",
        ],
    }
    (OUT / "manga-local-summary.json").write_text(
        json.dumps(out, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps(out["analysis"], indent=2))


if __name__ == "__main__":
    main()
