#!/usr/bin/env python3
"""Acquire legacy SDSS observed spectra and independently integrate Dn4000."""

from pathlib import Path
import concurrent.futures, datetime, json, time, urllib.request
import numpy as np, pandas as pd
from astropy.io import fits
import extinction
from acquire import ROOT, WORK, OUT, sha, check_source


def get(row):
    pnum, mjd, fib = int(row.plateid), int(row.mjd), int(row.fiberid)
    name = f"spec-{pnum:04d}-{mjd}-{fib:04d}.fits"
    p = WORK / "spectra" / name
    u = f"https://data.sdss.org/sas/dr17/sdss/spectro/redux/26/spectra/lite/{pnum:04d}/{name}"
    r = {
        "file": "spectra/" + name,
        "url": u,
        "plateid": pnum,
        "mjd": mjd,
        "fiberid": fib,
        "access_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    }
    try:
        if not p.exists():
            for attempt in range(3):
                try:
                    b = urllib.request.urlopen(u, timeout=45).read()
                    p.with_suffix(".partial").write_bytes(b)
                    with fits.open(p.with_suffix(".partial")) as f:
                        assert f[1].name == "COADD"
                    p.with_suffix(".partial").replace(p)
                    break
                except Exception:
                    if attempt == 2:
                        raise
                    time.sleep(1)
            r["access"] = "downloaded"
        else:
            r["access"] = "reused"
        check_source(p)
        r.update(sha256=sha(p), bytes=p.stat().st_size)
    except Exception as e:
        r["error"] = str(e)
    return r


def pixel_weights(wave, lo, hi):
    # Integrate actual partially intersected pixels in rest wavelength; bins in
    # log wavelength are not equal-width wavelength samples.
    edges = np.r_[
        wave[0] - (wave[1] - wave[0]) / 2,
        (wave[1:] + wave[:-1]) / 2,
        wave[-1] + (wave[-1] - wave[-2]) / 2,
    ]
    return np.maximum(0, np.minimum(edges[1:], hi) - np.maximum(edges[:-1], lo))


def measure(p, z, ebv):
    with fits.open(p) as ff:
        d = ff[1].data
        wave = 10 ** np.asarray(d["loglam"], float)
        wr = wave / (1 + z)
        fl = np.asarray(d["flux"], float)
        iv = np.asarray(d["ivar"], float)
        mask = np.asarray(d["and_mask"])
        meta = ff[2].data
    good = np.isfinite(fl) & np.isfinite(iv) & (iv > 0) & (mask == 0)
    # D_n uses mean Fnu, proportional to wavelength^2 * Flambda. Common (1+z)
    # and speed-of-light normalizations cancel in the ratio.
    B = pixel_weights(wr, 3850, 3950)
    R = pixel_weights(wr, 4000, 4100)
    coverage = [float(np.sum(w * good) / 100) for w in [B, R]]
    if min(coverage) < 0.9:
        return {"valid": False, "coverage_blue_red": coverage}
    B = B * good
    R = R * good
    B /= B.sum()
    R /= R.sum()
    weighted = fl * wr**2
    out = {
        "valid": True,
        "coverage_blue_red": coverage,
        "blue_coverage_fraction": coverage[0],
        "red_coverage_fraction": coverage[1],
        "blue_valid_pixels": int(np.count_nonzero(B)),
        "red_valid_pixels": int(np.count_nonzero(R)),
        "spec_z": float(meta["Z"][0]),
    }
    for lab, attenuation in [
        ("raw", np.ones(len(wave))),
        ("MW_corrected", 10 ** (0.4 * extinction.odonnell94(wave, 3.1 * ebv, 3.1))),
    ]:
        f = weighted * attenuation
        wB = B * wr**2 * attenuation
        wR = R * wr**2 * attenuation
        b = B @ f
        r = R @ f
        variance = np.zeros(len(iv))
        variance[good] = 1 / iv[good]
        vb = np.sum(wB * wB * variance)
        vr = np.sum(wR * wR * variance)
        # Observed Fnu in microJy; rest-frame wavelength band integration is a
        # normalized weighting, not a rest-frame luminosity conversion.
        unit = (1 + z) ** 2 * 1e12 / 2.99792458e18
        out.update(
            {
                lab + "_blue_Fnu_uJy": b * unit,
                lab + "_red_Fnu_uJy": r * unit,
                lab + "_blue_variance_uJy2": vb * unit**2,
                lab + "_red_variance_uJy2": vr * unit**2,
                lab + "_blue_red_covariance_uJy2": 0.0,
            }
        )
        if b <= 0 or r <= 0:
            out["valid"] = False
            out["invalid_reason"] = "nonpositive band flux"
            continue
        out[lab + "_Dn4000"] = r / b
        out[lab + "_Dn4000_error"] = np.sqrt(vr / (b * b) + r * r * vb / (b**4))
    return out


def main():
    (WORK / "spectra").mkdir(exist_ok=True)
    p = WORK / "sdss-spectral-measurements.csv"
    input_sha = sha(p)
    d = pd.read_csv(
        p, dtype={"specObjID": str, "bestObjID": str, "photometric_objID": str}
    )
    d = d.drop_duplicates(["plateid", "mjd", "fiberid"])
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        records = list(pool.map(get, d.itertuples(index=False)))
    h = pd.read_csv(WORK / "mpajhu-selected-hosts.csv")
    h = h.rename(
        columns={"mpa_PLATEID": "plateid", "mpa_MJD": "mjd", "mpa_FIBERID": "fiberid"}
    ).drop_duplicates(["plateid", "mjd", "fiberid"])
    joined = d.merge(
        h[["plateid", "mjd", "fiberid", "mpa_E_BV_SFD"]],
        on=["plateid", "mjd", "fiberid"],
        validate="one_to_one",
    )
    bykey = {
        (int(r.plateid), int(r.mjd), int(r.fiberid)): r
        for r in joined.itertuples(index=False)
    }
    measurements = []
    for r in records:
        if "error" in r:
            continue
        key = (r["plateid"], r["mjd"], r["fiberid"])
        q = bykey[key]
        try:
            m = measure(WORK / r["file"], q.spectral_z, q.mpa_E_BV_SFD)
        except Exception as e:
            m = {"valid": False, "error": str(e)}
        measurements.append(dict(plateid=key[0], mjd=key[1], fiberid=key[2], **m))
    pd.DataFrame(measurements).to_csv(WORK / "remeasured-dn4000.csv", index=False)
    # Equation checks: flat Fnu -> 1; flat Flambda has analytic lambda^2 band ratio.
    w = np.linspace(3800, 4150, 35001)
    b = pixel_weights(w, 3850, 3950)
    r = pixel_weights(w, 4000, 4100)
    b /= b.sum()
    r /= r.sum()
    flatnu = (r @ np.ones(len(w))) / (b @ np.ones(len(w)))
    flatlambda = (r @ (w * w)) / (b @ (w * w))
    exact = (4100**3 - 4000**3) / (3950**3 - 3850**3)
    assert abs(flatnu - 1) < 1e-12 and abs(flatlambda - exact) < 1e-8
    (OUT / "spectra-acquisition.json").write_text(
        json.dumps(
            {
                "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                "code_sha256": sha(Path(__file__)),
                "input_sha256": input_sha,
                "input_rows": len(d),
                "band_mean_units": "Observed Fnu microJy in rest-wavelength-selected 3850-3950 and 4000-4100 Angstrom bands; zero blue-red covariance assumes supplied diagonal resampled-pixel ivar. SFD E(B-V) is unscaled, ODonnell94 Rv3.1 for MW_corrected; raw means permit other forward models.",
                "requested": len(records),
                "downloaded_or_reused": sum("sha256" in r for r in records),
                "records": records,
                "remeasured_sha256": sha(WORK / "remeasured-dn4000.csv"),
                "equation_checks": {
                    "flat_Fnu_ratio": flatnu,
                    "flat_Flambda_ratio": flatlambda,
                    "analytic_flat_Flambda_ratio": exact,
                },
                "limits": "Observed resampled flux/ivar; pixel covariance not supplied. Legacy3arcsec fibre. Ratios recomputed with and without foreground ODonnell94 correction; no stellar continuum/emission-line model subtraction or velocity-dispersion correction.",
            },
            indent=2,
        )
        + "\n"
    )
    print(
        json.dumps(
            {
                "requested": len(records),
                "valid_Dn4000": sum(x.get("valid", False) for x in measurements),
                "download_failures": sum("error" in r for r in records),
            }
        )
    )


if __name__ == "__main__":
    main()
