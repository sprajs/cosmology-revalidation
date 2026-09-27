#!/usr/bin/env python3
"""Fetch MPA-JHU diagnostics and native SDSS host fluxes for exact spectra."""

from pathlib import Path
import concurrent.futures, datetime, hashlib, io, json, time, urllib.request, urllib.parse
import pandas as pd
from acquire import ROOT, WORK, OUT, sha, check_source

ID_DTYPES = {
    k: "string"
    for k in ["specObjID", "bestObjID", "photometric_objID", "photometric_flags"]
}
ENDPOINT = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"
COLS = [
    "g.specObjID",
    "g.plateid",
    "g.mjd",
    "g.fiberid",
    "g.ra AS spectral_ra",
    "g.dec AS spectral_dec",
    "g.z AS spectral_z",
    "g.z_err",
    "g.reliable",
    "g.sn_median",
    "s.sciencePrimary",
    "s.bestObjID",
    "s.run2d",
    "s.run1d",
    "i.d4000_n",
    "i.d4000_n_err",
    "i.d4000_n_sub",
    "i.d4000_n_sub_err",
    "i.lick_hd_a",
    "i.lick_hd_a_err",
    "i.lick_hd_a_sub",
    "i.lick_hd_a_sub_err",
    "i.lick_mgb",
    "i.lick_mgb_err",
    "i.lick_fe5270",
    "i.lick_fe5270_err",
    "i.lick_fe5335",
    "i.lick_fe5335_err",
]
for line in ["h_alpha", "h_beta", "oiii_5007", "nii_6584"]:
    COLS += ["l." + line + s for s in ["_flux", "_flux_err", "_eqw", "_eqw_err"]]
COLS += [
    "p.objID AS photometric_objID",
    "p.ra AS photometric_ra",
    "p.dec AS photometric_dec",
    "p.flags AS photometric_flags",
    "p.clean AS photometric_clean",
    "p.type AS photometric_type",
    "p.petroRad_r",
    "p.psfFwhm_r",
]
for b in "ugriz":
    for prefix in [
        "modelFlux",
        "modelFluxIvar",
        "cModelFlux",
        "cModelFluxIvar",
        "fiberFlux",
        "fiberFluxIvar",
        "extinction",
    ]:
        COLS.append("p." + prefix + "_" + b)


def fetch(item):
    k, rows = item
    conditions = [
        "(g.plateid="
        + str(int(r[0]))
        + " AND g.mjd="
        + str(int(r[1]))
        + " AND g.fiberid="
        + str(int(r[2]))
        + ")"
        for r in rows
    ]
    sql = (
        "SELECT "
        + ",".join(COLS)
        + " FROM galSpecInfo g JOIN galSpecIndx i ON i.specObjID=g.specObjID JOIN galSpecLine l ON l.specObjID=g.specObjID JOIN SpecObjAll s ON s.specObjID=g.specObjID LEFT JOIN PhotoObjAll p ON p.objID=s.bestObjID WHERE "
        + " OR ".join(conditions)
    )
    stem = (
        WORK
        / "sdss-queries"
        / ("query-" + hashlib.sha256(sql.encode()).hexdigest()[:20])
    )
    stem.parent.mkdir(exist_ok=True)
    stem.with_suffix(".sql").write_text(sql + "\n")
    p = stem.with_suffix(".csv")
    r = {
        "batch": k,
        "query_file": stem.name,
        "requested_spectra": len(rows),
        "sql_sha256": sha(stem.with_suffix(".sql")),
        "endpoint": ENDPOINT,
        "access_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    }
    try:
        if not p.exists():
            u = ENDPOINT + "?" + urllib.parse.urlencode({"cmd": sql, "format": "csv"})
            b = None
            for attempt in range(3):
                try:
                    b = urllib.request.urlopen(
                        urllib.request.Request(
                            u,
                            headers={
                                "User-Agent": "cosmology-revalidation public galaxy spectrum audit"
                            },
                        ),
                        timeout=90,
                    ).read()
                    break
                except Exception:
                    if attempt == 2:
                        raise
                    time.sleep(2)
            d = pd.read_csv(io.BytesIO(b), comment="#", dtype=ID_DTYPES)
            if "d4000_n" not in d:
                raise ValueError(
                    "No spectral-column CSV returned: "
                    + b[:250].decode(errors="replace")
                )
            p.write_bytes(b)
            r["access"] = "downloaded"
        else:
            r["access"] = "reused"
        check_source(p)
        d = pd.read_csv(p, comment="#", dtype=ID_DTYPES)
        r.update(rows=len(d), sha256=sha(p))
        return r, d
    except Exception as e:
        r["error"] = str(e)
        return r, None


def main():
    p = WORK / "mpajhu-selected-hosts.csv"
    d = pd.read_csv(p)
    keys = (
        d[["mpa_PLATEID", "mpa_MJD", "mpa_FIBERID"]]
        .drop_duplicates()
        .sort_values(["mpa_PLATEID", "mpa_MJD", "mpa_FIBERID"])
        .to_numpy()
    )
    items = [(i // 20, keys[i : i + 20]) for i in range(0, len(keys), 20)]
    records = []
    frames = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        for r, q in pool.map(fetch, items):
            records.append(r)
            if q is not None:
                frames.append(q)
            print(json.dumps(r), flush=True)
    if frames:
        cat = pd.concat(frames, ignore_index=True)
        for field in ID_DTYPES:
            assert cat[field].dropna().str.fullmatch(r"[0-9]+").all(), field
        assert cat.specObjID.is_unique
        cat.to_csv(WORK / "sdss-spectral-measurements.csv", index=False)
    else:
        cat = pd.DataFrame()
    report = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "code_sha256": sha(Path(__file__)),
        "host_crosswalk_sha256": sha(p),
        "records": records,
        "csv_read_dtype": ID_DTYPES,
        "identifier_validation": "Exact decimal digit strings from original query bytes, including nullable photometric identifiers and flags; never reconstructed from floats. Consumers must use csv_read_dtype.",
        "requested_unique_spectra": len(keys),
        "returned_rows": len(cat),
        "output_sha256": (
            sha(WORK / "sdss-spectral-measurements.csv") if len(cat) else None
        ),
        "units": "Native model/cmodel/fiberFlux columns in SDSS nanomaggies; IVAR columns in inverse flux squared. No logarithm/asinh magnitude conversion applied. Fibre diameter3arcsec for SDSS legacy spectra. MPA-JHU published spectral processing still model-dependent.",
    }
    (OUT / "sdss-acquisition.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
