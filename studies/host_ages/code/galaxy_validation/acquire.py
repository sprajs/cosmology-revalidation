#!/usr/bin/env python3
"""Acquire independent spectral catalogues, recording every HTTP result and hash."""

from pathlib import Path
import argparse, concurrent.futures, datetime, hashlib, json, subprocess

ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "studies/host_ages/results/galaxy_validation"
WORK = ROOT / ".work/galaxy-validation"
URLS = {
    "R19-ReadMe": "https://cdsarc.cds.unistra.fr/ftp/J/ApJ/874/32/ReadMe",
    "R19-table7.dat": "https://cdsarc.cds.unistra.fr/ftp/J/ApJ/874/32/table7.dat",
    "titan-host-properties.csv": "https://raw.githubusercontent.com/SterlingYM/age-of-titans/9b9ed9e5faf8b2f6267c0dce8cea97704d9b95ea/SFH_Iyer/host_props_with_SN_age_good_Mar18.csv",
    "R19-campbell_global.tsv": "https://raw.githubusercontent.com/benjaminrose/MC-Age/92713be96a89da991fe53bffcc596a5c0942fc37/data/campbell_global.tsv",
    "gal_info_dr7_v5_2.fit.gz": "https://wwwmpa.mpa-garching.mpg.de/SDSS/DR7/Data/gal_info_dr7_v5_2.fit.gz",
    "gal_indx_dr7_v5_2.fit.gz": "https://wwwmpa.mpa-garching.mpg.de/SDSS/DR7/Data/gal_indx_dr7_v5_2.fit.gz",
    "gal_line_dr7_v5_2.fit.gz": "https://wwwmpa.mpa-garching.mpg.de/SDSS/DR7/Data/gal_line_dr7_v5_2.fit.gz",
    "manga-firefly-globalprop-v3_1_1-miles.fits": "https://data.sdss.org/sas/dr17/manga/spectro/firefly/v3_1_1/manga-firefly-globalprop-v3_1_1-miles.fits",
    "manga-firefly-globalprop-v3_1_1-mastar.fits": "https://data.sdss.org/sas/dr17/manga/spectro/firefly/v3_1_1/manga-firefly-globalprop-v3_1_1-mastar.fits",
    "mpajhu-methods.html": "https://www.sdss4.org/dr17/spectro/galaxy_mpajhu/",
    "firefly-methods.html": "https://www.sdss4.org/dr17/manga/manga-data/manga-firefly-value-added-catalog/",
    "mpajhu-dr7.html": "https://wwwmpa.mpa-garching.mpg.de/SDSS/DR7/",
}


def sha(p):
    with p.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def check_source(path):
    lockfile = Path(__file__).with_name("source-lock.json")
    if lockfile.exists():
        expected = json.loads(lockfile.read_text())["files"].get(
            str(path.relative_to(WORK))
        )
        if expected is not None and sha(path) != expected:
            raise ValueError("Source hash changed: " + str(path.relative_to(WORK)))


def get(name):
    p = WORK / name
    u = URLS[name]
    r = {
        "file": name,
        "url": u,
        "access_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    }
    if p.exists():
        r["access"] = "reused existing bytes"
    else:
        tmp = p.with_suffix(p.suffix + ".partial")
        cmd = [
            "curl",
            "--silent",
            "--show-error",
            "--fail",
            "--location",
            "--retry",
            "3",
            "--connect-timeout",
            "30",
            "--max-time",
            "900",
            "--output",
            str(tmp),
            u,
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode:
            r["error"] = proc.stderr.strip()
            return r
        tmp.replace(p)
        r["access"] = "downloaded"
    check_source(p)
    r.update(bytes=p.stat().st_size, sha256=sha(p))
    return r


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--full", action="store_true")
    q = a.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    names = (
        list(URLS)
        if q.full
        else [x for x in URLS if "gal_indx" not in x and "gal_line" not in x]
    )
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        records = list(pool.map(get, names))
    for r in records:
        print(json.dumps(r), flush=True)
    (OUT / "acquisition.json").write_text(
        json.dumps(
            {
                "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                "code_sha256": sha(Path(__file__)),
                "records": records,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
