#!/usr/bin/env python3
"""Header-only CALWF3 RAW→FLT prerequisite inventory; never reads FITS data."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from shutil import which

from astropy.io import fits

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
RAW_DIR = BASE / "dark_ramp_raw_pilot/files"
FLT_DIR = BASE / "pixel_acquisition"
REF_DIRS = (
    BASE / "quality_reference_audit/reference_files",
    BASE / "calwf3_variance_source",
    BASE / "pixel_acquisition/pam_validation",
)
ROOTS = ("icxoi1bcq", "icxoi4hgq")
REF_KEYS = (
    "BPIXTAB", "CCDTAB", "OSCNTAB", "CRREJTAB", "DARKFILE",
    "NLINFILE", "PFLTFILE", "DFLTFILE", "IMPHTTAB", "IDCTAB", "MDRIZTAB",
)
SWITCH_KEYS = (
    "DQICORR", "ZSIGCORR", "ZOFFCORR", "DARKCORR", "BLEVCORR",
    "NLINCORR", "FLATCORR", "CRCORR", "UNITCORR", "PHOTCORR",
    "DRIZCORR", "RPTCORR",
)


def digest(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            sha.update(chunk)
    return sha.hexdigest()


def main() -> None:
    objects = {}
    refs = {}
    for root in ROOTS:
        raw = RAW_DIR / f"{root}_raw.fits"
        flt = FLT_DIR / f"{root}_flt.fits"
        rh = fits.getheader(raw, 0)
        fh = fits.getheader(flt, 0)
        objects[root] = {
            "raw": {"path": str(raw.resolve()), "bytes": raw.stat().st_size, "sha256": digest(raw)},
            "flt": {"path": str(flt.resolve()), "bytes": flt.stat().st_size, "sha256": digest(flt)},
            "archive_cal_ver": fh.get("CAL_VER"),
            "archive_crds_ctx": fh.get("CRDS_CTX"),
            "raw_switches": {key: rh.get(key) for key in SWITCH_KEYS},
            "archive_flt_switches": {key: fh.get(key) for key in SWITCH_KEYS},
            "raw_refs": {key: rh.get(key) for key in REF_KEYS},
        }
        for key in REF_KEYS:
            value = rh.get(key)
            if not isinstance(value, str) or not value.startswith("iref$"):
                raise ValueError(f"{root}: {key} is not an explicit iref$ name: {value!r}")
            name = value.split("$", 1)[1]
            record = refs.setdefault(name, {"used_by": [], "local": []})
            record["used_by"].append({"root": root, "key": key})
    for name, record in refs.items():
        for directory in REF_DIRS:
            match = directory / name
            if match.is_file():
                record["local"].append({"path": str(match.resolve()), "bytes": match.stat().st_size, "sha256": digest(match)})
    head = json.loads((HERE / "head-results.json").read_text())
    for row in head["rows"]:
        refs[row["name"]]["official_head"] = row
    result = {
        "purpose": "header/runtime feasibility only; no SCI/ERR/DQ arrays read",
        "objects": objects,
        "refs": refs,
        "binaries": {name: which(name) for name in ("calwf3.e", "wf3ir.e", "cmake", "gcc", "gfortran", "make", "pkg-config")},
        "head_results_sha256": digest(HERE / "head-results.json"),
    }
    out = HERE / "inventory.json"
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(out)
    print(digest(out))


if __name__ == "__main__":
    main()
