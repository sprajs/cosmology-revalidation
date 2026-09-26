#!/usr/bin/env python3
"""Independent, metadata-only audit of original DES RAISIN SIMLIB LIBID 11."""
import csv
import hashlib
import json
import math
from pathlib import Path
import struct
from collections import Counter

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
SOURCE = ROOT / "runs/research_2026_09_26/raisin_sign_source/sim/simlibs/DES_RAISIN.simlib"
CONFIG = ROOT / "runs/research_2026_09_26/raisin_profile_solver_review/csp-filter-interpretation/native-design/timing-followon-design/prospective-engineering/generation/original/sim.input"


def f32(x):
    return struct.unpack("f", struct.pack("f", x))[0]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_libid():
    inside = False
    meta = {}
    rows = []
    for lineno, line in enumerate(SOURCE.read_text().splitlines(), 1):
        if line.strip() == "LIBID: 11":
            if inside:
                raise ValueError("duplicate LIBID")
            inside = True
            continue
        if line.strip() == "END_LIBID: 11":
            inside = False
            break
        if not inside:
            continue
        words = line.split()
        if not words:
            continue
        if words[0] == "RA:":
            meta["ra"] = float(words[1])
            meta["dec"] = float(words[3])
            meta["nobs_header"] = int(words[5])
        elif words[0] == "MWEBV:":
            meta["mwebv"] = float(words[1])
            meta["pixsize"] = float(words[3])
        elif words[0] == "REDSHIFT:":
            meta["source_z"] = float(words[1])
            meta["source_peak"] = float(words[3])
        elif words[0] == "FIELD:":
            meta["field"] = words[1]
        elif words[0] == "S:":
            if len(words) != 13:
                raise ValueError((lineno, len(words), line))
            values = list(map(float, words[4:]))
            rows.append({"line": lineno, "ordinal": len(rows) + 1,
                         "mjd_text": words[1], "mjd": float(words[1]),
                         "exposure_id": int(words[2]), "band": words[3],
                         "gain": values[0], "rdnoise": values[1], "skysig": values[2],
                         "psf1": values[3], "psf2": values[4], "psfrat": values[5],
                         "zpavg": values[6], "zperr": values[7], "mag": values[8],
                         "raw_line": line})
    if inside or len(rows) != meta.get("nobs_header"):
        raise ValueError("LIBID11 block/end/NOBS mismatch")
    return meta, rows


def read_config():
    d = {}
    for line in CONFIG.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if ":" in line:
            key, value = line.split(":", 1)
            d[key.strip()] = value.strip()
    for key, value in {"USE_SIMLIB_REDSHIFT": "1", "USE_SIMLIB_PEAKMJD": "1", "VEL_CMBAPEX": "0", "GENSIGMA_VPEC": "0", "VPEC_ERR": "0"}.items():
        if d.get(key) != value:
            raise ValueError(("prospective configuration changed", key, d.get(key)))
    if any(k.startswith("HOSTLIB_") for k in d):
        raise ValueError("new HOSTLIB override needs redshift review")
    return d


def range_stat(rows, field):
    x = [r[field] for r in rows]
    return {"min": min(x), "max": max(x), "nonfinite": sum(not math.isfinite(v) for v in x),
            "negative": sum(v < 0 for v in x), "zero": sum(v == 0 for v in x)}


def main():
    meta, rows = read_libid()
    read_config()
    if meta["source_z"] != .453 or meta["source_peak"] != 57707.8:
        raise ValueError("LIBID11 source metadata changed")
    if len(rows) != 735 or sorted(r["exposure_id"] for r in rows) != list(range(735)):
        raise ValueError("original physical rows/IDs changed")
    if meta["pixsize"] <= 0 or not math.isfinite(meta["pixsize"]):
        raise ValueError("pixel scale invalid")
    zhel = f32(meta["source_z"])  # source-assigned CMB, equal to HEL only for verified prospective VEL/VPEC=0
    peak = f32(meta["source_peak"])
    if peak != 57707.80078125:
        raise ValueError("R4 peak mismatch")
    mult = Counter((r["mjd_text"], r["band"]) for r in rows)
    mult_mjd = Counter(r["mjd_text"] for r in rows)
    for r in rows:
        r["phase_source"] = (r["mjd"] - meta["source_peak"]) / (1 + meta["source_z"])
        r["phase_nominal"] = (r["mjd"] - peak) / (1 + zhel)
        r["phase_peak_minus4"] = (r["mjd"] - (peak - 4)) / (1 + zhel)
        r["phase_peak_plus4"] = (r["mjd"] - (peak + 4)) / (1 + zhel)
        r["selected_nominal"] = -15 <= r["phase_nominal"] <= 45
        r["selected_minus4"] = -15 <= r["phase_peak_minus4"] <= 45
        r["selected_plus4"] = -15 <= r["phase_peak_plus4"] <= 45
        r["selected_any_in_domain"] = r["phase_peak_plus4"] <= 45 and r["phase_peak_minus4"] >= -15
        r["selected_all_in_domain"] = r["phase_peak_plus4"] >= -15 and r["phase_peak_minus4"] <= 45
        r["mjd_band_multiplicity"] = mult[(r["mjd_text"], r["band"])]
        r["mjd_multiplicity"] = mult_mjd[r["mjd_text"]]
    selected = [r for r in rows if r["selected_nominal"]]
    counts = {k: dict(Counter(r["band"] for r in rows if r[k])) for k in ("selected_nominal", "selected_minus4", "selected_plus4", "selected_any_in_domain", "selected_all_in_domain")}
    if counts["selected_nominal"] != {"J": 3, "H": 3, "g": 27, "r": 28, "i": 28, "z": 28}:
        raise ValueError(("reported 117 phase support mismatch", counts["selected_nominal"]))
    if any(not math.isfinite(r[f]) for r in rows for f in ("mjd", "gain", "rdnoise", "skysig", "psf1", "psf2", "psfrat", "zpavg", "zperr", "mag")):
        raise ValueError("nonfinite SIMLIB metadata")
    for name in ("gain", "psf1", "zpavg"):
        if any(r[name] <= 0 for r in rows):
            raise ValueError(("expected positive SIMLIB metadata", name))
    for name in ("rdnoise", "skysig", "psf2", "psfrat", "zperr"):
        if any(r[name] < 0 for r in rows):
            raise ValueError(("expected nonnegative SIMLIB metadata", name))
    summary = {"source": meta, "rows": len(rows), "redshift_helio_conditional": zhel,
               "peak_r4": peak, "nominal_peak_shift": peak - meta["source_peak"],
               "counts": counts,
               "selected_nominal_phase_extrema": [min(r["phase_nominal"] for r in selected), max(r["phase_nominal"] for r in selected)],
               "nominal_selected_phase_over_peak_domain_extrema": [min(r["phase_peak_plus4"] for r in selected), max(r["phase_peak_minus4"] for r in selected)],
               "fit_guard_minus20_plus70_all_nominal_selected_over_peak_domain": all(-20 <= r["phase_peak_plus4"] and r["phase_peak_minus4"] <= 70 for r in selected),
               "unique_exposure_ids": len(set(r["exposure_id"] for r in rows)),
               "unique_mjd_band": len(mult), "duplicated_mjd_band_keys": sum(v > 1 for v in mult.values()),
               "max_mjd_band_multiplicity": max(mult.values()), "unique_mjd": len(mult_mjd),
               "max_mjd_multiplicity": max(mult_mjd.values()),
               "metadata_ranges": {name: range_stat(rows, name) for name in ("mjd", "gain", "rdnoise", "skysig", "psf1", "psf2", "psfrat", "zpavg", "zperr", "mag")}}
    HERE.mkdir(parents=True, exist_ok=True)
    with (HERE / "physical-epochs.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    (HERE / "result.json").write_text(json.dumps(summary, indent=2) + "\n")
    paths = [SOURCE, CONFIG, HERE / "protocol.json", Path(__file__)]
    (HERE / "manifest.json").write_text(json.dumps({str(p.relative_to(ROOT)): digest(p) for p in paths}, indent=2) + "\n")


if __name__ == "__main__":
    main()
