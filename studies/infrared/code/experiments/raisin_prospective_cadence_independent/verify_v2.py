#!/usr/bin/env python3
"""Additive metadata-only verification of corrected prospective SIMLIB v2."""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
P = ROOT / "runs/research_2026_09_26/raisin_profile_solver_review/csp-filter-interpretation/native-design/timing-followon-design/prospective-engineering"


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def keyvals(path):
    out = {}
    for line in path.read_text().splitlines():
        s = line.split("#", 1)[0].strip()
        if ":" in s:
            k, v = s.split(":", 1)
            out[k.strip()] = v.strip()
    return out


def main():
    old = P / "inputs/cadence.simlib"
    new = P / "inputs-v2/cadence.simlib"
    manifest = P / "inputs-v2/manifest.json"
    meta = json.loads(manifest.read_text())
    for rel, expected in meta["files"].items():
        if sha(ROOT / rel) != expected:
            raise ValueError(("v2 manifest hash mismatch", rel))
    lines = new.read_text().splitlines()
    source = list(csv.DictReader((HERE / "physical-epochs.csv").open()))
    expected = [r["raw_line"] for r in source if r["selected_nominal"] == "True"]
    actual = [line for line in lines if line.startswith("S:")]
    if actual != expected:
        raise ValueError(("v2 exposure byte/order mismatch", len(actual), len(expected)))
    if len(actual) != 117:
        raise ValueError("v2 exposure count")
    exact = {
        "SURVEY:": "SURVEY: DES      FILTERS: grizJH      TELESCOPE: CTIO",
        "PIXSIZE:": "PIXSIZE: 0.27",
        "PSF_UNIT:": "PSF_UNIT: ARCSEC_FWHM",
        "NLIBID:": "NLIBID: 1",
        "LIBID:": "LIBID: 11",
        "RA:": "RA: 54.520306    DECL: -29.391666   NOBS: 117",
        "MWEBV:": "MWEBV: 0.006   PIXSIZE: 0.270",
        "REDSHIFT:": "REDSHIFT: 0.45300000905990601   PEAKMJD: 57707.80078125",
        "FIELD:": "FIELD: DES16E1dcx  # CCDS: [42]",
        "END_LIBID:": "END_LIBID: 11",
        "END_OF_SIMLIB:": "END_OF_SIMLIB: 1 ENTRIES",
    }
    for prefix, full in exact.items():
        found = [line for line in lines if line.startswith(prefix)]
        if found != [full]:
            raise ValueError(("v2 header mismatch", prefix, found))
    branches = {}
    for branch in ("original", "ledger", "noiseless"):
        f = P / f"generation-v2/{branch}/sim.input"
        d = keyvals(f)
        wanted = {"SIMLIB_FILE": "../../inputs-v2/cadence.simlib", "USE_SIMLIB_PEAKMJD": "1",
                  "USE_SIMLIB_REDSHIFT": "1", "GENRANGE_TREST": "-15 45", "GENSIGMA_VPEC": "0",
                  "VPEC_ERR": "0", "VEL_CMBAPEX": "0"}
        if any(d.get(k) != v for k, v in wanted.items()) or any(k.startswith("HOSTLIB_") for k in d):
            raise ValueError(("v2 frame/phase controls changed", branch))
        branches[branch] = {"sha256": sha(f), "relevant": {k: d[k] for k in wanted}}
    result = {"status": "PASS", "v2_manifest_sha256": sha(manifest), "v2_simlib_sha256": sha(new),
              "preserved_malformed_v1_sha256": sha(old), "v2_exposures": len(actual),
              "source_exposure_lines_exact": True, "header_exact": exact, "branches": branches,
              "heliocentric_phase_note": "v11_04d epoch_rest uses zHEL; v2 VEL_CMBAPEX=0 and GENSIGMA_VPEC=0 imply zHEL=zCMB=private header 0.45300000905990601, conditional on source path and no host override"}
    (HERE / "v2-result.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
