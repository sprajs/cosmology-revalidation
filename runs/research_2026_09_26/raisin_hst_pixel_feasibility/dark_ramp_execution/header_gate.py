#!/usr/bin/env python3
"""Validate frozen 14-RAW structure and chronological first-eight metadata only."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from astropy.io import fits

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
DESIGN = json.loads((BASE / "pixel_design/variance-diagnostic-review/calwf3_read_matrix/dark-quality-design/experiment-protocol.json").read_text())
PILOT = BASE / "dark_ramp_raw_pilot/files"
FROZEN_TIMES = DESIGN["time_operator"]["t_seconds_nominal"]


def hash_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def path_for(root: str) -> Path:
    p = HERE / "raw" / f"{root}_raw.fits"
    return p if p.exists() else PILOT / f"{root}_raw.fits"


def main() -> None:
    results = []
    for plan in DESIGN["cohort"]:
        root = plan["root"]
        p = path_for(root)
        row = {"root": root, "visit": plan["visit"], "path": str(p.resolve()), "pass": False}
        results.append(row)
        if not p.exists():
            row["error"] = "missing input"
            continue
        row["bytes"] = p.stat().st_size
        row["sha256"] = hash_file(p)
        try:
            with fits.open(p, memmap=True, do_not_scale_image_data=True) as hdus:
                hdr = hdus[0].header
                if len(hdus) != 81 or hdr["NSAMP"] != 16 or hdr["SAMP_SEQ"] != "SPARS50" or hdr["SUBARRAY"]:
                    raise ValueError("RAW count/pattern/subarray mismatch")
                for key in ("FILTER", "CCDAMP", "CCDGAIN", "EXPFLAG", "BPIXTAB"):
                    got = str(hdr[key]).strip()
                    want = str(plan[key]).strip()
                    if got != want:
                        raise ValueError(f"{key}: {got!r} != {want!r}")
                zero = hdus["SCI", 16].header
                if float(zero["SAMPTIME"]) != 0:
                    raise ValueError("missing zero-read coordinate")
                reads = []
                for index, version in enumerate(range(15, 8, -1)):
                    sci = hdus["SCI", version].header
                    samp = hdus["SAMP", version].header
                    clock = hdus["TIME", version].header
                    dq = hdus["DQ", version].header
                    t = float(sci["SAMPTIME"])
                    if t != FROZEN_TIMES[index] or float(clock["PIXVALUE"]) != t:
                        raise ValueError(f"time mismatch EXTVER {version}")
                    if int(samp["PIXVALUE"]) != index + 2:
                        raise ValueError(f"SAMP count mismatch EXTVER {version}")
                    if (sci["NAXIS1"], sci["NAXIS2"], sci["BITPIX"], sci["BSCALE"], sci["BZERO"], sci["LTV1"], sci["LTV2"]) != (1024, 1024, 16, 1, 32768, 5, 5):
                        raise ValueError(f"SCI encoding mismatch EXTVER {version}")
                    if dq["NAXIS"] == 0:
                        dq_meta = {"mode": "constant", "pixelvalue": int(dq["PIXVALUE"])}
                    elif dq["NAXIS"] == 2 and (dq["NAXIS1"], dq["NAXIS2"]) == (1024, 1024):
                        dq_meta = {"mode": "image", "bitpix": int(dq["BITPIX"])}
                    else:
                        raise ValueError(f"DQ encoding mismatch EXTVER {version}")
                    reads.append({"extver": version, "samptime": t, "dq": dq_meta})
                row["first7_nonzero_reads"] = reads
                row["pass"] = True
        except Exception as exc:
            row["error"] = repr(exc)
    out = {"purpose": "metadata-only; no RAW SCI or DQ array opened",
           "design_sha256": hash_file(BASE / "pixel_design/variance-diagnostic-review/calwf3_read_matrix/dark-quality-design/experiment-protocol.json"),
           "all_pass": all(r["pass"] for r in results), "records": results}
    (HERE / "header-gate.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"all_pass": out["all_pass"], "passed": sum(r["pass"] for r in results), "count": len(results)}, indent=2))


if __name__ == "__main__":
    main()
