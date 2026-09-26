#!/usr/bin/env python3
"""Frozen RAW dark pair extractor. Only --synthetic is allowed before root release."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
from astropy.io import fits

from fixed_mask import R0, R1, R2, footprint

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
DESIGN_PATH = BASE / "pixel_design/variance-diagnostic-review/calwf3_read_matrix/dark-quality-design/experiment-protocol.json"
OPERATORS_PATH = BASE / "pixel_design/variance-diagnostic-review/calwf3_read_matrix/dark-quality-design/fixed-operators.npz"
ACTIVE_START, ACTIVE_STOP = 5, 1019
ROOTS = ("idbx43p7q", "idp247tnq")
CCDTAB_PATH = BASE / "quality_reference_audit/reference_files/t2c16200i_ccd.fits"
QUADRANTS = {
    "B": (slice(5, 512), slice(5, 512)),
    "C": (slice(5, 512), slice(512, 1019)),
    "A": (slice(512, 1019), slice(5, 512)),
    "D": (slice(512, 1019), slice(512, 1019)),
}
MAX_LIVE_ARRAY_BYTES = 268435456
MAX_NEW_OUTPUT_BYTES = 100000000


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def validate_freeze() -> tuple[dict, np.ndarray, np.ndarray, np.ndarray, dict]:
    freeze = json.loads((HERE / "execution-freeze.json").read_text())
    if sha(Path(__file__)) != freeze["executor_sha256"]:
        raise RuntimeError("executor source changed after freeze")
    if sha(DESIGN_PATH) != freeze["design_sha256"]:
        raise RuntimeError("design changed after freeze")
    if sha(OPERATORS_PATH) != freeze["operators_sha256"]:
        raise RuntimeError("operator file changed after freeze")
    if sha(HERE / "fixed-union-bad.npy") != freeze["mask_sha256"]:
        raise RuntimeError("fixed mask changed after freeze")
    if sha(HERE / "header-gate.json") != freeze["header_gate_sha256"]:
        raise RuntimeError("header gate changed after freeze")
    header_gate = json.loads((HERE / "header-gate.json").read_text())
    if not header_gate["all_pass"] or len(header_gate["records"]) != 14:
        raise RuntimeError("14-RAW header gate failed")
    if {x["root"]: x["sha256"] for x in header_gate["records"]} != freeze["raw_sha256"]:
        raise RuntimeError("RAW input hash ledger changed")
    bad = np.load(HERE / "fixed-union-bad.npy", allow_pickle=False)
    if bad.shape != (1024, 1024) or bad.dtype != bool:
        raise RuntimeError("unexpected fixed mask schema")
    with np.load(OPERATORS_PATH, allow_pickle=False) as op:
        times, powers, h = op["times_seconds"].copy(), op["powers"].copy(), op["h"].copy()
    design = json.loads(DESIGN_PATH.read_text())
    if not np.array_equal(times, np.asarray(design["time_operator"]["t_seconds_nominal"])):
        raise RuntimeError("frozen nominal times differ")
    if not np.array_equal(powers, np.asarray(design["time_operator"]["powers"])) or h.shape != (6, 7):
        raise RuntimeError("frozen projection operator differs")
    for row in header_gate["records"]:
        if not Path(row["path"]).exists() or sha(Path(row["path"])) != row["sha256"]:
            raise RuntimeError(f"RAW file changed: {row['root']}")
    return design, bad, powers, h, {x["root"]: x for x in header_gate["records"]}


def check_release() -> None:
    path = HERE / "score-release.json"
    if not path.exists():
        raise SystemExit("No root score-release.json; REAL RAW SCI scoring forbidden")
    release = json.loads(path.read_text())
    freeze_hash = sha(HERE / "execution-freeze.json")
    if release.get("approved") is not True or release.get("execution_freeze_sha256") != freeze_hash:
        raise SystemExit("release does not match exact frozen executor")


def load_ccdtab_gains() -> dict[str, float]:
    with fits.open(CCDTAB_PATH, memmap=True) as hdus:
        rows = [row for row in hdus["CCD"].data
                if str(row["CCDAMP"]).strip() == "ABCD"
                and int(row["CCDCHIP"]) == 1 and float(row["CCDGAIN"]) == 2.5
                and int(row["BINAXIS1"]) == 1 and int(row["BINAXIS2"]) == 1]
        if len(rows) != 1 or int(rows[0]["AMPX"]) != 512 or int(rows[0]["AMPY"]) != 512:
            raise RuntimeError("CCDTAB gain/boundary row is not unique")
        return {amp: float(rows[0][f"ATODGN{amp}"]) for amp in "ABCD"}


def decode_unsigned_sci(hdu) -> np.ndarray:
    hdr = hdu.header
    if (hdr.get("BITPIX"), hdr.get("BZERO"), hdr.get("BSCALE"), hdr.get("NAXIS1"), hdr.get("NAXIS2")) != (16, 32768, 1, 1024, 1024):
        raise RuntimeError("SCI FITS unsigned encoding mismatch")
    raw = np.asarray(hdu.data)
    if raw.dtype.kind != "i" or raw.dtype.itemsize != 2:
        raise RuntimeError("SCI is not raw signed 16-bit storage")
    # Mandatory promotion before subtraction: FITS BZERO=32768 is physical unsigned DN.
    out = raw.astype(np.int32) + np.int32(32768)
    if np.min(out) < 0 or np.max(out) > 65535:
        raise RuntimeError("unsigned decode outside 16-bit range")
    return out


def read_pair_difference(earlier: Path, later: Path) -> tuple[np.ndarray, list[dict], np.ndarray, np.ndarray]:
    diff = np.empty((7, 1024, 1024), dtype=np.float64)
    digital = []
    endpoint_any = np.zeros((1024, 1024), dtype=bool)
    dq_any = np.zeros_like(endpoint_any)
    with fits.open(earlier, memmap=True, do_not_scale_image_data=True) as old, fits.open(later, memmap=True, do_not_scale_image_data=True) as new:
        for k, version in enumerate(range(15, 8, -1)):
            a = decode_unsigned_sci(new["SCI", version])
            b = decode_unsigned_sci(old["SCI", version])
            diff[k] = a.astype(np.float64) - b.astype(np.float64)
            endpoint_any |= (a == 0) | (a == 65535) | (b == 0) | (b == 65535)
            active = np.s_[ACTIVE_START:ACTIVE_STOP, ACTIVE_START:ACTIVE_STOP]
            row = {"read": k + 1, "extver": version,
                   "earlier_zero": int(np.count_nonzero(b[active] == 0)),
                   "earlier_max": int(np.count_nonzero(b[active] == 65535)),
                   "later_zero": int(np.count_nonzero(a[active] == 0)),
                   "later_max": int(np.count_nonzero(a[active] == 65535))}
            for label, hdus in (("earlier", old), ("later", new)):
                dq = hdus["DQ", version]
                if dq.header["NAXIS"] == 0:
                    row[f"{label}_dq_nonzero"] = int(dq.header["PIXVALUE"] != 0) * (1014 * 1014)
                    row[f"{label}_dq_value"] = int(dq.header["PIXVALUE"])
                    if int(dq.header["PIXVALUE"]) != 0:
                        dq_any[:] = True
                else:
                    data = np.asarray(dq.data).astype(np.uint16)
                    row[f"{label}_dq_nonzero"] = int(np.count_nonzero(data[active]))
                    dq_any |= data != 0
            digital.append(row)
    if diff.nbytes + 2 * 1024 * 1024 * (4 + 8) > MAX_LIVE_ARRAY_BYTES:
        raise RuntimeError("live array budget exceeded")
    return diff, digital, endpoint_any, dq_any


def matrix_score(data: np.ndarray, selected: np.ndarray, h: np.ndarray) -> dict:
    # data shape (7, ny, nx); selected is fixed mask over the same rectangle.
    vectors = data[:, selected]
    n = int(vectors.shape[1])
    if n == 0:
        return {"n": 0, "mean": [None] * 7, "gamma": None,
                "direct": [None] * len(h), "matrix": [None] * len(h)}
    gamma = (vectors @ vectors.T) / (2.0 * n)
    projections = h @ vectors
    direct = np.mean(np.square(projections), axis=1) / 2.0
    matrix = np.einsum("pi,ij,pj->p", h, gamma, h)
    gap = np.abs(direct - matrix)
    if np.any(gap > np.maximum(1e-12, 1e-10 * np.maximum(np.abs(direct), np.abs(matrix)))):
        raise RuntimeError(f"direct/matrix contraction mismatch: {float(np.max(gap))}")
    return {"n": n, "mean": np.mean(vectors, axis=1).tolist(),
            "gamma": gamma, "direct": direct, "matrix": matrix, "max_gap": float(np.max(gap))}


def spatial_weights(bad: np.ndarray, coverage: list[dict]) -> list[dict]:
    extent = math.ceil(R2 + 0.5)
    side = extent * 2 + 1
    ap = np.zeros((side, side), dtype=np.float64)
    bg = np.zeros_like(ap)
    small, inner, outer = (footprint(r) for r in (R0, R1, R2))
    for source, target in ((small, ap), (outer, bg)):
        start = (side - source.shape[0]) // 2
        target[start:start + source.shape[0], start:start + source.shape[1]] += source
    start = (side - inner.shape[0]) // 2
    bg[start:start + inner.shape[0], start:start + inner.shape[1]] -= inner
    weights = []
    for row in coverage:
        cy, cx = row["y"], row["x"]
        good = ~bad[cy - extent:cy + extent + 1, cx - extent:cx + extent + 1]
        frac_a = float(np.sum(good * ap) / np.sum(ap))
        frac_b = float(np.sum(good * bg) / np.sum(bg))
        if abs(frac_a - row["aperture_coverage"]) > 1e-12 or abs(frac_b - row["annulus_coverage"]) > 1e-12:
            raise RuntimeError("spatial coverage changed after freeze")
        if row["eligible"]:
            w = good * ap - (np.sum(good * ap) / np.sum(good * bg)) * good * bg
            if abs(float(np.sum(w))) > 1e-8:
                raise RuntimeError("spatial background balance failed")
            weights.append({**row, "extent": extent, "w": w})
    return weights


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def synthetic() -> None:
    freeze = json.loads((HERE / "execution-freeze.json").read_text())
    if sha(Path(__file__)) != freeze["executor_sha256"]:
        raise RuntimeError("source changed after freeze")
    with np.load(OPERATORS_PATH, allow_pickle=False) as op:
        h = op["h"]
    # Integer slopes, common zero and single cosmic-ray step retained.
    a = np.array([[[100 + 2 * t, 400 + 4 * t], [200 + t, 300 + 3 * t]] for t in range(7)], dtype=np.float64)
    b = a + np.array([[[0, 0], [0, 0]] if t < 3 else [[0, 11], [0, 0]] for t in range(7)], dtype=np.float64)
    result = matrix_score(b - a, np.ones((2, 2), dtype=bool), h)
    assert result["n"] == 4 and result["max_gap"] < 1e-12
    # Unsigned-before-cast check at the FITS wrap boundary, plus constant-offset cancellation.
    physical = np.array([0, 65535], dtype=np.int32)
    stored = (physical - 32768).astype(np.int16)
    decoded = stored.astype(np.int32) + 32768
    assert np.array_equal(decoded, physical)
    assert float(decoded[1]) - float(decoded[0]) == 65535.0
    gains = load_ccdtab_gains()
    assert list(gains) == list("ABCD") and all(2.0 < value < 3.0 for value in gains.values())
    report = {"pass": True, "source_sha256": sha(Path(__file__)),
              "design_sha256": sha(DESIGN_PATH), "operator_sha256": sha(OPERATORS_PATH),
              "matrix_max_gap": result["max_gap"], "unsigned_boundary": decoded.tolist(),
              "source_ccdtab_gain": gains,
              "real_sci_read": False}
    (HERE / "executor-synthetic.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


def score() -> None:
    started = time.monotonic()
    check_release()
    design, bad, powers, h, sources = validate_freeze()
    fixed = json.loads((HERE / "fixed-mask-result.json").read_text())
    if fixed["fixed_union_mask_sha256"] != sha(HERE / "fixed-union-bad.npy"):
        raise RuntimeError("mask result mismatch")
    sites = spatial_weights(bad, fixed["coverage"])
    gains = load_ccdtab_gains()
    n_pair = len(design["pairs"])
    q_gamma = np.full((n_pair, 2, 4, 7, 7), np.nan)
    q_mean = np.full((n_pair, 2, 4, 7), np.nan)
    block_gamma = np.full((n_pair, 16, 16, 7, 7), np.nan)
    block_mean = np.full((n_pair, 16, 16, 7), np.nan)
    quadrants = list(QUADRANTS)
    quadrant_rows, block_rows, aperture_rows, digital_rows, contributions, gates = [], [], [], [], [], []
    for ipair, pair in enumerate(design["pairs"]):
        if time.monotonic() - started > design["resources_analysis"]["full_analysis_seconds"]:
            raise RuntimeError("cumulative analysis-time cap")
        old = Path(sources[pair["earlier"]]["path"])
        new = Path(sources[pair["later"]]["path"])
        D, dig, endpoint_any, dq_any = read_pair_difference(old, new)
        digital_rows.extend({"pair": pair["pair"], **row} for row in dig)
        primary = np.einsum("t,tyx->yx", h[0], D)
        active_good = np.zeros((1024, 1024), dtype=bool)
        active_good[ACTIVE_START:ACTIVE_STOP, ACTIVE_START:ACTIVE_STOP] = True
        active_good &= ~bad
        total = float(np.sum(np.square(primary[active_good])) / 2)
        for label, condition in (("digital_endpoint", endpoint_any), ("raw_DQ", dq_any)):
            selection = active_good & condition
            subtotal = float(np.sum(np.square(primary[selection])) / 2)
            contributions.append({"pair": pair["pair"], "kind": pair["kind"],
                                  "category": label, "selected_pixels": int(np.count_nonzero(selection)),
                                  "total_masked_pixels": int(np.count_nonzero(active_good)),
                                  "half_squared_slope_contribution": subtotal,
                                  "half_squared_slope_all_masked": total,
                                  "fraction_of_total": None if total == 0 else subtotal / total})
        del primary
        for iq, (quad, (ys, xs)) in enumerate(QUADRANTS.items()):
            slab = D[:, ys, xs]
            for iscope, (scope, mask) in enumerate((("masked", ~bad[ys, xs]), ("geometry_only", np.ones(slab.shape[1:], dtype=bool)))):
                result = matrix_score(slab, mask, h)
                q_gamma[ipair, iscope, iq] = result["gamma"]
                q_mean[ipair, iscope, iq] = result["mean"]
                gates.append({"pair": pair["pair"], "region": quad, "scope": scope, "max_gap": result["max_gap"]})
                for power, direct, matrix in zip(powers, result["direct"], result["matrix"]):
                    quadrant_rows.append({"pair": pair["pair"], "visit": pair["visit"], "kind": pair["kind"],
                                          "quadrant": quad, "scope": scope, "n": result["n"],
                                          "power": float(power), "dn2_per_nominal_s2": float(direct),
                                          "matrix_dn2_per_nominal_s2": float(matrix),
                                          "electron2_per_nominal_s2": float(direct) * gains[quad] ** 2})
        for by in range(16):
            ys = slice(max(5, by * 64), min(1019, (by + 1) * 64))
            for bx in range(16):
                xs = slice(max(5, bx * 64), min(1019, (bx + 1) * 64))
                result = matrix_score(D[:, ys, xs], ~bad[ys, xs], h)
                if result["n"]:
                    block_gamma[ipair, by, bx] = result["gamma"]
                    block_mean[ipair, by, bx] = result["mean"]
                block_rows.append({"pair": pair["pair"], "visit": pair["visit"],
                                   "block_x": bx, "block_y": by, "n": result["n"],
                                   "primary_dn2_per_nominal_s2": None if not result["n"] else float(result["direct"][0]),
                                   "max_gate_gap": result.get("max_gap")})
        for site in sites:
            cx, cy, ex = site["x"], site["y"], site["extent"]
            z = np.einsum("tyx,yx->t", D[:, cy - ex:cy + ex + 1, cx - ex:cx + ex + 1], site["w"])
            slopes = h @ z
            for power, slope in zip(powers, slopes):
                aperture_rows.append({"pair": pair["pair"], "visit": pair["visit"], "kind": pair["kind"],
                                      "x": cx, "y": cy, "power": float(power),
                                      "signed_dn_per_nominal_s": float(slope),
                                      "half_squared_dn2_per_nominal_s2": float(slope * slope / 2),
                                      "aperture_coverage": site["aperture_coverage"],
                                      "annulus_coverage": site["annulus_coverage"]})
        del D, endpoint_any, dq_any
        # Preserve every completed pair if the cumulative bound or a later
        # input fails. These files contain observed outcomes only after release.
        np.savez_compressed(HERE / "partial-read-matrix.npz",
                            quadrant_gamma=q_gamma[:ipair + 1], quadrant_mean=q_mean[:ipair + 1],
                            block_gamma=block_gamma[:ipair + 1], block_mean=block_mean[:ipair + 1],
                            quadrant_names=np.array(quadrants), powers=powers)
        (HERE / "partial-score.json").write_text(json.dumps({
            "completed_pairs": [x["pair"] for x in design["pairs"][:ipair + 1]],
            "elapsed_seconds": time.monotonic() - started,
            "quadrant_rows": quadrant_rows, "block_rows": block_rows,
            "aperture_rows": aperture_rows, "digital_rows": digital_rows,
            "digital_contributions": contributions, "arithmetic_gates": gates,
        }) + "\n")
        if ((HERE / "partial-score.json").stat().st_size
                + (HERE / "partial-read-matrix.npz").stat().st_size > MAX_NEW_OUTPUT_BYTES):
            raise RuntimeError("partial output size cap exceeded")
        if time.monotonic() - started > design["resources_analysis"]["full_analysis_seconds"]:
            raise RuntimeError("cumulative analysis-time cap after completed pair")
    np.savez_compressed(HERE / "read-matrix.npz", quadrant_gamma=q_gamma, quadrant_mean=q_mean,
                        block_gamma=block_gamma, block_mean=block_mean,
                        quadrant_names=np.array(quadrants), powers=powers)
    write_csv(HERE / "pair-quadrant-power.csv", quadrant_rows)
    write_csv(HERE / "normal-search-sensitivity.csv",
              [row for row in quadrant_rows if row["kind"] == "prespecified_sensitivity_one_pair"])
    write_csv(HERE / "block-pair-power.csv", block_rows)
    write_csv(HERE / "aperture-pair-power.csv", aperture_rows)
    write_csv(HERE / "raw-dq-digital-ledger.csv", digital_rows)
    write_csv(HERE / "raw-dq-digital-contributions.csv", contributions)
    summary = {"elapsed_seconds": time.monotonic() - started, "pairs": n_pair,
               "site_support": len(sites), "gates": gates,
               "max_direct_matrix_gap": max(x["max_gap"] for x in gates),
               "source_sha256": sha(Path(__file__)), "score_release_sha256": sha(HERE / "score-release.json")}
    (HERE / "arithmetic-gates.json").write_text(json.dumps(summary, indent=2) + "\n")
    output_bytes = sum((HERE / name).stat().st_size for name in
                       ("read-matrix.npz", "pair-quadrant-power.csv", "normal-search-sensitivity.csv",
                        "block-pair-power.csv", "aperture-pair-power.csv", "raw-dq-digital-ledger.csv",
                        "raw-dq-digital-contributions.csv", "arithmetic-gates.json",
                        "partial-read-matrix.npz", "partial-score.json"))
    if output_bytes > MAX_NEW_OUTPUT_BYTES:
        raise RuntimeError("new output size cap exceeded")
    print(json.dumps({"status": "scored", "pairs": n_pair, "seconds": summary["elapsed_seconds"],
                      "output_bytes": output_bytes}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("synthetic", "score"))
    args = parser.parse_args()
    if args.mode == "synthetic":
        synthetic()
    else:
        score()


if __name__ == "__main__":
    main()
