"""Rebuild frozen HST geometry and masks from the eight calibrated exposures."""

from __future__ import annotations
import json, math, time, zipfile
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import convolve
from lib.hst_geometry import curvature_gate, make_operator, tangent, untangent
from lib.paths import DATA
from lib.records import sha256 as digest, write_json

BASE = DATA / "hst/exposures"
PAM = DATA / "hst/pam.fits"
RA0, DEC0 = 8.9685, -43.35812
NPIX = 1014
OPERATOR_FIELDS_CAP_BYTES = 50_000_000
DEFAULTS = {}
HERE = None


def verify_freeze():
    return json.loads((HERE / "input-freeze.json").read_text())


def read_design():
    return (
        {"lattice_spacing_arcsec": 6},
        {
            "secondary": {
                "aperture_min_PAMweighted_coverage": 0.9,
                "annulus_min_native_area_coverage": 0.75,
            }
        },
    )


def exposure_metadata():
    out = []
    for group, code in [("search", "icxoi1"), ("template", "icxoi4")]:
        candidates = sorted(BASE.glob(code + "*_flt.fits"))
        assert len(candidates) == 4
        one = []
        for path in candidates:
            with fits.open(
                path, memmap=True, lazy_load_hdus=True, do_not_scale_image_data=True
            ) as hdus:
                h = hdus[0].header
                s = hdus["SCI"].header
                assert h["FILTER"] == "F160W" and h["DETECTOR"] == "IR"
                assert h["SUBARRAY"] is False and h["SUBTYPE"] == "FULLIMAG"
                assert [s["NAXIS1"], s["NAXIS2"]] == [NPIX, NPIX]
                assert [s["LTV1"], s["LTV2"]] == [0.0, 0.0]
                assert [s["LTM1_1"], s["LTM2_2"]] == [1.0, 1.0]
                assert h["IDCTAB"] == "iref$w3m18525i_idc.fits"
                one.append(
                    {
                        "path": str(path),
                        "filename": path.name,
                        "sha256": digest(path),
                        "visit": group,
                        "expstart": float(h["EXPSTART"]),
                        "photflam": float(h["PHOTFLAM"]),
                        "wcs": WCS(s, naxis=2),
                    }
                )
        one.sort(key=lambda r: (r["expstart"], r["filename"]))
        for ordinal, r in enumerate(one, 1):
            r["ordinal"] = ordinal
            r["design"] = ordinal in (1, 3)
        out += one
    return out


def tangent_bounds_and_grid(exps, spacing):
    # WCS-only bounds; candidate annulus ring is checked against every FITS
    # footprint with a two-native-pixel edge guard before any data values.
    bounds = []
    side = np.linspace(3, 1012, 33)
    edge = np.concatenate(
        (
            np.column_stack((side, np.full_like(side, 3))),
            np.column_stack((side, np.full_like(side, 1012))),
            np.column_stack((np.full_like(side, 3), side)),
            np.column_stack((np.full_like(side, 1012), side)),
        )
    )
    for e in exps:
        rd = e["wcs"].all_pix2world(edge, 1)
        x, y = tangent(rd[:, 0], rd[:, 1], RA0, DEC0)
        bounds.append((min(x), max(x), min(y), max(y)))
    xlo, xhi = max(q[0] for q in bounds), min(q[1] for q in bounds)
    ylo, yhi = max(q[2] for q in bounds), min(q[3] for q in bounds)
    if xhi <= xlo or yhi <= ylo:
        raise ValueError("no common WCS bounding rectangle")
    grid = []
    for iy in range(math.ceil(ylo / spacing), math.floor(yhi / spacing) + 1):
        for ix in range(math.ceil(xlo / spacing), math.floor(xhi / spacing) + 1):
            x, y = ix * spacing, iy * spacing
            if math.hypot(x, y) <= 5:
                continue
            ra, dec = untangent(x, y, RA0, DEC0)
            theta = np.linspace(0, 2 * math.pi, 64, endpoint=False)
            xr = x + 2 * np.cos(theta)
            yr = y + 2 * np.sin(theta)
            rr, dd = untangent(xr, yr, RA0, DEC0)
            ringworld = np.column_stack((rr, dd))
            inside = True
            for e in exps:
                pix = e["wcs"].all_world2pix(ringworld, 1)
                if not np.all(np.isfinite(pix)) or not np.all(
                    (pix >= 3) & (pix <= 1012)
                ):
                    inside = False
                    break
            if inside:
                tx = min(3, int(4 * (x - xlo) / (xhi - xlo)))
                ty = min(3, int(4 * (y - ylo) / (yhi - ylo)))
                grid.append(
                    {
                        "id": len(grid),
                        "east": x,
                        "north": y,
                        "ra": float(ra),
                        "dec": float(dec),
                        "tile": 4 * ty + tx,
                    }
                )
    return (xlo, xhi, ylo, yhi), grid


def make_source_mask(path, valid_policy="all9"):
    """Design SCI only. Invalid 3x3 boxes are not source detections."""
    with fits.open(path, memmap=True, do_not_scale_image_data=True) as h:
        sci, err, dq = h["SCI"].data, h["ERR"].data, h["DQ"].data
        valid = np.isfinite(sci) & np.isfinite(err) & (err > 0) & (dq == 0)
        background = np.zeros((NPIX, NPIX), np.float32)
        for y in range(0, NPIX, 64):
            for x in range(0, NPIX, 64):
                sub = sci[y : y + 64, x : x + 64]
                vv = valid[y : y + 64, x : x + 64]
                med = float(np.median(sub[vv])) if vv.any() else 0.0
                background[y : y + 64, x : x + 64] = med
        residual = np.where(valid, sci - background, 0.0).astype(np.float32)
        err2 = np.where(valid, err * err, 0.0).astype(np.float32)
        kernel = np.ones((3, 3), np.float32)
        sums = convolve(residual, kernel, mode="constant", cval=0.0)
        variances = convolve(err2, kernel, mode="constant", cval=0.0)
        good9 = (
            convolve(
                valid.astype(np.uint8), kernel.astype(np.uint8), mode="constant", cval=0
            )
            == 9
        )
        mask = good9 & (np.abs(sums) > 5 * np.sqrt(variances))
        return (
            mask,
            background,
            {
                "valid_pixels": int(valid.sum()),
                "source_pixels": int(mask.sum()),
                "unscored_3x3_boxes_invalid_quality_or_edge": int((~good9).sum()),
            },
        )


def near_source(wcs, mask, candidate):
    center = wcs.all_world2pix([[candidate["ra"], candidate["dec"]]], 1)[0]
    cx, cy = center
    rad = math.ceil(2.4 / 0.11) + 4
    xmin, xmax = max(1, int(cx) - rad), min(NPIX, int(cx) + rad)
    ymin, ymax = max(1, int(cy) - rad), min(NPIX, int(cy) + rad)
    iy, ix = np.nonzero(mask[ymin - 1 : ymax, xmin - 1 : xmax])
    if len(ix) == 0:
        return False
    pix = np.column_stack((ix + xmin, iy + ymin))
    corners = []
    for dx, dy in [(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)]:
        rd = wcs.all_pix2world(pix + [dx, dy], 1)
        east, north = tangent(rd[:, 0], rd[:, 1], RA0, DEC0)
        corners.append(
            np.column_stack((east - candidate["east"], north - candidate["north"]))
        )
    polygon = np.stack(corners, axis=1)
    # Distance from the lattice center to the actual WCS-mapped flagged
    # pixel quadrilateral, rather than merely to its center.
    end = np.roll(polygon, -1, axis=1)
    vec = end - polygon
    t = np.clip(-np.sum(polygon * vec, axis=2) / np.sum(vec * vec, axis=2), 0, 1)
    closest = polygon + t[:, :, None] * vec
    d = np.sqrt(np.sum(closest * closest, axis=2)).min(axis=1)
    cross = polygon[:, :, 0] * vec[:, :, 1] - polygon[:, :, 1] * vec[:, :, 0]
    inside = np.all(cross >= 0, axis=1) | np.all(cross <= 0, axis=1)
    d[inside] = 0
    return bool(np.any(d <= 2.4))


def operator_validity(path, op):
    with fits.open(path, memmap=True, do_not_scale_image_data=True) as h:
        flat = op.flat
        sci, err, dq = (h[n].data.ravel()[flat] for n in ["SCI", "ERR", "DQ"])
        v = np.isfinite(sci) & np.isfinite(err) & (err > 0) & (dq == 0)
        bits = {}
        for d in dq[~v]:
            bits[str(int(d))] = bits.get(str(int(d)), 0) + 1
        # FITS BITPIX16 is signed here. Cast through uint16 to decode its
        # actual 16 stored bits, without sign-extending flagged bit 15.
        dqbits = dq.astype(np.uint16).astype(np.uint64)
        bit_counts = {
            str(1 << bit): int(np.count_nonzero((dqbits & (1 << bit)) != 0))
            for bit in range(16)
            if np.any((dqbits & (1 << bit)) != 0)
        }
        full_a = float(np.dot(op.aperture, op.pam))
        good_a = float(np.dot(op.aperture * op.pam, v))
        full_b = float(np.sum(op.annulus))
        good_b = float(np.dot(op.annulus, v))
        return v, {
            "strict": bool(np.all(v)),
            "aperture_coverage": good_a / full_a if full_a > 0 else -1,
            "annulus_coverage": good_b / full_b if full_b > 0 else -1,
            "dq_values_at_invalid": bits,
            "dq_bit_counts": bit_counts,
            "nonfinite_sci": int(np.count_nonzero(~np.isfinite(sci))),
            "nonfinite_err": int(np.count_nonzero(~np.isfinite(err))),
            "nonpositive_err": int(np.count_nonzero(np.isfinite(err) & (err <= 0))),
            "n_invalid": int((~v).sum()),
        }


def synthetic_operator_gate(op, v, kappa):
    a, b, p = op.aperture, op.annulus, op.pam
    vv = v.astype(float)
    A, B = float(np.dot(vv * a, p)), float(np.dot(vv, b))
    if A <= 0 or B <= 0:
        return False
    w = kappa * (vv * a * p - (A / B) * vv * b)
    if abs(w.sum()) > 1e-10 * max(1, np.sum(abs(w))):
        return False
    x = 2 + 0.0001 * (op.flat % NPIX) + 0.0002 * (op.flat // NPIX)
    err = 1 + 0.0001 * (op.flat % NPIX)
    direct = kappa * (np.dot(vv * a * p, x) - A * np.dot(vv * b, x) / B)
    variance = kappa**2 * np.dot((vv * a * p - A / B * vv * b) ** 2, err**2)
    inject = np.zeros(len(w))
    j = int(np.argmax(abs(w)))
    inject[j] = -0.375  # signed known test pulse, independent of FITS SCI
    injection_change = np.dot(w, x + inject) - np.dot(w, x)
    return (
        abs(np.dot(w, x) - direct) <= 1e-10 * max(1, abs(direct))
        and abs(np.dot(w * w, err * err) - variance) <= 1e-10 * max(1, variance)
        and abs(injection_change - w[j] * inject[j])
        <= 1e-10 * max(1, abs(injection_change))
        and abs(np.dot(2 * w, x) - 2 * np.dot(w, x))
        <= 1e-10 * max(1, abs(np.dot(w, x)))
        and abs(kappa * np.dot(vv * a * p, np.ones(len(w))) - kappa * A) < 1e-9
    )


def stage_a():
    start = time.monotonic()
    frozen = verify_freeze()
    p, a = read_design()
    exps = exposure_metadata()
    with fits.open(PAM, memmap=True) as pamfile:
        pam = pamfile[1].data
        assert pam.shape == (NPIX, NPIX)
        grid11 = [50, 141, 233, 324, 416, 507, 598, 690, 781, 873, 964]
        curvature = [curvature_gate(e["wcs"], RA0, DEC0, grid11) for e in exps]
        if max(curvature) > 1e-4:
            raise ValueError("corner-bilinear WCS curvature gate failed")
        bounds, grid = tangent_bounds_and_grid(exps, p["lattice_spacing_arcsec"])
        masks, medians, maskmeta = {}, {}, {}
        for i, e in enumerate(exps):
            if e["design"]:
                masks[i], medians[i], maskmeta[e["filename"]] = make_source_mask(
                    e["path"]
                )
        np.savez_compressed(
            HERE / "source_masks.npz", **{f"exposure_{i}": z for i, z in masks.items()}
        )
        rows, op_rows = [], []
        flats, aps, anns, ps = [], [], [], []
        total = 0

        def save_arrays(path):
            # Write one concatenated field at a time. np.savez_compressed
            # would materialize all four concatenations simultaneously.
            with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as z:
                for name, parts, dtype in [
                    ("flat", flats, np.int32),
                    ("aperture", aps, np.float32),
                    ("annulus", anns, np.float32),
                    ("pam", ps, np.float32),
                ]:
                    field = np.concatenate(parts) if parts else np.array([], dtype)
                    with z.open(name + ".npy", "w", force_zip64=True) as dst:
                        np.lib.format.write_array(dst, field, allow_pickle=False)
                    del field

        def stop_stage_a(reason, current=None):
            if current is not None:
                current["reason"] = reason
                rows.append(current)
            save_arrays(HERE / "operators.partial.npz")
            partial = {
                "status": "stage_a_stopped_before_heldout_aperture_sums",
                "completed_candidates": len(rows),
                "total_candidates": len(grid),
                "reason": reason,
                "rows": rows,
                "operator_rows": op_rows,
                "retained_operator_field_bytes": total * 16,
                "source_masks_sha256": digest(HERE / "source_masks.npz"),
                "operators_partial_sha256": digest(HERE / "operators.partial.npz"),
                "mask_metadata": maskmeta,
            }
            (HERE / "stage_a_partial.json").write_text(
                json.dumps(partial, indent=2) + "\n"
            )
            raise RuntimeError(
                "Stage A stopped at a numerical or resource gate: " + reason
            )

        for c in grid:
            reason = None
            if any(near_source(exps[i]["wcs"], masks[i], c) for i in masks):
                reason = "design_source_mask"
            row = dict(c, reason=reason)
            if reason is None:
                strict, secondary = True, True
                qualities = []
                temp = []
                for i, e in enumerate(exps):
                    try:
                        op = make_operator(
                            e["wcs"],
                            c["east"],
                            c["north"],
                            RA0,
                            DEC0,
                            pam,
                            0.4,
                            (1.2, 2.0),
                        )
                        v, q = operator_validity(e["path"], op)
                    except (ValueError, AssertionError, FloatingPointError) as exc:
                        stop_stage_a(
                            "operator_or_quality_numeric:"
                            + repr(exc)
                            + ":candidate="
                            + str(c["id"])
                            + ":exposure="
                            + e["filename"],
                            row,
                        )
                    q["resolution"] = op.resolution
                    q["numerical_closure"] = op.closure
                    q["filename"] = e["filename"]
                    if not synthetic_operator_gate(
                        op, np.ones(len(v), bool), e["photflam"] / exps[0]["photflam"]
                    ):
                        stop_stage_a(
                            "synthetic_full_operator_failure:candidate="
                            + str(c["id"])
                            + ":exposure="
                            + e["filename"],
                            row,
                        )
                    if (
                        q["aperture_coverage"] > 0
                        and q["annulus_coverage"] > 0
                        and not synthetic_operator_gate(
                            op, v, e["photflam"] / exps[0]["photflam"]
                        )
                    ):
                        stop_stage_a(
                            "synthetic_masked_operator_failure:candidate="
                            + str(c["id"])
                            + ":exposure="
                            + e["filename"],
                            row,
                        )
                    strict &= q["strict"]
                    secondary &= (
                        q["aperture_coverage"] >= 0.90 and q["annulus_coverage"] >= 0.75
                    )
                    qualities.append(q)
                    temp.append((i, op))
                if reason is None:
                    row["strict"] = bool(strict)
                    row["secondary"] = bool(secondary)
                    row["qualities"] = qualities
                    # Design-only background stratum; never read test SCI here.
                    bg = []
                    for i in masks:
                        center = exps[i]["wcs"].all_world2pix([[c["ra"], c["dec"]]], 1)[
                            0
                        ]
                        x, y = [int(round(z)) for z in center]
                        bg.append(float(medians[i][y - 1, x - 1]))
                    row["design_background_mean"] = float(np.mean(bg))
                    if strict or secondary:
                        proposed = total + sum(len(op.flat) for _, op in temp)
                        if proposed * 16 > OPERATOR_FIELDS_CAP_BYTES:
                            stop_stage_a(
                                "100MB_working_array_cap:next_candidate="
                                + str(c["id"])
                                + ":proposed_operator_bytes="
                                + str(proposed * 16),
                                row,
                            )
                        for i, op in temp:
                            n = len(op.flat)
                            op_rows.append(
                                {
                                    "candidate": c["id"],
                                    "exposure": i,
                                    "start": total,
                                    "end": total + n,
                                }
                            )
                            flats.append(op.flat)
                            aps.append(op.aperture.astype(np.float32))
                            anns.append(op.annulus.astype(np.float32))
                            ps.append(op.pam.astype(np.float32))
                            total += n
            rows.append(row)
            if time.monotonic() - start > 120:
                stop_stage_a("120-second_stage_a_benchmark_cap")
        # The full design arrays are no longer needed while concatenating.
        masks.clear()
        medians.clear()
        save_arrays(HERE / "operators.npz")
    strict = [r for r in rows if r.get("strict")]
    secondary = [r for r in rows if r.get("secondary")]
    out = {
        "status": "stage_a_complete_no_heldout_aperture_sums",
        "elapsed_seconds": time.monotonic() - start,
        "input_freeze_sha256": digest(HERE / "input-freeze.json"),
        "bounds_east_north_arcsec": bounds,
        "curvature_max_native_pixels": max(curvature),
        "curvature_by_exposure": curvature,
        "mask_metadata": maskmeta,
        "n_geometry_grid": len(grid),
        "n_strict": len(strict),
        "n_secondary": len(secondary),
        "strict_tiles": sorted(set(r["tile"] for r in strict)),
        "secondary_tiles": sorted(set(r["tile"] for r in secondary)),
        "rows": rows,
        "operator_rows": op_rows,
        "exposures": [{k: v for k, v in e.items() if k != "wcs"} for e in exps],
        "source_masks_sha256": digest(HERE / "source_masks.npz"),
        "operators_sha256": digest(HERE / "operators.npz"),
    }
    (HERE / "stage_a.json").write_text(json.dumps(out, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: out[k]
                for k in [
                    "status",
                    "elapsed_seconds",
                    "n_geometry_grid",
                    "n_strict",
                    "n_secondary",
                    "strict_tiles",
                    "secondary_tiles",
                ]
            },
            indent=2,
        )
    )


def run(out, config):
    global HERE
    HERE = out
    freeze = {
        "files": [
            {"filename": p.name, "sha256": digest(p)}
            for p in sorted(BASE.glob("*_flt.fits"))
        ],
        "pam_sha256": digest(PAM),
    }
    write_json(out / "input-freeze.json", freeze)
    stage_a()
    result = json.loads((out / "stage_a.json").read_text())
    original = json.loads((DATA / "hst/stage_a.json").read_text())
    identity = lambda s: [
        (r["id"], r.get("strict", False), r.get("secondary", False)) for r in s["rows"]
    ]
    if identity(result) != identity(original):
        raise RuntimeError("Rebuilt geometry support differs from frozen design")
    with (
        np.load(out / "operators.npz") as rebuilt,
        np.load(DATA / "hst/operators.npz") as frozen,
    ):
        exact = {k: np.array_equal(rebuilt[k], frozen[k]) for k in rebuilt.files}
    if not all(exact.values()):
        raise RuntimeError("Rebuilt aperture operators differ from frozen arrays")
    return {
        "geometry_candidates": result["n_geometry_grid"],
        "strict_support": result["n_strict"],
        "secondary_support": result["n_secondary"],
        "exact_operator_arrays": exact,
        "scope": "Rebuild of the fixed geometry, design-exposure source masks and quality support. Held-out aperture brightness sums are evaluated separately by hst-repeat.",
    }
