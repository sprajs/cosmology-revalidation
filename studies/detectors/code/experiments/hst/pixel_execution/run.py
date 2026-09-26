"""Two-stage current-product WFC3/IR signed repeat-noise experiment.

Stage A creates geometry, design masks and support only. Stage B requires a
separate root release and is the first code path that sums held-out SCI values.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import time
import zipfile
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import convolve

from geometry import curvature_gate, make_operator, tangent, untangent

HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "pixel_acquisition"
DESIGN = HERE.parent / "pixel_design"
PAM = BASE / "pam_validation" / "ir_wfc3_map.fits"
RA0, DEC0 = 8.9685, -43.35812
NPIX = 1014
# New in-memory arrays: keep retained operator fields below 50 MB so the
# design masks/medians, one concatenated field, and one bounded 128x128
# quadrature chunk remain below the declared 100 MB working-array cap.
OPERATOR_FIELDS_CAP_BYTES = 50_000_000


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for z in iter(lambda: f.read(1 << 20), b""):
            h.update(z)
    return h.hexdigest()


def verify_freeze():
    frozen = json.loads((HERE / "input-freeze.json").read_text())
    for record in frozen["files"]:
        path = Path(record["path"])
        if digest(path) != record["sha256"] or path.stat().st_size != record["bytes"]:
            raise ValueError("input/source freeze mismatch: " + str(path))
    return frozen


def read_design():
    p = json.loads((DESIGN / "protocol.json").read_text())
    a = json.loads((DESIGN / "masked-secondary-amendment.json").read_text())
    assert p["design_exposure_ordinals"] == [1, 3]
    assert p["measurement_exposure_ordinals"] == [2, 4]
    assert a["secondary"]["aperture_min_PAMweighted_coverage"] == .9
    assert a["secondary"]["annulus_min_native_area_coverage"] == .75
    return p, a


def exposure_metadata():
    out = []
    for group, code in [("search", "icxoi1"), ("template", "icxoi4")]:
        candidates = sorted(BASE.glob(code + "*_flt.fits"))
        assert len(candidates) == 4
        one = []
        for path in candidates:
            with fits.open(path, memmap=True, lazy_load_hdus=True,
                           do_not_scale_image_data=True) as hdus:
                h = hdus[0].header
                s = hdus["SCI"].header
                assert h["FILTER"] == "F160W" and h["DETECTOR"] == "IR"
                assert h["SUBARRAY"] is False and h["SUBTYPE"] == "FULLIMAG"
                assert [s["NAXIS1"], s["NAXIS2"]] == [NPIX, NPIX]
                assert [s["LTV1"], s["LTV2"]] == [0., 0.]
                assert [s["LTM1_1"], s["LTM2_2"]] == [1., 1.]
                assert h["IDCTAB"] == "iref$w3m18525i_idc.fits"
                one.append({"path": str(path), "filename": path.name, "sha256": digest(path),
                            "visit": group, "expstart": float(h["EXPSTART"]),
                            "photflam": float(h["PHOTFLAM"]), "wcs": WCS(s, naxis=2)})
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
    edge = np.concatenate((np.column_stack((side, np.full_like(side, 3))),
                           np.column_stack((side, np.full_like(side, 1012))),
                           np.column_stack((np.full_like(side, 3), side)),
                           np.column_stack((np.full_like(side, 1012), side))))
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
                if not np.all(np.isfinite(pix)) or not np.all((pix >= 3) & (pix <= 1012)):
                    inside = False
                    break
            if inside:
                tx = min(3, int(4 * (x - xlo) / (xhi - xlo)))
                ty = min(3, int(4 * (y - ylo) / (yhi - ylo)))
                grid.append({"id": len(grid), "east": x, "north": y,
                             "ra": float(ra), "dec": float(dec), "tile": 4 * ty + tx})
    return (xlo, xhi, ylo, yhi), grid


def make_source_mask(path, valid_policy="all9"):
    """Design SCI only. Invalid 3x3 boxes are not source detections."""
    with fits.open(path, memmap=True, do_not_scale_image_data=True) as h:
        sci, err, dq = h["SCI"].data, h["ERR"].data, h["DQ"].data
        valid = np.isfinite(sci) & np.isfinite(err) & (err > 0) & (dq == 0)
        background = np.zeros((NPIX, NPIX), np.float32)
        for y in range(0, NPIX, 64):
            for x in range(0, NPIX, 64):
                sub = sci[y:y+64, x:x+64]
                vv = valid[y:y+64, x:x+64]
                med = float(np.median(sub[vv])) if vv.any() else 0.
                background[y:y+64, x:x+64] = med
        residual = np.where(valid, sci - background, 0.).astype(np.float32)
        err2 = np.where(valid, err * err, 0.).astype(np.float32)
        kernel = np.ones((3, 3), np.float32)
        sums = convolve(residual, kernel, mode="constant", cval=0.)
        variances = convolve(err2, kernel, mode="constant", cval=0.)
        good9 = convolve(valid.astype(np.uint8), kernel.astype(np.uint8),
                         mode="constant", cval=0) == 9
        mask = good9 & (np.abs(sums) > 5 * np.sqrt(variances))
        return mask, background, {"valid_pixels": int(valid.sum()), "source_pixels": int(mask.sum()),
                                  "unscored_3x3_boxes_invalid_quality_or_edge": int((~good9).sum())}


def near_source(wcs, mask, candidate):
    center = wcs.all_world2pix([[candidate["ra"], candidate["dec"]]], 1)[0]
    cx, cy = center
    rad = math.ceil(2.4 / .11) + 4
    xmin, xmax = max(1, int(cx) - rad), min(NPIX, int(cx) + rad)
    ymin, ymax = max(1, int(cy) - rad), min(NPIX, int(cy) + rad)
    iy, ix = np.nonzero(mask[ymin-1:ymax, xmin-1:xmax])
    if len(ix) == 0:
        return False
    pix = np.column_stack((ix + xmin, iy + ymin))
    corners=[]
    for dx,dy in [(-.5,-.5),(.5,-.5),(.5,.5),(-.5,.5)]:
        rd=wcs.all_pix2world(pix+[dx,dy],1)
        east,north=tangent(rd[:,0],rd[:,1],RA0,DEC0)
        corners.append(np.column_stack((east-candidate["east"],north-candidate["north"])))
    polygon=np.stack(corners,axis=1)
    # Distance from the lattice center to the actual WCS-mapped flagged
    # pixel quadrilateral, rather than merely to its center.
    end=np.roll(polygon,-1,axis=1)
    vec=end-polygon
    t=np.clip(-np.sum(polygon*vec,axis=2)/np.sum(vec*vec,axis=2),0,1)
    closest=polygon+t[:,:,None]*vec
    d=np.sqrt(np.sum(closest*closest,axis=2)).min(axis=1)
    cross=polygon[:,:,0]*vec[:,:,1]-polygon[:,:,1]*vec[:,:,0]
    inside=np.all(cross>=0,axis=1)|np.all(cross<=0,axis=1)
    d[inside]=0
    return bool(np.any(d<=2.4))


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
        bit_counts = {str(1 << bit): int(np.count_nonzero((dqbits & (1 << bit)) != 0))
                      for bit in range(16) if np.any((dqbits & (1 << bit)) != 0)}
        full_a = float(np.dot(op.aperture, op.pam))
        good_a = float(np.dot(op.aperture * op.pam, v))
        full_b = float(np.sum(op.annulus))
        good_b = float(np.dot(op.annulus, v))
        return v, {"strict": bool(np.all(v)), "aperture_coverage": good_a / full_a if full_a > 0 else -1,
                   "annulus_coverage": good_b / full_b if full_b > 0 else -1,
                   "dq_values_at_invalid": bits, "dq_bit_counts": bit_counts,
                   "nonfinite_sci": int(np.count_nonzero(~np.isfinite(sci))),
                   "nonfinite_err": int(np.count_nonzero(~np.isfinite(err))),
                   "nonpositive_err": int(np.count_nonzero(np.isfinite(err) & (err <= 0))),
                   "n_invalid": int((~v).sum())}


def synthetic_operator_gate(op, v, kappa):
    a, b, p = op.aperture, op.annulus, op.pam
    vv = v.astype(float)
    A, B = float(np.dot(vv * a, p)), float(np.dot(vv, b))
    if A <= 0 or B <= 0:
        return False
    w = kappa * (vv * a * p - (A / B) * vv * b)
    if abs(w.sum()) > 1e-10 * max(1, np.sum(abs(w))):
        return False
    x = 2 + .0001 * (op.flat % NPIX) + .0002 * (op.flat // NPIX)
    err = 1 + .0001 * (op.flat % NPIX)
    direct = kappa * (np.dot(vv * a * p, x) - A * np.dot(vv * b, x) / B)
    variance = kappa**2 * np.dot((vv * a * p - A / B * vv * b)**2, err**2)
    inject = np.zeros(len(w))
    j = int(np.argmax(abs(w)))
    inject[j] = -.375  # signed known test pulse, independent of FITS SCI
    injection_change = np.dot(w, x + inject) - np.dot(w, x)
    return (abs(np.dot(w, x) - direct) <= 1e-10 * max(1, abs(direct)) and
            abs(np.dot(w*w, err*err) - variance) <= 1e-10 * max(1, variance) and
            abs(injection_change - w[j] * inject[j]) <= 1e-10 * max(1, abs(injection_change)) and
            abs(np.dot(2*w, x) - 2*np.dot(w, x)) <= 1e-10 * max(1, abs(np.dot(w, x))) and
            abs(kappa * np.dot(vv * a * p, np.ones(len(w))) - kappa * A) < 1e-9)


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
                masks[i], medians[i], maskmeta[e["filename"]] = make_source_mask(e["path"])
        np.savez_compressed(HERE / "source_masks.npz", **{f"exposure_{i}": z for i,z in masks.items()})
        rows, op_rows = [], []
        flats, aps, anns, ps = [], [], [], []
        total = 0
        def save_arrays(path):
            # Write one concatenated field at a time. np.savez_compressed
            # would materialize all four concatenations simultaneously.
            with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as z:
                for name, parts, dtype in [("flat", flats, np.int32),
                                           ("aperture", aps, np.float32),
                                           ("annulus", anns, np.float32),
                                           ("pam", ps, np.float32)]:
                    field = np.concatenate(parts) if parts else np.array([], dtype)
                    with z.open(name + ".npy", "w", force_zip64=True) as dst:
                        np.lib.format.write_array(dst, field, allow_pickle=False)
                    del field

        def stop_stage_a(reason, current=None):
            if current is not None:
                current["reason"] = reason
                rows.append(current)
            save_arrays(HERE / "operators.partial.npz")
            partial = {"status":"stage_a_stopped_before_heldout_aperture_sums",
                "completed_candidates":len(rows), "total_candidates":len(grid),
                "reason":reason, "rows":rows, "operator_rows":op_rows,
                "retained_operator_field_bytes":total * 16,
                "source_masks_sha256":digest(HERE / "source_masks.npz"),
                "operators_partial_sha256":digest(HERE / "operators.partial.npz"),
                "mask_metadata":maskmeta}
            (HERE / "stage_a_partial.json").write_text(json.dumps(partial, indent=2)+"\n")
            raise RuntimeError("Stage A stopped at a numerical or resource gate: " + reason)

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
                        op = make_operator(e["wcs"], c["east"], c["north"], RA0, DEC0,
                                           pam, .4, (1.2, 2.0))
                        v, q = operator_validity(e["path"], op)
                    except (ValueError, AssertionError, FloatingPointError) as exc:
                        stop_stage_a("operator_or_quality_numeric:" + repr(exc) +
                                     ":candidate=" + str(c["id"]) + ":exposure=" + e["filename"], row)
                    q["resolution"] = op.resolution
                    q["numerical_closure"] = op.closure
                    q["filename"] = e["filename"]
                    if not synthetic_operator_gate(op, np.ones(len(v), bool), e["photflam"] / exps[0]["photflam"]):
                        stop_stage_a("synthetic_full_operator_failure:candidate=" +
                                     str(c["id"]) + ":exposure=" + e["filename"], row)
                    if (q["aperture_coverage"] > 0 and q["annulus_coverage"] > 0 and
                            not synthetic_operator_gate(op, v, e["photflam"] / exps[0]["photflam"])):
                        stop_stage_a("synthetic_masked_operator_failure:candidate=" +
                                     str(c["id"]) + ":exposure=" + e["filename"], row)
                    strict &= q["strict"]
                    secondary &= q["aperture_coverage"] >= .90 and q["annulus_coverage"] >= .75
                    qualities.append(q)
                    temp.append((i, op))
                if reason is None:
                    row["strict"] = bool(strict)
                    row["secondary"] = bool(secondary)
                    row["qualities"] = qualities
                    # Design-only background stratum; never read test SCI here.
                    bg = []
                    for i in masks:
                        center = exps[i]["wcs"].all_world2pix([[c["ra"], c["dec"]]], 1)[0]
                        x, y = [int(round(z)) for z in center]
                        bg.append(float(medians[i][y-1, x-1]))
                    row["design_background_mean"] = float(np.mean(bg))
                    if strict or secondary:
                        proposed = total + sum(len(op.flat) for _, op in temp)
                        if proposed * 16 > OPERATOR_FIELDS_CAP_BYTES:
                            stop_stage_a("100MB_working_array_cap:next_candidate=" +
                                         str(c["id"]) + ":proposed_operator_bytes=" +
                                         str(proposed * 16), row)
                        for i, op in temp:
                            n = len(op.flat)
                            op_rows.append({"candidate": c["id"], "exposure": i,
                                            "start": total, "end": total+n})
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
    out = {"status": "stage_a_complete_no_heldout_aperture_sums", "elapsed_seconds": time.monotonic()-start,
           "input_freeze_sha256":digest(HERE / "input-freeze.json"),
           "bounds_east_north_arcsec": bounds, "curvature_max_native_pixels": max(curvature),
           "curvature_by_exposure": curvature, "mask_metadata": maskmeta,
           "n_geometry_grid": len(grid), "n_strict": len(strict), "n_secondary": len(secondary),
           "strict_tiles": sorted(set(r["tile"] for r in strict)),
           "secondary_tiles": sorted(set(r["tile"] for r in secondary)),
           "rows": rows, "operator_rows": op_rows,
           "exposures": [{k:v for k,v in e.items() if k != "wcs"} for e in exps],
           "source_masks_sha256": digest(HERE / "source_masks.npz"),
           "operators_sha256": digest(HERE / "operators.npz")}
    (HERE / "stage_a.json").write_text(json.dumps(out, indent=2)+"\n")
    print(json.dumps({k: out[k] for k in ["status","elapsed_seconds","n_geometry_grid","n_strict","n_secondary",
                                            "strict_tiles","secondary_tiles"]}, indent=2))


def _make_weight(op, valid, kappa):
    v = valid.astype(float)
    A = float(np.dot(v * op["aperture"], op["pam"]))
    B = float(np.dot(v, op["annulus"]))
    if A <= 0 or B <= 0:
        raise ValueError("nonpositive retained aperture/annulus area")
    w = kappa * (v * op["aperture"] * op["pam"] - (A / B) * v * op["annulus"])
    if abs(w.sum()) > 1e-9 * max(1, np.sum(abs(w))):
        raise ValueError("constant-SCI cancellation failed")
    return w


def _score_pair(pair, seed, bootstrap_n):
    d = np.array([x["d"] for x in pair], float)
    V = np.array([x["V"] for x in pair], float)
    tile = np.array([x["tile"] for x in pair], int)
    east = np.array([x["east"] for x in pair], float)
    north = np.array([x["north"] for x in pair], float)
    bg = np.array([x["design_background_mean"] for x in pair], float)
    assert np.all(np.isfinite(V)) and np.all(V > 0)

    def metrics(q):
        dd, vv = d[q], V[q]
        zz = dd / np.sqrt(vv)
        intercept = np.sum(dd / vv) / np.sum(1 / vv)
        return {"n": len(dd), "weighted_intercept": float(intercept),
                "signed_standardized_mean": float(np.mean(zz)),
                "uncentered_mean_z2": float(np.mean(zz**2)),
                "centered_variance_ratio": float(np.sum((dd - intercept)**2 / vv) / (len(dd)-1)) if len(dd)>1 else None,
                "sum_d": float(np.sum(dd)), "sum_d2": float(np.sum(dd**2)), "sum_V": float(np.sum(vv)),
                "sum_d2_minus_sum_V": float(np.sum(dd**2)-np.sum(vv)),
                "sum_d2_over_sum_V": float(np.sum(dd**2)/np.sum(vv))}

    overall = metrics(np.arange(len(d)))
    uniques = sorted(set(tile))
    tile_counts = {str(t): int(np.sum(tile==t)) for t in uniques}
    resampling_allowed = len(uniques) >= 12 and all(n>=3 for n in tile_counts.values())
    deleted=[]
    boot=[]
    if resampling_allowed:
        deleted = [{"removed_tile":t, **metrics(np.flatnonzero(tile != t))}
                   for t in uniques if np.sum(tile != t) > 1]
        rng = np.random.default_rng(seed)
        for rep in range(bootstrap_n):
            draw = rng.choice(uniques, size=len(uniques), replace=True)
            idx = np.concatenate([np.flatnonzero(tile == t) for t in draw])
            if len(idx) > 1:
                boot.append({"replicate":rep, **metrics(idx)})
    overall["tile_delete_one_centered_ratio_range"] = [float(np.min([z["centered_variance_ratio"] for z in deleted])),
                                                          float(np.max([z["centered_variance_ratio"] for z in deleted]))] if deleted else None
    overall["tile_bootstrap_centered_ratio_2p5_97p5"] = np.quantile(
        [z["centered_variance_ratio"] for z in boot], [.025,.975]).tolist() if boot else None
    overall["tile_bootstrap_valid_replicates"] = len(boot)
    overall["tile_counts"] = tile_counts
    overall["tile_resampling_allowed"] = resampling_allowed
    overall["tile_delete_one"] = deleted
    overall["tile_bootstrap_replicates"] = boot
    # Design-only background terciles, deterministic quantile boundaries.
    cuts = np.quantile(bg, [1/3, 2/3])
    strata = []
    for k in range(3):
        q = np.flatnonzero((bg <= cuts[0] if k == 0 else
                            (bg > cuts[0]) & (bg <= cuts[1]) if k == 1 else bg > cuts[1]))
        strata.append(metrics(q) if len(q) else {"n": 0})
    overall["design_background_tercile_edges"] = cuts.tolist()
    overall["design_background_strata"] = strata
    dist = np.hypot(east[:, None]-east[None, :], north[:, None]-north[None, :])
    intercept = overall["weighted_intercept"]
    products = (d[:, None]-intercept) * (d[None, :]-intercept)
    upper = np.triu(np.ones(dist.shape, bool), 1)
    sep = []
    for lo, hi in [(6,12),(12,24),(24,float("inf"))]:
        idx = upper & (dist >= lo) & (dist < hi)
        sep.append({"lo_arcsec":lo,"hi_arcsec":hi,"n_pairs":int(idx.sum()),
                    "mean_centered_pair_product":float(products[idx].mean()) if idx.any() else None})
    overall["separation_pair_products"] = sep
    overall["spatial_tile_strata"] = [{"tile":int(t), **metrics(np.flatnonzero(tile==t))}
                                       for t in uniques]
    overall["z_min_max"] = [float(np.min(d/np.sqrt(V))), float(np.max(d/np.sqrt(V)))]
    return overall


def _write_csv(path, rows):
    if not rows:
        path.write_text("")
        return
    cols = list(rows[0])
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=cols)
        writer.writeheader()
        writer.writerows(rows)


def stage_b(release_file):
    start = time.monotonic()
    frozen = verify_freeze()
    release = json.loads(release_file.read_text())
    cap = release.get("stage_b_seconds_cap")
    if not isinstance(cap,(int,float)) or cap<=0:
        raise PermissionError("root Stage B release must set a positive seconds cap")
    prior = HERE / "stage_a.json"
    if not release.get("authorized_stage_b") or release.get("stage_a_sha256") != digest(prior):
        raise PermissionError("root Stage B release missing or does not match frozen Stage A")
    state = json.loads(prior.read_text())
    assert state["status"] == "stage_a_complete_no_heldout_aperture_sums"
    assert digest(HERE / "operators.npz") == state["operators_sha256"]
    assert digest(HERE / "source_masks.npz") == state["source_masks_sha256"]
    p, amend = read_design()
    exps = exposure_metadata()
    assert [e["sha256"] for e in exps] == [e["sha256"] for e in state["exposures"]]
    with np.load(HERE / "operators.npz", allow_pickle=False) as archive:
        arr = {k: archive[k] for k in ["flat","aperture","annulus","pam"]}
    bykey = {(r["candidate"],r["exposure"]):r for r in state["operator_rows"]}
    cohorts = {"strict":[r for r in state["rows"] if r.get("strict")],
               "secondary":[r for r in state["rows"] if r.get("secondary")]}
    score_ledger = []
    for branch, cohort in cohorts.items():
        if len(cohort) < 30:
            continue
        for candidate in cohort:
            if time.monotonic()-start > cap:
                _write_csv(HERE / "signed_exposure_partial.csv", score_ledger)
                (HERE/"stage_b_partial.json").write_text(json.dumps({"phase":"primary_exposure_ledger",
                    "completed_exposure_rows":len(score_ledger),"reason":"root_stage_b_seconds_cap",
                    "seconds_cap":cap,
                    "signed_exposure_partial_sha256":digest(HERE / "signed_exposure_partial.csv")},indent=2)+"\n")
                raise TimeoutError("Stage B root resource cap before primary completion")
            for i, e in enumerate(exps):
                z = bykey[(candidate["id"],i)]
                lo,hi = z["start"],z["end"]
                op = {k:v[lo:hi] for k,v in arr.items()}
                with fits.open(e["path"], memmap=True, do_not_scale_image_data=True) as f:
                    flat = op["flat"]
                    sci = f["SCI"].data.ravel()[flat]
                    err = f["ERR"].data.ravel()[flat]
                    dq = f["DQ"].data.ravel()[flat]
                    v = np.isfinite(sci) & np.isfinite(err) & (err>0) & (dq==0)
                    if branch == "strict" and not np.all(v):
                        raise ValueError("strict support changed since Stage A")
                    q = candidate["qualities"][i]
                    A = np.dot(v * op["aperture"],op["pam"]) / np.dot(op["aperture"],op["pam"])
                    B = np.dot(v,op["annulus"]) / np.sum(op["annulus"])
                    if abs(A-q["aperture_coverage"])>1e-6 or abs(B-q["annulus_coverage"])>1e-6:
                        raise ValueError("support coverage changed since Stage A")
                    kappa = e["photflam"] / exps[0]["photflam"]
                    w = _make_weight(op, v if branch == "secondary" else np.ones_like(v), kappa)
                    fsum = float(np.dot(w[v],sci[v]))
                    variance = float(np.dot(w[v]**2,err[v]**2))
                    if not np.isfinite(variance) or variance <= 0 or not np.isfinite(fsum):
                        raise ValueError("invalid measured flux or variance")
                    score_ledger.append({"branch":branch,"candidate":candidate["id"],"visit":e["visit"],
                        "ordinal":e["ordinal"],"filename":e["filename"],"east":candidate["east"],
                        "north":candidate["north"],"tile":candidate["tile"],
                        "design_background_mean":candidate["design_background_mean"],
                        "aperture_coverage":float(A),"annulus_coverage":float(B),"flux":fsum,"V":variance})
    _write_csv(HERE / "signed_exposure_ledger.csv", score_ledger)
    keyed = {(r["branch"],r["candidate"],r["visit"],r["ordinal"]):r for r in score_ledger}
    output = {"status":"stage_b_scored_after_root_release","release_file":str(release_file),
              "release_sha256":digest(release_file),"branches":{}}
    pair_rows = []
    template_rows = []
    Hcommon = np.array([[1,0,-.5,-.5],[0,1,-.5,-.5]])
    Hseparate = np.array([[1,0,-1,0],[0,1,0,-1]])
    for branch, cohort in cohorts.items():
        tilecounts={str(t):sum(r["tile"]==t for r in cohort) for t in sorted(set(r["tile"] for r in cohort))}
        entry = {"support_n":len(cohort),"support_tiles":sorted(set(r["tile"] for r in cohort)),
                 "tile_counts":tilecounts,
                 "primary_adequacy":"full" if len(cohort)>=100 and len(tilecounts)>=12 and all(n>=3 for n in tilecounts.values())
                 else "limited" if len(cohort)>=30 else "stop_below_30"}
        if len(cohort)<30:
            output["branches"][branch]=entry
            continue
        byvisit = {}
        for visit in ["search","template"]:
            pairs=[]
            for c in cohort:
                p2,p4=[keyed[(branch,c["id"],visit,j)] for j in [2,4]]
                row={"branch":branch,"visit":visit,"candidate":c["id"],"east":c["east"],
                     "north":c["north"],"tile":c["tile"],
                     "design_background_mean":c["design_background_mean"],
                     "d":p2["flux"]-p4["flux"],"V":p2["V"]+p4["V"]}
                row["z"]=row["d"]/math.sqrt(row["V"])
                pairs.append(row)
            byvisit[visit]=_score_pair(pairs,p["bootstrap"]["seed"],p["bootstrap"]["spatial_tile_resamples"])
            pair_rows+=pairs
        entry["visits"]=byvisit
        entry["count_weighted_pooled_uncentered_mean_z2"] = float(np.average(
            [byvisit[v]["uncentered_mean_z2"] for v in ["search","template"]],
            weights=[byvisit[v]["n"] for v in ["search","template"]]))
        # Four-exposure common-template algebra, retaining H D H^T.
        for c in cohort:
            s2,s4,t2,t4=[keyed[(branch,c["id"],visit,j)] for visit,j in
                            [("search",2),("search",4),("template",2),("template",4)]]
            x=np.array([q["flux"] for q in [s2,s4,t2,t4]])
            var=np.array([q["V"] for q in [s2,s4,t2,t4]])
            common=Hcommon@x; separate=Hseparate@x
            cc=Hcommon@np.diag(var)@Hcommon.T
            sc=Hseparate@np.diag(var)@Hseparate.T
            delta=float(common[0]*common[1]-separate[0]*separate[1])
            decomposition=float(.5*(x[0]-x[1])*(x[3]-x[2])+.25*(x[2]-x[3])**2)
            if abs(delta-decomposition)>1e-8*max(1,abs(delta)):
                raise AssertionError("common-template algebra failed")
            template_rows.append({"branch":branch,"candidate":c["id"],"east":c["east"],"north":c["north"],
                "S2":x[0],"S4":x[1],"T2":x[2],"T4":x[3],
                "common1":common[0],"common2":common[1],"separate1":separate[0],"separate2":separate[1],
                "common_cov_offdiag":cc[0,1],"common_cov_diag1":cc[0,0],"common_cov_diag2":cc[1,1],
                "separate_cov_offdiag":sc[0,1],"crossproduct_difference":delta,
                "predicted_common_cov_offdiag":.25*(var[2]+var[3]),"decomposition":decomposition,
                "same_template_difference_minus_search_pair":float((common[0]-common[1])-(x[0]-x[1]))})
        entry["common_template"]={"n":len(cohort),
            "mean_crossproduct_difference":float(np.mean([r["crossproduct_difference"] for r in template_rows if r["branch"]==branch])),
            "mean_predicted_offdiag":float(np.mean([r["predicted_common_cov_offdiag"] for r in template_rows if r["branch"]==branch]))}
        output["branches"][branch]=entry
    _write_csv(HERE / "signed_pair_ledger.csv",pair_rows)
    _write_csv(HERE / "common_template_ledger.csv",template_rows)
    # Durable primary checkpoint before dependent sensitivities: a later
    # numerical or runtime gate cannot erase the already frozen main result.
    primary_checkpoint={**output,"elapsed_seconds_to_primary":time.monotonic()-start,
        "signed_exposure_ledger_sha256":digest(HERE / "signed_exposure_ledger.csv"),
        "signed_pair_ledger_sha256":digest(HERE / "signed_pair_ledger.csv"),
        "common_template_ledger_sha256":digest(HERE / "common_template_ledger.csv")}
    (HERE / "stage_b_primary.json").write_text(json.dumps(primary_checkpoint,indent=2)+"\n")
    # Dependent strict-primary sensitivities; frozen locations and original
    # source masks, with quality support reported afresh at each radius.
    strict = cohorts["strict"]
    sensitivity_rows = []
    sensitivity = {}
    if len(strict) >= 30:
        reverse = {}
        for visit in ["search","template"]:
            pairs=[]
            for c in strict:
                p1,p3=[keyed[("strict",c["id"],visit,j)] for j in [1,3]]
                row={"branch":"reverse_same_original_mask","visit":visit,"candidate":c["id"],
                     "east":c["east"],"north":c["north"],"tile":c["tile"],
                     "design_background_mean":c["design_background_mean"],
                     "d":p1["flux"]-p3["flux"],"V":p1["V"]+p3["V"]}
                row["z"]=row["d"]/math.sqrt(row["V"])
                pairs.append(row)
            reverse[visit]=_score_pair(pairs,p["bootstrap"]["seed"],p["bootstrap"]["spatial_tile_resamples"])
            sensitivity_rows += pairs
        sensitivity["split_reversal_on_original_mask"] = {
            "interpretation":"dependent: scored design exposures on original design-selected mask; no new independent source mask",
            "visits":reverse}
        with fits.open(PAM, memmap=True) as pfile:
            pam=pfile[1].data
            for radius in p["sensitivities"]["aperture_radii_arcsec"]:
                rows_by_visit={"search":[],"template":[]}
                failures=[]
                for c in strict:
                    per={}
                    valid_all=True
                    for i,e in enumerate(exps):
                        op=make_operator(e["wcs"],c["east"],c["north"],RA0,DEC0,pam,radius,(1.2,2.0))
                        v,q=operator_validity(e["path"],op)
                        if not q["strict"]:
                            valid_all=False
                            failures.append({"candidate":c["id"],"filename":e["filename"],
                                             "reason":"radius_specific_invalid_pixel","invalid":q["n_invalid"]})
                        if not synthetic_operator_gate(op,np.ones(len(v),bool),e["photflam"]/exps[0]["photflam"]):
                            raise ValueError("radius synthetic operator closure failure")
                        if e["ordinal"] not in (2,4):
                            continue
                        with fits.open(e["path"],memmap=True,do_not_scale_image_data=True) as h:
                            sci=h["SCI"].data.ravel()[op.flat]
                            err=h["ERR"].data.ravel()[op.flat]
                            # A candidate is scored only if all eight exposures
                            # pass; retain these values locally until then.
                            if q["strict"]:
                                w=_make_weight({"aperture":op.aperture,"annulus":op.annulus,
                                                "pam":op.pam},np.ones(len(v),bool),
                                               e["photflam"]/exps[0]["photflam"])
                                per[(e["visit"],e["ordinal"])]=(float(np.dot(w,sci)),
                                                                    float(np.dot(w*w,err*err)))
                    if not valid_all:
                        continue
                    for visit in ["search","template"]:
                        f2,v2=per[(visit,2)]
                        f4,v4=per[(visit,4)]
                        row={"branch":f"radius_{radius}","visit":visit,"candidate":c["id"],
                             "east":c["east"],"north":c["north"],"tile":c["tile"],
                             "design_background_mean":c["design_background_mean"],
                             "d":f2-f4,"V":v2+v4}
                        row["z"]=row["d"]/math.sqrt(row["V"])
                        rows_by_visit[visit].append(row)
                        sensitivity_rows.append(row)
                    if time.monotonic()-start > cap:
                        (HERE/"stage_b_partial.json").write_text(json.dumps({"phase":"declared_sensitivities",
                            "primary_checkpoint_sha256":digest(HERE/"stage_b_primary.json"),
                            "reason":"root_stage_b_seconds_cap","seconds_cap":cap},indent=2)+"\n")
                        raise TimeoutError("Stage B root resource cap; primary checkpoint preserved")
                key=f"radius_{radius}"
                sensitivity[key]={"base_locations":len(strict),"retained_locations":len(rows_by_visit["search"]),
                                  "failures":failures,"visits":{v:_score_pair(rows_by_visit[v],p["bootstrap"]["seed"],
                                       p["bootstrap"]["spatial_tile_resamples"]) for v in rows_by_visit
                                       if len(rows_by_visit[v])>=30}}
    _write_csv(HERE / "sensitivity_pairs.csv",sensitivity_rows)
    output["declared_sensitivities"] = sensitivity
    output["elapsed_seconds"]=time.monotonic()-start
    (HERE / "stage_b.json").write_text(json.dumps(output,indent=2)+"\n")
    print(json.dumps({k:v for k,v in output.items() if k!="branches"},indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=["a", "b"])
    parser.add_argument("--release-file", type=Path)
    args = parser.parse_args()
    if args.stage == "a":
        stage_a()
    else:
        if args.release_file is None:
            raise SystemExit("Stage B requires root release file")
        stage_b(args.release_file)


if __name__ == "__main__":
    main()
