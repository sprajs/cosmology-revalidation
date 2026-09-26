"""PAM versus native FLT header-WCS geometry; no science arrays are read."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


HERE = Path(__file__).resolve().parent
ACQ = HERE.parent
protocol = json.loads((HERE / "protocol.json").read_text())
grid = protocol["grid_fits_1_based_x_y"]
points = [(x, y) for y in grid for x in grid]
center = (507, 507)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def unit_vectors(world_deg):
    ra, dec = np.deg2rad(world_deg[:, 0]), np.deg2rad(world_deg[:, 1])
    return np.column_stack((np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)))


def angular_area(wcs, xy, h):
    p = np.asarray(xy, float)
    steps = np.array([[0, 0], [h, 0], [-h, 0], [0, h], [0, -h]], float)
    vectors = unit_vectors(wcs.all_pix2world(p[None, :] + steps, 1))
    dx = (vectors[1] - vectors[2]) / (2 * h)
    dy = (vectors[3] - vectors[4]) / (2 * h)
    return abs(float(np.dot(vectors[0], np.cross(dx, dy))))


with fits.open(HERE / "ir_wfc3_map.fits", memmap=True) as f:
    assert len(f) == 2 and f[1].header["NAXIS1"] == 1014 and f[1].header["NAXIS2"] == 1014
    assert f[0].header["FILETYPE"] == f[1].header["FILETYPE"] == "PIXEL AREA MAP"
    pam = f[1].data
    pcenter = float(pam[center[1] - 1, center[0] - 1])
    pvalues = np.asarray([float(pam[y - 1, x - 1]) for x, y in points])
    assert pcenter > 0 and np.all(np.isfinite(pvalues)) and np.all(pvalues > 0)
    pnorm = pvalues / pcenter
    pam_info = {"shape": list(pam.shape), "center_value": pcenter, "grid_min": float(pvalues.min()),
                "grid_max": float(pvalues.max()), "units": "dimensionless normalized pixel-area factor per STScI WFC3 handbook 9.1.2"}

with fits.open(HERE / "w3m18525i_idc.fits", memmap=True) as f:
    idc = {"filename_header": f[0].header.get("FILENAME"), "filetype": f[0].header.get("FILETYPE"),
           "useafter": f[0].header.get("USEAFTER"), "pedigree": f[0].header.get("PEDIGREE"),
           "table_rows": f[1].header.get("NAXIS2"), "table_columns": f[1].columns.names}

rows = []
for path in sorted(ACQ.glob("*_flt.fits")):
    with fits.open(path, memmap=True, lazy_load_hdus=True, do_not_scale_image_data=True) as f:
        primary, sci = f[0].header, f["SCI"].header
        assert primary["SUBARRAY"] is False and primary["SUBTYPE"] == "FULLIMAG"
        assert sci["NAXIS1"] == sci["NAXIS2"] == 1014
        assert [sci[k] for k in ["LTV1", "LTV2"]] == [0, 0]
        assert [sci[k] for k in ["LTM1_1", "LTM2_2"]] == [1, 1]
        assert primary["IDCTAB"] == "iref$w3m18525i_idc.fits"
        assert "TAN-SIP" in sci["CTYPE1"] and "TAN-SIP" in sci["CTYPE2"]
        wcs = WCS(sci, naxis=2)
        values = {}
        for h in protocol["finite_difference_steps_pixels"]:
            acenter = angular_area(wcs, center, h)
            areas = np.array([angular_area(wcs, p, h) for p in points])
            norm = areas / acenter
            diff = norm / pnorm - 1
            values[str(h)] = {"center_area_sr": acenter, "median_absolute_relative_discrepancy": float(np.median(abs(diff))),
                              "max_absolute_relative_discrepancy": float(np.max(abs(diff))),
                              "grid_relative_discrepancy": diff.tolist(), "grid_normalized_jacobian": norm.tolist()}
        hfull, hhalf = [str(x) for x in protocol["finite_difference_steps_pixels"]]
        full = np.array(values[hfull]["grid_normalized_jacobian"])
        half = np.array(values[hhalf]["grid_normalized_jacobian"])
        step_change = abs(half / full - 1)
        rows.append({"filename": path.name, "file_sha256": sha256(path), "wcsname": sci.get("WCSNAME"),
                     "sipname": primary.get("SIPNAME"), "detector": primary.get("DETECTOR"),
                     "ltv": [sci["LTV1"], sci["LTV2"]], "h": values,
                     "half_step_max_relative_change": float(step_change.max()),
                     "half_step_median_relative_change": float(np.median(step_change))})

result = {"official_pam": pam_info, "official_idctab": idc, "grid_x_y": grid, "n_points_per_flt": len(points),
          "files": rows,
          "all_flts_median_abs_relative_discrepancy": float(np.median([v for r in rows for v in r["h"]["0.5"]["grid_relative_discrepancy"]], axis=None)),
          "all_flts_max_abs_relative_discrepancy": float(max(r["h"]["0.5"]["max_absolute_relative_discrepancy"] for r in rows)),
          "all_flts_max_half_step_change": float(max(r["half_step_max_relative_change"] for r in rows))}
# Replace the aggregate signed median with the declared absolute-relative median.
result["all_flts_median_abs_relative_discrepancy"] = float(np.median(
    [abs(v) for r in rows for v in r["h"]["0.5"]["grid_relative_discrepancy"]]))
(HERE / "result.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({k: v for k, v in result.items() if k.startswith("all_flts")}, indent=2))
