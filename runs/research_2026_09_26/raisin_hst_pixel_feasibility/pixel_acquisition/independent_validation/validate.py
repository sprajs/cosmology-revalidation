"""Metadata-only verification of the ten frozen MAST FITS products.

Never access HDU.data; integrity checks stream file bytes and run fitsverify.
"""
from __future__ import annotations

import hashlib
import json
import pathlib
import re
import subprocess

from astropy.io import fits


BASE = pathlib.Path(__file__).resolve().parent.parent
OUT = pathlib.Path(__file__).resolve().parent
parent_protocol = json.loads((BASE / "protocol.json").read_text())
parent_rows = json.loads((BASE / "result.json").read_text())
expected = parent_protocol["filenames"]
assert len(expected) == len(set(expected)) == 10
assert {x["filename"] for x in parent_rows} == set(expected)


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


rows = []
for filename in expected:
    parent = next(x for x in parent_rows if x["filename"] == filename)
    sidecar = json.loads((BASE / (filename + ".acquisition.json")).read_text())
    path = BASE / filename
    digest = sha256(path)
    v = subprocess.run(["fitsverify", "-e", str(path)], text=True, capture_output=True, check=False)
    with fits.open(path, memmap=True, lazy_load_hdus=True,
                   do_not_scale_image_data=True) as hdus:
        primary = hdus[0].header
        hdu_layout = []
        for hdu in hdus:
            hdr = hdu.header
            hdu_layout.append({
                "name": hdu.name,
                "ver": hdr.get("EXTVER"),
                "xtension": hdr.get("XTENSION", "PRIMARY"),
                "bitpix": hdr.get("BITPIX"),
                "naxis": hdr.get("NAXIS"),
                "shape_header": [hdr.get(f"NAXIS{i}") for i in range(1, hdr.get("NAXIS", 0) + 1)],
                "bunit": hdr.get("BUNIT"),
                "checksum_present": "CHECKSUM" in hdr,
                "datasum_present": "DATASUM" in hdr,
            })
        sci = hdus["SCI"].header
        pkeys = ["FILENAME", "ROOTNAME", "ASN_ID", "TARGNAME", "PROPOSID", "INSTRUME", "DETECTOR",
                 "FILTER", "EXPSTART", "EXPEND", "EXPTIME", "CAL_VER", "DRIZCORR", "NLINCORR",
                 "PHOTFLAM", "PHOTPLAM", "IDCTAB", "PFLTFILE", "DFLTFILE", "DGEOFILE", "NPOLFILE",
                 "D2IMFILE", "PAMFILE", "NDRIZIM", "MDRIZSKY", "DRIZSCAL", "DRIZPIXF"]
        skeys = ["WCSAXES", "WCSNAME", "CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2",
                 "CD1_1", "CD1_2", "CD2_1", "CD2_2", "A_ORDER", "B_ORDER", "AP_ORDER", "BP_ORDER",
                 "CPDIS1", "CPDIS2", "D2IMDIS1", "D2IMDIS2", "IDCTAB"]
        driz_inputs = [primary.get(f"D{i:03d}DATA") for i in range(1, 5)] if "_drz" in filename else []
        driz_other = [{k: primary.get(f"D{i:03d}{k}") for k in ["VER", "GEOM", "PIXF", "KERN", "COEF", "WKEY"]}
                      for i in range(1, 5)] if driz_inputs else []
        row = {
            "filename": filename,
            "uri_match": parent["uri"] == "mast:HST/product/" + filename,
            "size": path.stat().st_size,
            "size_parent_match": path.stat().st_size == parent["bytes"],
            "size_http_match": path.stat().st_size == int(parent["headers"]["Content-Length"]),
            "size_sidecar_match": path.stat().st_size == sidecar["bytes"],
            "sha256": digest,
            "sha256_parent_match": digest == parent["sha256"],
            "sha256_sidecar_match": digest == sidecar["sha256"],
            "fitsverify_errors_only_exit": v.returncode,
            "fitsverify_errors_only_tail": v.stdout[-700:],
            "hdu_layout": hdu_layout,
            "primary": {k: primary.get(k) for k in pkeys},
            "science_wcs": {k: sci.get(k) for k in skeys},
            "drizzle_inputs": driz_inputs,
            "drizzle_parameters": driz_other,
            "embedded_wcs_extensions": [h.name for h in hdus if h.name in ("HDRLET", "WCSCORR")],
        }
        rows.append(row)

roots = {r["filename"] for r in rows}
for row in rows:
    if row["drizzle_inputs"]:
        parsed = [re.sub(r"\[.*", "", x or "") for x in row["drizzle_inputs"]]
        row["drizzle_input_filenames"] = parsed
        row["four_input_membership_exact"] = (len(set(parsed)) == 4 and
                                               all(p in roots for p in parsed))
        row["input_association_match"] = all(
            next(x for x in rows if x["filename"] == p)["primary"]["ASN_ID"] == row["primary"]["ASN_ID"]
            for p in parsed if p in roots)

summary = {
    "n_expected": len(expected),
    "n_verified": len(rows),
    "all_identity_checks": all(all(r[k] for k in ["uri_match", "size_parent_match", "size_http_match",
                                                   "size_sidecar_match", "sha256_parent_match", "sha256_sidecar_match"])
                               for r in rows),
    "all_fitsverify_error_only_clean": all(r["fitsverify_errors_only_exit"] == 0 for r in rows),
    "embedded_checksum_hdus": sum(h["checksum_present"] for r in rows for h in r["hdu_layout"]),
    "embedded_datasum_hdus": sum(h["datasum_present"] for r in rows for h in r["hdu_layout"]),
    "drizzle_four_input_exact": all(r["four_input_membership_exact"] and r["input_association_match"]
                                    for r in rows if r["drizzle_inputs"]),
    "total_bytes": sum(r["size"] for r in rows),
}
(OUT / "result.json").write_text(json.dumps({"summary": summary, "files": rows}, indent=2) + "\n")
print(json.dumps(summary, indent=2))
