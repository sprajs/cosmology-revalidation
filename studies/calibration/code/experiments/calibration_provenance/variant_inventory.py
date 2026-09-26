"""Read-only, pickle-free inventory of released DES calibration/SALT3 variants."""
from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import zipfile
from pathlib import Path

import numpy as np
from numpy.lib import format as npformat

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
DES = ROOT / "sources/repos/des-science__DES-SN5YR@1.3"
SYS = DES / "2_LCFIT_MODEL/SALT3.DES5YR-SYS"
PANT = ROOT / "sources/repos/PantheonPlusSH0ES__DataRelease@7fc6805/Pantheon+_Data/2_CALIBRATION/FRAGILISTIC_COVARIANCE.npz"
PANT_UNSUFFIXED = ROOT / "sources/repos/PantheonPlusSH0ES__DataRelease/Pantheon+_Data/2_CALIBRATION/FRAGILISTIC_COVARIANCE.npz"
KCOR = ROOT / "phase2/official/inputs/SNDATA_ROOT/kcor/DES/DES-SN5YR"
FILTER = ROOT / "phase2/official/inputs/SNDATA_ROOT/filters/DES/DES-SN3YR_DECam"
FITOPTS = DES / "7_PIPPIN_FILES/base_files/lcfit/fitopts.yml"
FITNML = DES / "7_PIPPIN_FILES/base_files/lcfit/lcfit_desSMP_5yr.nml"
SOURCES = [
    DES / "7_PIPPIN_FILES/D5yr_analysis.yml",
    SYS / "SUBMIT.INFO",
    SYS / "DOCUMENTATION.README",
    SYS / "SCRIPTS_TRAIN.tar.gz",
    ROOT / "phase2/official/build/SNANA-2fe0f56/src/genmag_SALT2.c",
    ROOT / "phase2/official/build/SNANA-2fe0f56/src/genmag_SEDtools.c",
    ROOT / "phase2/official/build/SNANA-2fe0f56/src/snlc_fit.car",
    ROOT / "phase2/official/build/SNANA-2fe0f56/util/create_covariance.py",
]


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def parse_shifts(path: Path, colon: bool) -> list[tuple[str, str, str, str]]:
    rows = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) != 4:
            continue
        key = fields[0].rstrip(":") if colon else fields[0]
        if key in {"MAGSHIFT", "WAVESHIFT", "LAMSHIFT"}:
            rows.append((key, fields[1], fields[2], fields[3]))
    return rows


def safe_npz(path: Path) -> tuple[dict, np.ndarray, np.ndarray]:
    headers = {}
    with zipfile.ZipFile(path) as z:
        for member in z.infolist():
            assert member.filename.endswith(".npy") and "/" not in member.filename
            stream = io.BytesIO(z.read(member.filename))
            version = npformat.read_magic(stream)
            if version == (1, 0):
                shape, fortran, dtype = npformat.read_array_header_1_0(stream)
            else:
                shape, fortran, dtype = npformat.read_array_header_2_0(stream)
            assert not dtype.hasobject
            headers[member.filename] = {"shape": list(shape), "dtype": str(dtype),
                                        "fortran_order": bool(fortran)}
    with np.load(path, allow_pickle=False) as z:
        assert set(z.files) == {"cov", "labels"}
        cov = z["cov"].copy()
        labels = z["labels"].astype(str)
    assert cov.shape == (102, 102) and labels.shape == (102,)
    assert np.max(np.abs(cov - cov.T)) < 1e-18 and np.isfinite(cov).all()
    return headers, cov, labels


def main() -> None:
    fitopts = FITOPTS.read_text()
    fitnml = FITNML.read_text()
    submit = (SYS / "SUBMIT.INFO").read_text()
    assert "KCOR_FILE          = '$SNDATA_ROOT/kcor/DES/DES-SN5YR/calib_DES-SN5YR_DES.fits'" in fitnml
    assert "FITMODEL_NAME       = '$SNDATA_ROOT/models/SALT3/SALT3.DES5YR'" in fitnml
    kcor_input = (KCOR / "calib_DES-SN5YR_DES.input").read_text()
    assert "SURVEY: DES" in kcor_input
    for band in "griz":
        assert f"FILTER:  DES-{band}   DECam_{band}.dat" in kcor_input

    headers, cov, labels = safe_npz(PANT)
    assert sha(PANT) == sha(PANT_UNSUFFIXED)
    indices = {"DES3YR": [26, 27, 28, 29], "DES5YR": [30, 31, 32, 33]}
    for survey, ix in indices.items():
        assert [labels[j].strip() for j in ix] == [f"{survey} {b}" for b in "griz"]

    variants = []
    flat = []
    inputs = [PANT, PANT_UNSUFFIXED, FITOPTS, FITNML,
              KCOR / "calib_DES-SN5YR_DES.fits.gz",
              KCOR / "calib_DES-SN5YR_DES.input"] + SOURCES
    for band in "griz":
        inputs.append(FILTER / f"DECam_{band}.dat")
    for i in range(1, 10):
        shiftlist = SYS / "calib_shiftlists" / f"calib_saltshaker_{i-1}.txt"
        modeldir = SYS / f"SALT3.MODEL{i:03d}"
        info = modeldir / "SALT3.INFO"
        shifts = parse_shifts(shiftlist, False)
        embedded = parse_shifts(info, True)
        assert shifts == embedded, f"shift-list/INFO mismatch for {i}"
        assert len(shifts) in {203, 205}
        training = f"TRAINOPT{i:03d}"
        assert re.search(rf"\['{training}', None, 'SHIFTLIST_FILE\s+[^']*calib_saltshaker_{i-1}\.txt'\]", submit)
        assert re.search(rf"cal_{i}: 0\.3 FITMODEL_NAME\s+'\$SNDATA_ROOT/models/SALT3/SALT3.DES5YR-SYS/SALT3.MODEL{i:03d}'", fitopts)
        des_shift = {}
        for band in "griz":
            m = [float(v) for key, survey, filt, v in shifts
                 if key == "MAGSHIFT" and survey == "DES" and filt == f"DES-{band}"]
            w = [float(v) for key, survey, filt, v in shifts
                 if key in {"WAVESHIFT", "LAMSHIFT"} and survey == "DES" and filt == f"DES-{band}"]
            assert len(m) == len(w) == 1
            des_shift[band] = {"magshift_mag": m[0], "waveshift_angstrom": w[0]}
            flat.append({"pippin_label": f"cal_{i}", "training_label": training,
                         "model": f"SALT3.MODEL{i:03d}", "shiftlist": shiftlist.name,
                         "band": band, **des_shift[band], "pippin_weight": 0.3})
        files = {p.name: sha(p) for p in sorted(modeldir.iterdir()) if p.is_file()}
        assert "salt3_template_0.dat.gz" in files and "SALT3.INFO" in files
        inputs.extend([shiftlist] + [modeldir / name for name in files])
        variants.append({
            "pippin_label": f"cal_{i}", "pippin_weight": 0.3,
            "training_label": training, "model": f"SALT3.MODEL{i:03d}",
            "shiftlist": rel(shiftlist), "shiftlist_sha256": sha(shiftlist),
            "shift_entry_count": len(shifts), "embedded_info_entries_identical": True,
            "model_directory": rel(modeldir), "model_files_sha256": files,
            "des": des_shift,
        })
    result = {
        "schema": "des_calibration_salt_variant_inventory_v1",
        "scope": "Read-only source inventory; no fit, residual optimization or inferred prior",
        "source_sign_convention": {
            "MAGSHIFT": "SNANA adds the signed value (mag) to FILTER_SEDMODEL.magprimary and recomputes its model ZP",
            "WAVESHIFT": "SNANA adds the signed value (Angstrom) to each filter wavelength grid point",
            "applicability": "SALT3.INFO entries are parsed and applied if enabled; exact survey/full filter base must match loaded KCOR filter",
        },
        "nominal_model": rel(DES / "2_LCFIT_MODEL/SALT3.DES5YR"),
        "variants": variants,
        "fragilistic_covariance": {
            "path": rel(PANT), "sha256": sha(PANT),
            "duplicate_path": rel(PANT_UNSUFFIXED), "duplicate_equal": True,
            "members": headers, "dimension": 102,
            "max_symmetry_difference": float(np.max(np.abs(cov - cov.T))),
            "units": "not encoded in NPZ; paper describes calibration/zero-point uncertainty covariance (mag^2 interpretation)",
            "des_indices_zero_based": indices,
            "des_labels": {k: [str(labels[j]).strip() for j in v] for k, v in indices.items()},
            "des5yr_cov_block_numeric": cov[np.ix_(indices["DES5YR"], indices["DES5YR"])].tolist(),
            "all_labels": [str(x).strip() for x in labels],
            "wavelength_shift_covariance_in_npz": False,
            "exact_draw_seed_linkage_to_shiftlists": "not supplied in inspected release assets",
        },
        "kcor_prerequisite": {
            "released_fit_declares": "$SNDATA_ROOT/kcor/DES/DES-SN5YR/calib_DES-SN5YR_DES.fits",
            "local_compressed_fits": rel(KCOR / "calib_DES-SN5YR_DES.fits.gz"),
            "local_input": rel(KCOR / "calib_DES-SN5YR_DES.input"),
            "local_filter_files": [rel(FILTER / f"DECam_{b}.dat") for b in "griz"],
            "status": "DES nominal compressed KCOR/filter assets and matching DES filter names available locally; original fit execution logs and upstream training draw/seed map not established",
        },
    }
    inventory = OUT / "des_salt3_calibration_variants.json"
    inventory.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    csv_path = OUT / "des_salt3_calibration_variants.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(flat[0]))
        writer.writeheader()
        writer.writerows(flat)
    manifest = {
        "schema": "des_calibration_salt_variant_inventory_manifest_v1",
        "generator": rel(Path(__file__).resolve()),
        "generator_sha256": sha(Path(__file__).resolve()),
        "inputs_sha256": {rel(p): sha(p) for p in sorted(set(inputs))},
        "outputs_sha256": {p.name: sha(p) for p in (inventory, csv_path)},
        "checks": {"nine_variants": len(variants) == 9,
                   "all_shiftlists_equal_embedded_info": all(v["embedded_info_entries_identical"] for v in variants),
                   "four_des_bands_each": len(flat) == 36,
                   "pantheon_duplicate_equal": True},
    }
    (OUT / "des_salt3_calibration_variants_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
