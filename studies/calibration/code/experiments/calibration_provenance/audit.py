"""Read-only, pickle-free audit of local Dovekie calibration artifacts."""
from __future__ import annotations

import gzip
import hashlib
import io
import json
import tarfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from astropy.io import fits
from numpy.lib import format as npformat
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
DOV = ROOT / "sources/repos/bap37__Dovekie"
DES = ROOT / "sources/repos/des-science__DES-SN5YR@1.3"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(lambda: f.read(1024 * 1024), b""):
            h.update(part)
    return h.hexdigest()


def array_headers(path: Path) -> dict:
    result = {}
    with zipfile.ZipFile(path) as archive:
        for member in archive.infolist():
            assert member.filename.endswith(".npy") and "/" not in member.filename
            stream = io.BytesIO(archive.read(member.filename))
            version = npformat.read_magic(stream)
            if version == (1, 0):
                shape, fortran, dtype = npformat.read_array_header_1_0(stream)
            elif version in ((2, 0), (3, 0)):
                shape, fortran, dtype = npformat.read_array_header_2_0(stream)
            else:
                raise ValueError(version)
            assert not dtype.hasobject, (path, member.filename, dtype)
            result[member.filename] = {
                "shape": list(shape), "dtype": str(dtype), "fortran_order": fortran,
                "uncompressed_bytes": member.file_size,
            }
    return result


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    names = ("DOVEKIE.V9.npz", "DOVEKIE.V9.2.WD.npz", "DOVEKIE_COV_V9.0.npz")
    files = {n: {"sha256": sha(DOV / n), "members": array_headers(DOV / n)}
             for n in names}
    with np.load(DOV / names[0], allow_pickle=False) as z:
        v9 = z["samples"].copy()
        labels = z["labels"].astype(str)
        surveys = z["surveys_for_chisq"].astype(str)
    with np.load(DOV / names[1], allow_pickle=False) as z:
        wd = z["samples"].copy()
        wdlabels = z["labels"].astype(str)
        wdsurveys = z["surveys_for_chisq"].astype(str)
    with np.load(DOV / names[2], allow_pickle=False) as z:
        covariance = z["cov"].copy()
        covlabels = z["labels"].astype(str)
    assert v9.shape == wd.shape == (7000, 52)
    assert covariance.shape == (52, 52)
    assert np.array_equal(labels, wdlabels) and np.array_equal(surveys, wdsurveys)
    assert np.array_equal(covlabels, np.char.replace(np.char.replace(
        labels, "_", " "), "offset", "O"))
    assert np.isfinite(v9).all() and np.isfinite(wd).all() and np.isfinite(covariance).all()

    mean_v9 = v9.mean(axis=0)
    mean_wd = wd.mean(axis=0)
    cov_v9 = np.cov(v9, rowvar=False)
    cov_wd = np.cov(wd, rowvar=False)
    rows = pd.DataFrame({
        "index": np.arange(len(labels)), "label": labels, "covariance_label": covlabels,
        "v9_mean_mag": mean_v9, "v9_sd_mag": v9.std(axis=0, ddof=1),
        "wd_mean_mag": mean_wd, "wd_sd_mag": wd.std(axis=0, ddof=1),
        "wd_minus_v9_mean_mag": mean_wd - mean_v9,
        "stored_cov_diag_mag2": np.diag(covariance),
    })
    rows.to_csv(OUT / "offsets.csv", index=False)
    largest = rows.iloc[np.argmax(np.abs(rows.wd_minus_v9_mean_mag))]
    configs = yaml.safe_load((DOV / "DOVEKIE_DEFS.yml").read_text())
    des_indices = np.array([int(np.where(labels == f"DES-{band}_offset")[0][0])
                            for band in "griz"])
    projector = np.eye(4) - np.ones((4, 4)) / 4
    contrasts = np.array([[1., -1., 0., 0.], [0., -1., 1., 0.],
                          [0., -1., 0., 1.]])
    des_means_v9, des_means_wd = mean_v9[des_indices], mean_wd[des_indices]
    des_cov_v9 = covariance[np.ix_(des_indices, des_indices)]
    des_cov_wd = cov_wd[np.ix_(des_indices, des_indices)]
    vectors_path = OUT / "des_fourband_vectors_cov.npz"
    np.savez(vectors_path,
             bands=np.array(list("griz")),
             v9_mean_mag=des_means_v9, wd_mean_mag=des_means_wd,
             wd_minus_v9_mean_mag=des_means_wd - des_means_v9,
             v9_cov_mag2=des_cov_v9, wd_cov_mag2=des_cov_wd,
             gray_projector=projector,
             v9_grayfree_mean_mag=projector @ des_means_v9,
             wd_grayfree_mean_mag=projector @ des_means_wd,
             v9_grayfree_cov_mag2=projector @ des_cov_v9 @ projector,
             wd_grayfree_cov_mag2=projector @ des_cov_wd @ projector,
             gminus_r_iminus_r_zminus_r_contrast_matrix=contrasts,
             v9_relative_to_r_mag=contrasts @ des_means_v9,
             wd_relative_to_r_mag=contrasts @ des_means_wd,
             v9_relative_to_r_cov_mag2=contrasts @ des_cov_v9 @ contrasts.T,
             wd_relative_to_r_cov_mag2=contrasts @ des_cov_wd @ contrasts.T)

    raw = pd.read_csv(DOV / "output_observed_apermags/DES_observed.csv")
    av = pd.read_csv(DOV / "output_observed_apermags+AV/DES_observed.csv")
    separation, raw_index = cKDTree(raw[["RA", "DEC"]]).query(av[["RA", "DEC"]])
    magdiff = np.max(np.abs(raw.iloc[raw_index][["DES-g", "DES-r", "DES-i", "DES-z"]].to_numpy()
                            - av[["DES-g", "DES-r", "DES-i", "DES-z"]].to_numpy()))
    assert len(raw) == 2000 and len(av) == 1255
    assert len(set(raw_index)) == len(av) and np.all(np.diff(raw_index) > 0)
    assert separation.max() < 0.0001 and magdiff < 0.0001

    salt_archive = DOV / "SALT3.DOVEKIE.tar.gz"
    release_salt = DES / "2_LCFIT_MODEL/SALT3.DES5YR"
    salt_files = {}
    with tarfile.open(salt_archive, "r:gz") as archive:
        options = json.load(archive.extractfile("SALT3.DOVEKIE/options.json"))
        for name in ("SALT3.INFO", "salt3_template_0.dat", "salt3_color_correction.dat"):
            dov_bytes = archive.extractfile("SALT3.DOVEKIE/" + name).read()
            des_path = release_salt / (name if name == "SALT3.INFO" else name + ".gz")
            des_bytes = des_path.read_bytes() if name == "SALT3.INFO" else gzip.open(des_path, "rb").read()
            salt_files[name] = {
                "dovekie_sha256": hashlib.sha256(dov_bytes).hexdigest(),
                "des_release_sha256": hashlib.sha256(des_bytes).hexdigest(),
                "equal": dov_bytes == des_bytes,
            }
        testing_log = archive.extractfile("SALT3.DOVEKIE/testing.log").read().decode(
            "utf-8", errors="replace")
    log_has_no_shift = "No calibration shift file provided, continuing" in testing_log
    des_kcor_path = ROOT / "phase2/official/inputs/SNDATA_ROOT/kcor/DES/DES-SN5YR/calib_DES-SN5YR_DES.fits.gz"
    with fits.open(des_kcor_path, memmap=False) as hdus:
        zp = hdus["ZPoff"].data
        des_kcor_zpoff = {
            str(row["Filter Name"]).strip(): float(row["ZPoff(Primary)"])
            for row in zp
        }

    result = {
        "scope": "Safe local artifact/provenance audit; no model or dust-map execution",
        "audit_script_sha256": sha(Path(__file__)),
        "source_snapshot_commit": json.loads(
            (ROOT / "catalog/repositories/bap37__Dovekie.commit.json").read_text())["sha"],
        "pinned_paper_sha256": sha(ROOT / "papers/pdf/2506.05471v1.pdf"),
        "npz_files": files,
        "posterior": {
            "surveys": surveys.tolist(), "n_draws": len(v9), "n_offsets": len(labels),
            "v9_and_wd_labels_equal": True, "v9_and_wd_draws_identical": bool(np.array_equal(v9, wd)),
            "largest_absolute_mean_shift": {
                "label": largest.label, "wd_minus_v9_mag": float(largest.wd_minus_v9_mean_mag)},
            "des_offset_rows": rows.loc[rows.label.str.startswith("DES-")].to_dict("records"),
            "des_fourband_vectors_cov_file": str(vectors_path.relative_to(ROOT)),
            "des_fourband_vectors_cov_sha256": sha(vectors_path),
            "des_grayfree_v9_mmag": (1000 * projector @ des_means_v9).tolist(),
            "des_grayfree_wd_mmag": (1000 * projector @ des_means_wd).tolist(),
            "des_relative_to_r_v9_mmag": (1000 * contrasts @ des_means_v9).tolist(),
            "des_relative_to_r_wd_mmag": (1000 * contrasts @ des_means_wd).tolist(),
            "covariance_labels_are_transformed_v9_labels": True,
            "stored_cov_minus_sample_cov_v9_max_abs": float(np.max(np.abs(covariance - cov_v9))),
            "stored_cov_minus_sample_cov_wd_max_abs": float(np.max(np.abs(covariance - cov_wd))),
            "covariance_min_eigenvalue": float(np.linalg.eigvalsh(covariance)[0]),
        },
        "current_config": {
            "sha256": sha(DOV / "DOVEKIE_DEFS.yml"),
            "chainsfile": configs["chainsfile"],
            "surveys": configs["surveys_for_dovekie"],
            "matches_archived_surveys": configs["surveys_for_dovekie"] == surveys.tolist(),
            "whitedwarf_obs_loc": configs["whitedwarf_obs_loc"],
            "dustlaw": configs["dustlaw"],
        },
        "star_file_link": {
            "raw_sha256": sha(DOV / "output_observed_apermags/DES_observed.csv"),
            "av_sha256": sha(DOV / "output_observed_apermags+AV/DES_observed.csv"),
            "n_raw_rows": len(raw), "n_av_rows": len(av),
            "av_matches_unique_ordered_raw_subset": True,
            "max_coordinate_difference_deg": float(separation.max()),
            "max_des_magnitude_difference_mag": float(magdiff),
            "raw_index_min": int(raw_index.min()), "raw_index_max": int(raw_index.max()),
            "raw_ra_range_deg": [float(raw.RA.min()), float(raw.RA.max())],
            "raw_dec_range_deg": [float(raw.DEC.min()), float(raw.DEC.max())],
            "av_ra_range_deg": [float(av.RA.min()), float(av.RA.max())],
            "av_dec_range_deg": [float(av.DEC.min()), float(av.DEC.max())],
            "query_input_catalog_available": (DOV / "newcatalog/PS1_in_DES5YR.fits").exists(),
            "full_fgcm_catalog_available": (
                DOV / "rawstars/Y6A1_FGCM_V3_3_1_PSF_ALL_STARS.fits").exists(),
        },
        "salt_release_link": {
            "dovekie_archive_sha256": sha(salt_archive),
            "archive_options_calibrationshiftfile": options.get("calibrationshiftfile"),
            "archive_options_calibrationcovariance": options.get("calibrationcovariance"),
            "archive_log_says_no_calibration_shift_file": log_has_no_shift,
            "archive_vs_des_release_files": salt_files,
            "des_kcor_input_sha256": sha(
                ROOT / "phase2/official/inputs/SNDATA_ROOT/kcor/DES/DES-SN5YR/calib_DES-SN5YR_DES.input"),
            "des_kcor_fits_sha256": sha(des_kcor_path),
            "des_kcor_fits_zpoff_primary_mag": des_kcor_zpoff,
        },
        "source_sha256": {
            "dovekie_py": sha(DOV / "dovekie.py"),
            "prep_des_mcmc_py": sha(DOV / "PREP_DES_MCMC.py"),
            "query_ps1_des_py": sha(DOV / "queryps1inDES5yr.py"),
            "plotresults_py": sha(DOV / "plotresults.py"),
            "surfaces_py": sha(DOV / "surfaces-dovekie.py"),
        },
    }
    (OUT / "summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "cov_v9_max_abs": result["posterior"]["stored_cov_minus_sample_cov_v9_max_abs"],
        "n_raw": len(raw), "n_av": len(av), "star_link": True,
        "current_config_matches_archived_surveys": result["current_config"]["matches_archived_surveys"],
        "largest_mean_shift": result["posterior"]["largest_absolute_mean_shift"],
    }, indent=2))


if __name__ == "__main__":
    main()
