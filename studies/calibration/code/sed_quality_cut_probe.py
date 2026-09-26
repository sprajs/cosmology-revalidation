"""Frozen fixed-cohort DES colour/shape eligibility probe; no refits."""
from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/research_2026_09_26/sed_quality_cut_probe"
FREEZE = OUT / "freeze.json"
RELEASE = ROOT / "sources/repos/des-science__DES-SN5YR@1.3/4_DISTANCES_COVMAT/DES-SN5YR_HD+MetaData.csv"
BASE = ROOT / "runs/research_2026_09_26/astra_design/validation1020"
RESOLVED = ROOT / "runs/research_2026_09_26/sed_nonlinear_validation/resolved"
BINS = [(0.025, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.2)]
FIELDS = ["C1", "C2", "C3", "E1", "E2", "S1", "S2", "X1", "X2", "X3"]
ARMS = ["native", "observer", "sed"]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def unique(rows: list[dict[str, str]], key: str) -> dict[str, dict[str, str]]:
    result = {r[key]: r for r in rows}
    assert len(result) == len(rows), f"duplicate {key}"
    return result


def read_fitres(path: Path) -> dict[str, dict[str, str]]:
    names = None
    rows = []
    for line in path.open():
        if line.startswith("VARNAMES:"):
            names = line.split()[1:]
        elif line.startswith("SN:"):
            values = line.split()[1:]
            assert names is not None and len(values) == len(names)
            rows.append(dict(zip(names, values)))
    return unique(rows, "CID")


def status(x1: float, c: float) -> tuple[bool, bool, bool]:
    shape = -3.0 < x1 < 3.0
    color = -0.3 < c < 0.3
    return shape, color, shape and color


def zbin(z: float) -> str:
    for i, (lo, hi) in enumerate(BINS):
        if lo <= z < hi or (i == len(BINS) - 1 and lo <= z <= hi):
            return f"[{lo:g},{hi:g}{']' if i == len(BINS)-1 else ')'}"
    return "out_of_range"


def transition(a: bool, b: bool) -> str:
    return ("pass" if a else "fail") + "_to_" + ("pass" if b else "fail")


def failure_reason(shape: bool, color: bool) -> str:
    return "none" if shape and color else ("both" if not shape and not color else "x1" if not shape else "c")


def summarize(rows: list[dict]) -> dict:
    summary = {"objects": len(rows), "arms": {}, "released_inputs": {}, "transitions": {}, "by_field": {}, "by_zHEL_bin": {}}
    for arm in ARMS:
        summary["arms"][arm] = {
            "shape_pass": sum(r[f"{arm}_shape_pass"] for r in rows),
            "color_pass": sum(r[f"{arm}_color_pass"] for r in rows),
            "joint_pass": sum(r[f"{arm}_joint_pass"] for r in rows),
            "failure_reasons": dict(sorted(Counter(r[f"{arm}_failure_reason"] for r in rows).items())),
            "x1_boundary_equal": sum(abs(r[f"{arm}_x1"]) == 3.0 for r in rows),
            "c_boundary_equal": sum(abs(r[f"{arm}_c"]) == 0.3 for r in rows),
        }
    summary["arms"]["nominal_observed_diagnostic"] = {
        "shape_pass": sum(r["nominal_observed_shape_pass"] for r in rows),
        "color_pass": sum(r["nominal_observed_color_pass"] for r in rows),
        "joint_pass": sum(r["nominal_observed_joint_pass"] for r in rows),
        "max_abs_x1_offset_from_native": max(abs(r["nominal_observed_x1"]-r["native_x1"]) for r in rows),
        "max_abs_c_offset_from_native": max(abs(r["nominal_observed_c"]-r["native_c"]) for r in rows),
    }
    summary["released_inputs"] = {
        "shape_pass": sum(r["released_shape_pass"] for r in rows),
        "color_pass": sum(r["released_color_pass"] for r in rows),
        "joint_pass": sum(r["released_joint_pass"] for r in rows),
        "static_proxy_pass": sum(r["released_static_proxy_pass"] for r in rows),
        "static_proxy_failure_ids": [r["CID"] for r in rows if not r["released_static_proxy_pass"]],
        "x1_boundary_equal": sum(abs(r["released_x1"]) == 3.0 for r in rows),
        "c_boundary_equal": sum(abs(r["released_c"]) == 0.3 for r in rows),
    }
    for arm in ["observer", "sed"]:
        summary["transitions"][f"native_to_{arm}"] = dict(sorted(Counter(r[f"native_to_{arm}"] for r in rows).items()))
    summary["transitions"]["observer_to_sed"] = dict(sorted(Counter(r["observer_to_sed"] for r in rows).items()))
    summary["transitions"]["observer_sed_disagreement"] = sum(r["observer_joint_pass"] != r["sed_joint_pass"] for r in rows)
    for axis, values in [("by_field", FIELDS), ("by_zHEL_bin", [zbin(lo) for lo, _ in BINS] + ["out_of_range"])]:
        for value in values:
            subset = [r for r in rows if r["field" if axis == "by_field" else "zHEL_bin"] == value]
            summary[axis][value] = {
                "n": len(subset),
                "joint_pass": {arm: sum(r[f"{arm}_joint_pass"] for r in subset) for arm in ARMS},
                "native_to_observer": dict(sorted(Counter(r["native_to_observer"] for r in subset).items())),
                "native_to_sed": dict(sorted(Counter(r["native_to_sed"] for r in subset).items())),
                "observer_sed_disagreement": sum(r["observer_joint_pass"] != r["sed_joint_pass"] for r in subset),
            }
    return summary


def main() -> None:
    assert FREEZE.exists() and not (OUT / "result.json").exists()
    frozen = json.loads(FREEZE.read_text())
    assert frozen["status"] == "before_cut_outcomes" and frozen["cohort_rows"] == 1020
    for rel, digest in frozen["input_sha256"].items():
        assert sha(ROOT / rel) == digest, rel
    assert sha(OUT / "protocol.md") == frozen["input_sha256"]["docs/research-2026-09-26/sed-quality-cut-protocol.md"]
    resolved_manifest = json.loads((RESOLVED / "manifest.json").read_text())
    paired_key = "runs/research_2026_09_26/sed_nonlinear_validation/resolved/paired-responses.csv"
    assert resolved_manifest["outputs_sha256"][paired_key] == sha(ROOT / paired_key)

    cohort = read_csv(BASE / "cohort.csv")
    assert len(cohort) == 1020 and len({r["CID"] for r in cohort}) == 1020
    release = unique(read_csv(RELEASE), "CID")
    fitres = read_fitres(BASE / "fit.FITRES.TEXT")
    paired = [r for r in read_csv(RESOLVED / "paired-responses.csv") if r["target"] == "observed_flux"]
    pairmap = {(r["CID"], r["mode"]): r for r in paired}
    assert len(paired) == len(pairmap) == 2040
    assert set(pairmap) == {(r["CID"], m) for r in cohort for m in ["observer", "sed"]}
    assert set(fitres) == {r["CID"] for r in cohort}

    rows = []
    max_fitres_coordinate_difference = 0.0
    max_shift_identity_error = 0.0
    max_baseline_mode_difference = 0.0
    compound_fitres_field_ids = []
    for item in cohort:
        cid = item["CID"]
        released = release[cid]
        fitted = fitres[cid]
        assert int(released["IDSURVEY"]) == 10 and int(item["IDSURVEY"]) == 10
        assert item["field"] in FIELDS and item["field"] in fitted["FIELD"].split("+")
        if item["field"] != fitted["FIELD"]:
            compound_fitres_field_ids.append(cid)
        with np.load(BASE / "objectives" / f"objective_{cid}.npz", allow_pickle=False) as d:
            nominal = np.asarray(d["parameters_x0_x1_c_t0"], dtype=float)
            z_objective = float(d["zHEL"][0])
        assert nominal.shape == (4,) and np.all(np.isfinite(nominal))
        assert abs(z_objective - float(item["zHEL"])) < 1e-5
        fit_vector = np.array([float(fitted[k]) for k in ["x0", "x1", "c", "PKMJD"]])
        max_fitres_coordinate_difference = max(max_fitres_coordinate_difference, float(np.max(abs(nominal - fit_vector))))
        # FITRES prints the native fit coordinates after float32 storage;
        # objective NPZ preserves the unrounded values used by the exporter.
        fitres_rounding_tolerance = 0.5 * abs(np.spacing(fit_vector.astype(np.float32))).astype(float) + 1e-9
        assert np.all(abs(nominal - fit_vector) <= fitres_rounding_tolerance), cid

        parsed = {}
        for mode in ["observer", "sed"]:
            pair = pairmap[(cid, mode)]
            assert pair["paired_gate"] == "True"
            assert abs(float(pair["zHEL"]) - z_objective) < 1e-5
            base = np.array(json.loads(pair["baseline_theta"]), dtype=float)
            alt = np.array(json.loads(pair["alternative_theta"]), dtype=float)
            shift = np.array(json.loads(pair["actual_shift"]), dtype=float)
            assert base.shape == alt.shape == shift.shape == (4,)
            assert np.all(np.isfinite(base)) and np.all(np.isfinite(alt)) and np.all(np.isfinite(shift))
            err = float(np.max(abs(alt - base - shift)))
            max_shift_identity_error = max(max_shift_identity_error, err)
            assert err < 1e-10
            parsed[mode] = (base, alt)
        max_baseline_mode_difference = max(max_baseline_mode_difference, float(np.max(abs(parsed["observer"][0]-parsed["sed"][0]))))
        assert np.allclose(parsed["observer"][0], parsed["sed"][0], rtol=0, atol=1e-10)

        x = {
            "CID": cid, "field": item["field"], "zHEL": float(item["zHEL"]),
            "zHD": float(item["zHD"]), "zHEL_bin": zbin(float(item["zHEL"])),
            "released_x1": float(released["x1"]), "released_c": float(released["c"]),
            "native_x1": float(nominal[1]), "native_c": float(nominal[2]),
            "nominal_observed_x1": float(nominal[1] + parsed["observer"][0][1]),
            "nominal_observed_c": float(nominal[2] + parsed["observer"][0][2]),
            "observer_x1": float(nominal[1] + parsed["observer"][1][1]),
            "observer_c": float(nominal[2] + parsed["observer"][1][2]),
            "sed_x1": float(nominal[1] + parsed["sed"][1][1]),
            "sed_c": float(nominal[2] + parsed["sed"][1][2]),
            "released_x1ERR": float(released["x1ERR"]),
            "released_PKMJDERR": float(released["PKMJDERR"]),
            "released_cERR": float(released["cERR"]),
            "released_FITPROB": float(released["FITPROB"]),
        }
        for arm in ["released", "native", "nominal_observed", "observer", "sed"]:
            shape, color, joint = status(x[f"{arm}_x1"], x[f"{arm}_c"])
            x[f"{arm}_shape_pass"] = shape
            x[f"{arm}_color_pass"] = color
            x[f"{arm}_joint_pass"] = joint
            x[f"{arm}_failure_reason"] = failure_reason(shape, color)
        x["released_static_proxy_pass"] = (
            0 < x["released_x1ERR"] < 1
            and 0 < x["released_PKMJDERR"] < 2
            and 0 < x["released_cERR"] < 1.5
            and 0.001 < x["released_FITPROB"] < 1.1
            and 0.025 < x["zHD"] < 1.2
        )
        for arm in ["observer", "sed"]:
            x[f"native_to_{arm}"] = transition(x["native_joint_pass"], x[f"{arm}_joint_pass"])
        x["observer_to_sed"] = transition(x["observer_joint_pass"], x["sed_joint_pass"])
        rows.append(x)

    with (OUT / "per-object.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    result = summarize(rows)
    result["mapping_gates"] = {
        "max_abs_native_vs_fitres_coordinate": max_fitres_coordinate_difference,
        "max_abs_theta_shift_identity_error": max_shift_identity_error,
        "max_abs_observer_sed_nominal_theta_difference": max_baseline_mode_difference,
        "fitres_compound_field_ids": compound_fitres_field_ids,
        "resolved_paired_manifest_sha256": sha(RESOLVED / "manifest.json"),
    }
    result["scope"] = "Fixed 1020 previously selected DES objects; retrospective colour/shape-only eligibility under hypothetical nonlinear fitted means, no BBC or cosmology replay."
    (OUT / "result.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    (OUT / "executed_source.py").write_bytes(Path(__file__).read_bytes())
    manifest = {
        "protocol_sha256": sha(OUT / "protocol.md"),
        "freeze_sha256": sha(FREEZE),
        "executed_source_sha256": sha(OUT / "executed_source.py"),
        "result_sha256": sha(OUT / "result.json"),
        "per_object_sha256": sha(OUT / "per-object.csv"),
        "input_file_count": len(frozen["input_sha256"]),
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"objects": len(rows), "arms": result["arms"], "transitions": result["transitions"]}, indent=2))


if __name__ == "__main__":
    main()
