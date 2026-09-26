"""Independent, read-only source/asset check. No simulation or fitting."""
from __future__ import annotations

from collections import Counter
from pathlib import Path
import hashlib
import json
import re

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
BUNDLE = ROOT / "runs/research_2026_09_26/raisin_simulation_assets"
OUT = Path(__file__).resolve().parent
AUTHOR = BUNDLE / "author"
SIM_INPUT = ROOT / "runs/research_2026_09_26/raisin_sign_source/sim/inputs/DES/sim_DES_SNOOPY.input"
SRC = ROOT / "sources/repos/RickKessler__SNANA@v11_04k/src"
OIR = AUTHOR / "OIR.J19/OIR.INFO"
J22 = ROOT / "sources/repos/djones1040__RAISIN_DataRelease@a383c4b/model/OIR.J22/OIR.INFO"
BUNDLED_J19 = ROOT / "phase2/official/inputs/SNDATA_ROOT/models/OIR/OIR.J19/OIR.INFO"
LOGIC = ROOT / "phase2/official/inputs/SNDATA_ROOT/models/searcheff/SEARCHEFF_PIPELINE_LOGIC.DAT"
TREE = ROOT / "runs/research_2026_09_26/raisin_differential/mass-threshold/code-tree.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path, key: str) -> list[list[str]]:
    return [line.split()[1:] for line in path.read_text().splitlines()
            if line.lstrip().startswith(key + ":")]


def eff(path: Path) -> dict[str, list[tuple[float, float]]]:
    result: dict[str, list[tuple[float, float]]] = {}
    band = ""
    for line in path.read_text().splitlines():
        parts = line.split()
        if parts and parts[0] == "FILTER:":
            band = parts[1]
            assert band not in result
            result[band] = []
        elif parts and parts[0] == "SNR:":
            assert band
            result[band].append((float(parts[1]), float(parts[2])))
    return result


def main() -> None:
    acquisition = json.loads((BUNDLE / "acquisition.json").read_text())
    reported = json.loads((BUNDLE / "result.json").read_text())
    tree = json.loads(TREE.read_text())
    blobs = {x["path"]: x for x in tree["tree"] if x["type"] == "blob"}
    failures: list[str] = []

    assert acquisition["repository_commit"] == tree["sha"]
    assets = acquisition["files"]
    assert len(assets) == len({x["source_path"] for x in assets}) == 7
    blob_results = []
    for item in assets:
        path = ROOT / item["path"]
        data = path.read_bytes()
        got_blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        closure = (len(data) == item["bytes"] and sha(path) == item["sha256"]
                   and got_blob == item["git_blob"] == blobs[item["source_path"]]["sha"]
                   and len(data) == blobs[item["source_path"]]["size"])
        if not closure:
            failures.append("blob/size/hash mismatch: " + item["source_path"])
        blob_results.append({"source_path": item["source_path"], "sha256": sha(path),
                             "git_blob": got_blob, "bytes": len(data), "closed": closure})

    cfg = {}
    for raw in SIM_INPUT.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if ":" in line:
            key, val = line.split(":", 1)
            cfg[key.strip()] = val.strip()
    required_cfg = {
        "FLUXERRMODEL_FILE": "DES3YR_SIM_ERRORFUDGES.DAT",
        "GENMAG_SMEAR_MODELNAME": "OIR.J19",
        "SEARCHEFF_SPEC_FILE": "SEARCHEFF_SPEC_DES_Moller_G10_v7.DAT",
        "SEARCHEFF_PIPELINE_FILE": "SEARCHEFF_PIPELINE_RAISIN.DAT",
        "SEARCHEFF_PIPELINE_LOGIC_FILE": "SEARCHEFF_PIPELINE_LOGIC.DAT",
        "NEWMJD_DIF": "0.4",
        "APPLY_SEARCHEFF_OPT": "3",
    }
    for key, expected in required_cfg.items():
        actual = cfg.get(key, "")
        if (Path(actual).name if "/" in actual else actual) != expected:
            failures.append(f"input config {key}: {actual!r} != {expected!r}")

    errorfile = AUTHOR / "sim/inputs/DES/DES3YR_SIM_ERRORFUDGES.DAT"
    error_text = errorfile.read_text()
    maps = []
    for block in error_text.split("MAPNAME:")[1:]:
        head = block.splitlines()[0].strip()
        before_end = block.split("ENDMAP:", 1)[0]
        match = re.search(r"(?m)^BAND:\s*(\w+)\s+FIELD:\s*(\w+)", before_end)
        assert match
        values = [(float(a), float(b)) for a, b in re.findall(r"(?m)^ROW:\s*([\d.]+)\s+([\d.]+)", before_end)]
        maps.append({"name": head, "band": match.group(1), "field": match.group(2),
                     "rows": len(values), "min_scale": min(v[1] for v in values),
                     "max_scale": max(v[1] for v in values)})
    fieldgroups = {x[0]: x[1].split("+") for x in rows(errorfile, "DEFINE_FIELDGROUP")}
    if len(maps) != 8 or sum(m["rows"] for m in maps) != 71:
        failures.append("error map count/rows")
    if Counter((m["band"], m["field"]) for m in maps) != Counter(
            (b, f) for b in "griz" for f in ("DEEP", "SHALLOW")):
        failures.append("error map band/field coverage")
    error_range = [min(m["min_scale"] for m in maps), max(m["max_scale"] for m in maps)]
    if error_range != reported["error_scale_tabulated_range"]:
        failures.append("error range result mismatch")
    if fieldgroups != {"SHALLOW": ["E1", "E2", "S1", "S2", "C1", "C2", "X1", "X2"],
                       "DEEP": ["C3", "X3"]}:
        failures.append("error field-group mapping")

    des = AUTHOR / "sim/inputs/DES"
    pipe = eff(des / "SEARCHEFF_PIPELINE_RAISIN.DAT")
    ordinary = eff(des / "SEARCHEFF_PIPELINE_DES.DAT")
    pipe_check = {b: len(pipe[b]) == 33 and len(ordinary[b]) == 32
                  and pipe[b][:32] == ordinary[b] and pipe[b][-1] == (1_000_000.0, 1.0)
                  for b in "griz"}
    if set(pipe) != set("griz") or set(ordinary) != set("griz") or not all(pipe_check.values()):
        failures.append("pipeline mapping")
    nominal_spec = {float(a): float(b) for a, b in rows(des / "SEARCHEFF_SPEC_DES_Moller_G10_v7.DAT", "SPECEFF")}
    flat_spec = {float(a): float(b) for a, b in rows(AUTHOR / "sim/flatdist/DES/SEARCHEFF_SPEC_DES_Moller_G10_v7.DAT", "SPECEFF")}
    spec_check = (len(nominal_spec) == 99 and len(flat_spec) == 79
                  and all(nominal_spec[x] == y for x, y in flat_spec.items()))
    if not spec_check:
        failures.append("spectroscopic mapping")

    oir_bytes_equal = OIR.read_bytes() == BUNDLED_J19.read_bytes() == J22.read_bytes()
    if not oir_bytes_equal:
        failures.append("OIR author/bundled/release byte identity")
    sigma = np.array([float(x) for x in rows(OIR, "COLOR_SIGMA")[0]])
    corr = np.array([[float(x) for x in row] for row in rows(OIR, "COLOR_CORMAT")])
    scale = float(rows(OIR, "COLOR_SIGMA_SCALE")[0][0])
    cov = 1.3 * corr * (np.outer(sigma, sigma) + np.eye(7) * 1e-9)
    eig = np.linalg.eigvalsh(cov)
    if sigma.shape != (7,) or corr.shape != (7, 7) or scale != 1.0:
        failures.append("OIR schema")
    if not np.allclose(eig, reported["OIR_native_node_covariance_eigenvalues"], atol=1e-15, rtol=0):
        failures.append("OIR covariance result mismatch")

    flux_source = (SRC / "sntools_fluxErrModels.c").read_text()
    smear_source = (SRC / "sntools_genSmear.c").read_text()
    oir_init = smear_source.split("void init_genSmear_OIR", 1)[1].split("void get_genSmear_OIR", 1)[0]
    oir_read = smear_source.split("void read_OIR_INFO", 1)[1].split("void sort_OIR_BANDS", 1)[0]
    source_checks = {
        "map_default_both": "FLUXERRMAP[NMAP].MASK_APPLY   = 3" in flux_source,
        "map_applies_true": "MASK_APPLY & MASK_APPLY_SIM_FLUXERRMAP" in flux_source,
        "map_applies_reported": "MASK_APPLY & MASK_APPLY_DATA_FLUXERRMAP" in flux_source,
        "errscale_multiplies": "FLUXERR_OUT = (fluxErr * errModelVal)" in flux_source,
        "oir_diag_fudge": "CC += COV_DIAG_FUDGE" in oir_init,
        "oir_cov_scale": "COVAR1[N-1] = COVAR2[i][j] * COV_SCALE" in oir_init,
        "oir_sigma_scale_parsed": "GENSMEAR_OIR.COLOR_SIGMA_SCALE" in oir_read,
        "oir_sigma_scale_not_used_in_init": "GENSMEAR_OIR.COLOR_SIGMA_SCALE" not in oir_init,
    }
    if not all(source_checks.values()):
        failures.append("source-rule checks")
    logic_lines = [x for x in LOGIC.read_text().splitlines() if x.startswith("DES:")]
    if logic_lines != reported["DES_bundled_logic"]:
        failures.append("bundled detection logic")

    root_script = ROOT / "scripts/research_2026_09_26/raisin_simulation_assets.py"
    manifest = json.loads((BUNDLE / "manifest.json").read_text())["files_sha256"]
    manifest_failures = [p for p, digest in manifest.items() if sha(BUNDLE / p) != digest]
    if manifest_failures:
        failures.append("run manifest mismatch")
    if sha(root_script) != reported["inputs_sha256"][str(root_script.relative_to(ROOT))]:
        failures.append("reported source hash mismatch")
    if root_script.read_bytes() != (BUNDLE / "executed_source.py").read_bytes():
        failures.append("executed source snapshot mismatch")

    result = {
        "scope": "Independent read-only source/asset review; no simulation or historical execution claim",
        "pass": not failures,
        "failures": failures,
        "seven_blobs": blob_results,
        "nominal_config": {k: cfg.get(k) for k in required_cfg},
        "error_fieldgroups": fieldgroups,
        "error_maps": maps,
        "error_map_rows": sum(m["rows"] for m in maps),
        "error_scale_range": error_range,
        "pipeline_exact_shared32_plus_endpoint": pipe_check,
        "spectroscopic_rows": {"nominal": len(nominal_spec), "flatdist": len(flat_spec),
                               "shared": len(set(nominal_spec) & set(flat_spec)), "match": spec_check},
        "oir_byte_identity": {"author_J19": sha(OIR), "bundled_J19": sha(BUNDLED_J19),
                              "release_J22": sha(J22), "same_bytes": oir_bytes_equal},
        "oir_covariance_eigenvalues": eig.tolist(),
        "source_rule_checks": source_checks,
        "bundled_logic_DES": logic_lines,
        "root_manifest_entries": len(manifest), "root_manifest_failures": manifest_failures,
        "input_sha256": {str(p.relative_to(ROOT)): sha(p) for p in
                         [TREE, BUNDLE / "acquisition.json", BUNDLE / "result.json",
                          BUNDLE / "manifest.json", root_script, SRC / "sntools_fluxErrModels.c",
                          SRC / "sntools_genSmear.c", SIM_INPUT, OIR, BUNDLED_J19, J22, LOGIC]},
        "review_script_sha256": sha(Path(__file__)),
    }
    OUT.mkdir(exist_ok=True)
    (OUT / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"pass": result["pass"], "failures": failures,
                      "seven_blob_count": len(blob_results), "maps": len(maps)}, indent=2))


if __name__ == "__main__":
    main()
