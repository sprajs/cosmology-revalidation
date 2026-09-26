#!/usr/bin/env python3
"""Check restored source identities, workspace isolation and selected calculations.

This is not an end-to-end validation of all native, classifier or BayeSN studies.
"""
import ast
from datetime import datetime, timezone
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import warnings

ROOT = Path(__file__).resolve().parents[1]


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def main():
    manager = module(ROOT / "studies/manage.py", "study_manager")
    result = {"scope": __doc__.strip(), "identities": manager.verify()}
    paths = sorted({ROOT / r["path"] for r in manager.retained()})
    with warnings.catch_warnings(record=True) as notices:
        warnings.simplefilter("always", SyntaxWarning)
        for path in paths:
            if path.suffix == ".py":
                ast.parse(path.read_bytes(), filename=str(path))
            if path.suffix == ".sh":
                subprocess.run(["bash", "-n", str(path)], check=True)
    result["python_syntax_checked"] = sum(p.suffix == ".py" for p in paths)
    result["shell_syntax_checked"] = sum(p.suffix == ".sh" for p in paths)
    result["historical_string_escape_warnings"] = len(notices)
    name = "study-check-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%f")
    result["preparation"] = manager.prepare(name)
    work = manager.workspace(name)
    for invalid in ["../escape", "/tmp/escape"]:
        try:
            manager.within(work, invalid)
        except ValueError:
            pass
        else:
            raise AssertionError("Workspace path escaped")
    try:
        manager.prepare(name)
    except FileExistsError:
        result["existing_workspace_refused"] = True
    else:
        raise AssertionError("Existing workspace accepted")
    for number in (0, 1):
        manager.fetch(name, f"phase2/official/portable_pilot/assets/salt3_template_{number}.dat.gz", False)
    subprocess.run([sys.executable, str(work / "scripts/physics_audit/luminosity_check.py")],
                   cwd=work, check=True, stdout=subprocess.DEVNULL)
    luminosity = json.loads((work / "runs/physics_audit/luminosity-check.json").read_text())
    assert luminosity["passed"]
    result["luminosity"] = luminosity
    # Exercise a restored acquisition/preparation boundary on a real frozen RAW:
    # 16 reads -> first 8 reads, then an 8-read -> 8-read exact identity control.
    raw = next(r for r in manager.read("inputs.json")["files"] if r["path"].endswith("idp245t7q_raw.fits"))
    manager.fetch(name, raw["source_path"], False)
    prefix = module(work / "scripts/research_2026_09_26/prepare_calwf3_dark_prefix.py", "prefix")
    first = work / "first-eight.fits"
    control = work / "first-eight-control.fits"
    extraction = prefix.prefix(work / raw["source_path"], first)
    identity = prefix.prefix(first, control)
    assert extraction["source_nsamp"] == 16 and identity["null_whole_file_exact"]
    result["raw_prefix"] = {"source_sha256": extraction["source_sha256"],
                            "output_sha256": extraction["target_sha256"],
                            "read_data_blocks_verified": len(extraction["mapping"]),
                            "eight_read_identity_control": identity["null_whole_file_exact"]}
    # Corrupt only the isolated copy and verify it is rejected, without touching
    # the frozen input or attempting to redownload over the changed file.
    copy = work / raw["source_path"]
    with copy.open("r+b") as stream:
        stream.write(b"X")
    try:
        manager.fetch(name, raw["source_path"], False)
    except ValueError:
        result["changed_input_refused"] = True
    else:
        raise AssertionError("Changed input accepted")
    assert manager.digest(ROOT / raw["path"]) == raw["sha256"]
    result["frozen_input_unchanged"] = True
    shutil.rmtree(work)
    result["preparation"]["workspace_retained"] = False
    result["created_utc"] = datetime.now(timezone.utc).isoformat()
    result["status"] = "passed"
    result["code_sha256"] = {str(p.relative_to(ROOT)): manager.digest(p)
                             for p in [Path(__file__), ROOT / "studies/manage.py"]}
    out = ROOT / "validation/reports/studies.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: result[k] for k in ["status", "identities", "python_syntax_checked", "shell_syntax_checked", "raw_prefix"]}))


if __name__ == "__main__":
    main()
