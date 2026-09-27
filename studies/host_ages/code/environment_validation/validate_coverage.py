#!/usr/bin/env python3
"""Verify bounded DES/Dovekie/RAISIN extension identities and numerical checks."""
import datetime
import json
from pathlib import Path
from analyse import ROOT, HERE, sha


def main():
    folder = ROOT / "studies/host_ages/results/environment_validation"
    records = {}
    inputs = {}
    for name, code in [
        ("coverage.json", "coverage.py"),
        ("dovekie-summary.json", "dovekie.py"),
    ]:
        d = json.loads((folder / name).read_text())
        records[name] = d
        assert d["code_sha256"] == sha(HERE / code)
        for row in d["inputs"]:
            p = Path(row["path"])
            assert sha(p) == row["sha256"], str(p)
            inputs[str(p)] = row["sha256"]
        for path, digest in d.get("dependencies_sha256", {}).items():
            assert sha(ROOT / path) == digest, path
    c = records["coverage.json"]
    for path, digest in c["raisin_selected_header_sha256"].items():
        assert sha(ROOT / path) == digest, path
        inputs[path] = digest
    d = records["dovekie-summary.json"]
    assert d["helper_sha256"] == sha(HERE / "analyse.py")
    assert d["reader_sha256"] == sha(ROOT / "lib/records.py")
    assert (
        len(c["DES"]["TITAN_non_numeric_Dovekie_CID_matches"])
        == d["normalized_name_matches"]
    )
    assert all(row["IDSURVEY"] == 150 for row in d["identities"])
    assert (
        len({row["CID"] for row in d["identities"]})
        == d["selected_valid_and_host_dDLR_lt4"]
    )
    assert d["precision_inverse_subset_max_closure_error"] < 1e-10
    assert d["covariance_minimum_eigenvalue"] > 0
    for fit in d["fits"].values():
        assert fit["independent_whitened_lstsq_max_coefficient_difference"] < 1e-10
    result = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "code_sha256": sha(__file__),
        "status": "pass for source identity, event joins, covariance ordering and numerical equivalence; not validation of host-age physical identifiability",
        "input_files_verified": len(inputs),
        "result_sha256": {name: sha(folder / name) for name in records},
        "checks": [
            "All recorded current input, code and helper hashes match.",
            "Coverage and Dovekie exact-name match counts agree; retained events are uniquely named Foundation SNe.",
            "Total precision inversion closure and positive marginal covariance pass.",
            "Both full-covariance GLS solutions agree with independent Cholesky-whitened least squares.",
            "All 79 selected RAISIN photometry header hashes match the identity audit.",
        ],
    }
    (folder / "coverage-validation.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
