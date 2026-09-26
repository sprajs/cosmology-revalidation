"""Independent readback of saved DES cut table and source coordinates."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np

root = Path(__file__).resolve().parents[3]
out = Path(__file__).resolve().parent
rows = list(csv.DictReader((out / "per-object.csv").open()))
assert len(rows) == len({r["CID"] for r in rows}) == 1020
result = json.loads((out / "result.json").read_text())
paired = list(csv.DictReader((root / "runs/research_2026_09_26/sed_nonlinear_validation/resolved/paired-responses.csv").open()))
paired = {(r["CID"], r["mode"]): r for r in paired if r["target"] == "observed_flux"}
assert len(paired) == 2040

counts = Counter()
transitions = {"observer": Counter(), "sed": Counter()}
max_coordinate_error = 0.0
for row in rows:
    cid = row["CID"]
    with np.load(root / f"runs/research_2026_09_26/astra_design/validation1020/objectives/objective_{cid}.npz", allow_pickle=False) as data:
        nominal = data["parameters_x0_x1_c_t0"]
    for arm in ["native", "observer", "sed"]:
        if arm == "native":
            x1, c = nominal[1], nominal[2]
        else:
            theta = json.loads(paired[(cid, arm)]["alternative_theta"])
            x1, c = nominal[1] + theta[1], nominal[2] + theta[2]
        max_coordinate_error = max(max_coordinate_error, abs(x1 - float(row[f"{arm}_x1"])), abs(c - float(row[f"{arm}_c"])))
        passed = (-3 < x1 < 3) and (-0.3 < c < 0.3)
        counts[arm] += int(passed)
        if arm == "native":
            baseline = passed
        else:
            transitions[arm][("pass" if baseline else "fail") + "_to_" + ("pass" if passed else "fail")] += 1

assert max_coordinate_error < 1e-12
assert all(counts[a] == result["arms"][a]["joint_pass"] for a in counts)
assert all(dict(transitions[a]) == result["transitions"][f"native_to_{a}"] for a in transitions)
assert sum(v["n"] for v in result["by_field"].values()) == 1020
assert sum(v["n"] for v in result["by_zHEL_bin"].values()) == 1020
check = {
    "gate": True,
    "objects": 1020,
    "joint_pass": dict(counts),
    "baseline_transition_matrices": {a: dict(transitions[a]) for a in transitions},
    "max_recomputed_coordinate_error": max_coordinate_error,
    "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}
(out / "independent-verification.json").write_text(json.dumps(check, indent=2) + "\n")
print(json.dumps(check, indent=2))
