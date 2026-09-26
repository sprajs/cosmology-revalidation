#!/usr/bin/env python3
"""Post-score descriptive Γ/mean decomposition from saved sufficient arrays only."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROTOCOL = json.loads((HERE / "report-protocol.json").read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


for name, digest in PROTOCOL["inputs"].items():
    if sha(HERE / name) != digest:
        raise RuntimeError(f"saved score input changed: {name}")

with np.load(HERE / "read-matrix.npz", allow_pickle=False) as arrays:
    gamma = arrays["quadrant_gamma"]
    mean = arrays["quadrant_mean"]
    quadrants = list(arrays["quadrant_names"].astype(str))
with np.load(HERE.parent / "pixel_design/variance-diagnostic-review/calwf3_read_matrix/dark-quality-design/fixed-operators.npz", allow_pickle=False) as op:
    h = op["h"][0]
all_rows = list(csv.DictReader((HERE / "pair-quadrant-power.csv").open()))
pairs = list(dict.fromkeys(row["pair"] for row in all_rows))
assert len(pairs) == gamma.shape[0] == 8 and quadrants == list("BCAD")

rows = []
for ip, pair in enumerate(pairs):
    for iq, quadrant in enumerate(quadrants):
        source = next(row for row in all_rows if row["pair"] == pair
                      and row["quadrant"] == quadrant and row["scope"] == "masked"
                      and row["power"] == "0.0")
        g = float(h @ gamma[ip, 0, iq] @ h)
        signed = float(h @ mean[ip, 0, iq])
        coherent = signed * signed / 2
        centered = g - coherent
        recorded = float(source["dn2_per_nominal_s2"])
        if abs(g - recorded) > 1e-12 or centered < -1e-12:
            raise RuntimeError(f"Gamma/mean closure failed for {pair}/{quadrant}")
        rows.append({"pair": pair, "visit": source["visit"], "kind": source["kind"],
                     "quadrant": quadrant, "n_fixed_mask_pixels": int(source["n"]),
                     "half_second_moment_dn2_per_s2": g,
                     "signed_repeat_mean_slope_dn_per_s": signed,
                     "half_squared_repeat_mean_dn2_per_s2": coherent,
                     "centered_half_variance_dn2_per_s2": centered,
                     "coherent_fraction": coherent / g if g else None,
                     "saved_csv_gap": g - recorded})

with (HERE / "decomposition.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

quad_summary = {}
for q in quadrants:
    subset = [x for x in rows if x["quadrant"] == q and x["kind"] == "primary_all14"]
    quad_summary[q] = {"primary_pair_count": len(subset),
                       "second_moment_min": min(x["half_second_moment_dn2_per_s2"] for x in subset),
                       "second_moment_max": max(x["half_second_moment_dn2_per_s2"] for x in subset),
                       "coherent_fraction_min": min(x["coherent_fraction"] for x in subset),
                       "coherent_fraction_max": max(x["coherent_fraction"] for x in subset)}

digital = list(csv.DictReader((HERE / "raw-dq-digital-contributions.csv").open()))
if len(digital) != 16:
    raise RuntimeError("8-pair digital/DQ contribution ledger incomplete")
result = {"status": "saved-score arithmetic only; no RAW SCI read",
          "report_protocol_sha256": sha(HERE / "report-protocol.json"),
          "pairs": pairs, "n_decomposition_rows": len(rows),
          "max_saved_csv_gap": max(abs(x["saved_csv_gap"]) for x in rows),
          "quadrant_primary_ranges": quad_summary,
          "normal_search": [x for x in rows if x["pair"] == "normal_search"],
          "digital_and_DQ_contributions": digital}
(HERE / "descriptive-summary.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({"rows": len(rows), "max_gap": result["max_saved_csv_gap"],
                  "normal_search_rows": len(result["normal_search"])}, indent=2))
