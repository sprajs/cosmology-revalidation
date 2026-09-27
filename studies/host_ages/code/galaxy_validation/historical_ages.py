#!/usr/bin/env python3
"""Paired original/updated age ordering against the same independent spectra."""

from pathlib import Path
import argparse, json, datetime
import numpy as np, pandas as pd
from astropy.io import ascii
from scipy.stats import spearmanr
from acquire import ROOT, WORK, OUT, sha
from match import ARCHIVE


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", type=Path, default=ARCHIVE)
    args = ap.parse_args()
    q = pd.read_csv(WORK / "host-spectrum-age-ledger.csv", dtype={"specObjID": str})
    q = q[q.index_quality].drop_duplicates(["sample", "physical_host"])
    gp = args.archive / "data/host_ages/gupta2011"
    g = ascii.read(
        gp / "table2.dat", format="cds", readme=str(gp / "ReadMe")
    ).to_pandas()
    parsed = WORK / "parsed-r19"
    parsed.mkdir(exist_ok=True)
    (parsed / "table7.dat").write_bytes((WORK / "R19-table7.dat").read_bytes())
    r = ascii.read(
        parsed / "table7.dat", format="cds", readme=str(WORK / "R19-ReadMe")
    ).to_pandas()
    out = {
        "status": "Post-primary-inspection extension; tests age rank ordering, not absolute Gyr calibration",
        "samples": {},
    }
    for sample, table, col in [("G11", g, "Age"), ("R19", r, "logAg")]:
        table["source_id"] = sample + ":" + table.SNID.astype(str)
        d = q[q["sample"] == sample].merge(
            table[["source_id", col]], on="source_id", validate="one_to_one"
        )
        results = {}
        for field, sign in [("d4000_n", 1), ("lick_hd_a_sub", -1)]:
            vals = d[[col, "sed_age", field]].dropna().to_numpy()
            rng = np.random.default_rng(2026092802)

            def statistic(a):
                return np.array(
                    [sign * spearmanr(a[:, i], a[:, 2]).statistic for i in [0, 1]]
                )

            point = statistic(vals)
            boot = np.array(
                [
                    statistic(vals[rng.integers(len(vals), size=len(vals))])
                    for _ in range(2000)
                ]
            )
            results[field] = {
                "n": len(vals),
                "orientation": sign,
                "original_age_oriented_rho": float(point[0]),
                "updated_age_oriented_rho": float(point[1]),
                "updated_minus_original_rho": float(point[1] - point[0]),
                "paired_host_bootstrap_95_delta_rho": np.quantile(
                    boot[:, 1] - boot[:, 0], [0.025, 0.975]
                ).tolist(),
            }
        out["samples"][sample] = results
    out["provenance"] = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "code_sha256": sha(Path(__file__)),
        "design_sha256": sha(Path(__file__).with_name("historical-age-design.json")),
        "inputs": {
            p.name: sha(p)
            for p in [
                gp / "table2.dat",
                gp / "ReadMe",
                WORK / "R19-table7.dat",
                WORK / "R19-ReadMe",
                WORK / "host-spectrum-age-ledger.csv",
            ]
        },
    }
    out["limitations"] = [
        "Original and updated ages use related photometry and stellar assumptions; spectra validate ranking, not an unbiased absolute-age scale or age slope calibration.",
        "R19 CDS labels logAg/logAl are internally misleading relative to their numerical age entries; this analysis uses only rank order and never exponentiates them.",
        "Central fibre versus global SED aperture mismatch and selected small SDSS subset remain; bootstrap differences are conditional on measured summaries.",
    ]
    (OUT / "historical-age-summary.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
