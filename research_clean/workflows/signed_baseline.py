"""Extract and summarize the ten-object signed author precursor baseline."""

import json
import numpy as np
from lib.paths import DATA
from lib.records import write_rows

DEFAULTS = {"cut_days": [180, 365]}


def read(path):
    columns = None
    rows = []
    for lineno, line in enumerate(path.read_text().splitlines(), 1):
        a = line.split()
        if not a:
            continue
        if a[0] == "VARLIST:":
            columns = a[1:]
        elif a[0] == "OBS:":
            if columns is None or len(columns) != len(a) - 1:
                raise ValueError("Malformed signed photometry")
            item = dict(zip(columns, a[1:]))
            band = item.get("FLT", item.get("BAND"))
            if band not in ["g", "r", "i", "z"]:
                continue
            row = {
                "MJD": float(item["MJD"]),
                "flux": float(item["FLUXCAL"]),
                "error": float(item["FLUXCALERR"]),
                "band": band,
                "source_line": lineno,
            }
            if (
                not all(np.isfinite(row[k]) for k in ["MJD", "flux", "error"])
                or row["error"] <= 0
            ):
                raise ValueError("Invalid signed row; no silent clipping")
            rows.append(row)
    return rows


def run(out, cfg):
    if not cfg["cut_days"] or any(x <= 0 for x in cfg["cut_days"]):
        raise ValueError("Pre-explosion cut must be positive")
    cohort = json.loads((DATA / "signed/cohort.json").read_text())
    groups = []
    ledger = []
    for cut in cfg["cut_days"]:
        for case in cohort:
            raw = read(DATA / "signed" / case["raw_file"])
            selected = [r for r in raw if r["MJD"] < case["peak_header"] - cut]
            ledger.extend({"CID": case["CID"], "cut_days": cut, **r} for r in selected)
            for band in ["g", "r", "i", "z"]:
                rows = [r for r in selected if r["band"] == band]
                if len(rows) < 3:
                    raise ValueError("Insufficient object-band support")
                flux = np.array([r["flux"] for r in rows])
                error = np.array([r["error"] for r in rows])
                inv = 1 / error**2
                baseline = float(np.sum(inv * flux) / inv.sum())
                q = float(np.sum(((flux - baseline) / error) ** 2))
                groups.append(
                    {
                        "CID": case["CID"],
                        "band": band,
                        "cut_days": cut,
                        "rows": len(rows),
                        "negative": int((flux < 0).sum()),
                        "baseline": baseline,
                        "baseline_conditional_se": float(1 / np.sqrt(inv.sum())),
                        "Q": q,
                        "dof": len(rows) - 1,
                    }
                )
    write_rows(out / "signed_rows.csv", ledger)
    write_rows(out / "object_band_baselines.csv", groups)
    summary = {}
    for cut in cfg["cut_days"]:
        rows = [r for r in groups if r["cut_days"] == cut]
        q = sum(x["Q"] for x in rows)
        dof = sum(x["dof"] for x in rows)
        summary[str(cut)] = {
            "rows": sum(x["rows"] for x in rows),
            "groups": len(rows),
            "negative": sum(x["negative"] for x in rows),
            "Q": q,
            "dof": dof,
            "Q_per_dof": q / dof,
        }
    return {
        "objects": len(cohort),
        "cuts": summary,
        "scope": "Signed precursor photometry before fixed observer-frame pre-explosion cuts. Baselines are estimated per object/band; no row clipping, PHOTFLAG cut, subtraction from fitted SN data or error rescaling. Q/dof is descriptive under the diagonal-error reference, not proof of independent Gaussian noise.",
    }
