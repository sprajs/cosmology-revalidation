"""Summarize frozen NIR header shifts; does not launch a fit."""
from __future__ import annotations

from collections import Counter
from pathlib import Path
import csv
import hashlib
import json
import statistics

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/research_2026_09_26/raisin_nir_timing_sensitivity"
PROTO = OUT / "protocol-shortvpec.json"
CONDITIONS = ("baseline", "minus1", "minus0p5", "plus0p5", "plus1", "author_t0")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table(path: Path, key: str) -> dict[str, dict[str, str]]:
    cols = None
    out = {}
    for line in path.read_text().splitlines():
        if line.startswith("VARNAMES:"):
            cols = line.split()[1:]
        elif line.startswith(key + ":"):
            assert cols is not None
            vals = line.split()[1:]
            assert len(vals) == len(cols)
            row = dict(zip(cols, vals))
            assert row["CID"] not in out
            out[row["CID"]] = row
    return out


def accepted(path: Path) -> dict[str, Counter]:
    cols = None
    out = {}
    for line in path.read_text().splitlines():
        if line.startswith("VARNAMES:"):
            cols = line.split()[1:]
        elif line.startswith("OBS:"):
            assert cols is not None
            vals = line.split()[1:]
            assert len(vals) == len(cols)
            row = dict(zip(cols, vals))
            if row["DATAFLAG"] == "1":
                key = tuple(row[x] for x in ("MJD", "BAND", "FLUXCAL", "FLUXCAL_ERR"))
                out.setdefault(row["CID"], Counter())[key] += 1
    return out


def main() -> None:
    protocol = json.loads(PROTO.read_text())
    gate = json.loads((OUT / "baseline-gate.json").read_text())
    assert gate["pass"] and gate["protocol_sha256"] == sha(PROTO)
    cohort = protocol["cohort"]
    with (OUT / "header-ledger.csv").open(newline="") as f:
        headers = {(row["condition"], row["CID"]): row for row in csv.DictReader(f)}
    outputs = {}
    input_hashes = {str(PROTO.relative_to(ROOT)): sha(PROTO),
                    str((OUT / "baseline-gate.json").relative_to(ROOT)): sha(OUT / "baseline-gate.json"),
                    str((OUT / "header-ledger.csv").relative_to(ROOT)): sha(OUT / "header-ledger.csv"),
                    str(Path(__file__).relative_to(ROOT)): sha(Path(__file__))}
    for condition in CONDITIONS:
        work = OUT / "fits-shortvpec" / condition
        execution = json.loads((work / "execution.json").read_text())
        assert execution["returncode"] == 0 and execution["protocol_sha256"] == sha(PROTO)
        fit_path, plot_path = work / "fit.FITRES.TEXT", work / "fit.LCPLOT.TEXT"
        fit = table(fit_path, "SN")
        plot = accepted(plot_path)
        assert set(fit) == set(cohort)
        outputs[condition] = (fit, plot)
        for path in (work / "execution.json", fit_path, plot_path, work / "fit.log"):
            input_hashes[str(path.relative_to(ROOT))] = sha(path)
    baseline_fit, baseline_plot = outputs["baseline"]
    ledger = []
    failures = []
    for condition in CONDITIONS[1:]:
        fit, plot = outputs[condition]
        for cid in cohort:
            row, ref = fit[cid], baseline_fit[cid]
            header = headers[(condition, cid)]
            accepted_here, accepted_ref = plot.get(cid, Counter()), baseline_plot.get(cid, Counter())
            if row["ERRFLAG_FIT"] != "0":
                failures.append(f"{condition}/{cid}: ERRFLAG {row['ERRFLAG_FIT']}")
            if not (float(row["PKMJDERR"]) == 0 and float(row["STRETCH"]) == 1 and float(row["AV"]) == 0 and float(row["RV"]) == 1.518):
                failures.append(f"{condition}/{cid}: fixed parameter mismatch")
            if float(row["NDOF"]) != sum(accepted_here.values()) - 1:
                failures.append(f"{condition}/{cid}: NDOF acceptance mismatch")
            ledger.append({
                "condition": condition, "CID": cid,
                "header_original": header["original_header"], "header_new": header["new_header"],
                "header_delta_day": header["delta_day"],
                "baseline_peak": ref["PKMJD"], "fit_peak": row["PKMJD"],
                "baseline_DLMAG": ref["DLMAG"], "fit_DLMAG": row["DLMAG"],
                "delta_DLMAG": float(row["DLMAG"]) - float(ref["DLMAG"]),
                "baseline_DLMAGERR": ref["DLMAGERR"], "fit_DLMAGERR": row["DLMAGERR"],
                "baseline_FITCHI2": ref["FITCHI2"], "fit_FITCHI2": row["FITCHI2"],
                "baseline_NDOF": ref["NDOF"], "fit_NDOF": row["NDOF"],
                "baseline_accepted": sum(accepted_ref.values()), "fit_accepted": sum(accepted_here.values()),
                "accepted_multiset_equal": accepted_here == accepted_ref,
                "accepted_removed": sum((accepted_ref - accepted_here).values()),
                "accepted_added": sum((accepted_here - accepted_ref).values()),
                "ERRFLAG_FIT": row["ERRFLAG_FIT"],
            })
    path = OUT / "per-object-response.csv"
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=ledger[0].keys())
        writer.writeheader(); writer.writerows(ledger)
    summary = {}
    for condition in CONDITIONS[1:]:
        rows = [x for x in ledger if x["condition"] == condition]
        changes = [x["delta_DLMAG"] for x in rows]
        summary[condition] = {
            "n": len(rows), "header_delta_range_day": [min(float(x["header_delta_day"]) for x in rows),
                                                  max(float(x["header_delta_day"]) for x in rows)],
            "mean_delta_DLMAG": statistics.mean(changes),
            "median_delta_DLMAG": statistics.median(changes),
            "min_delta_DLMAG": min(changes), "max_delta_DLMAG": max(changes),
            "max_abs_delta_DLMAG": max(map(abs, changes)),
            "accepted_mask_changed_count": sum(not x["accepted_multiset_equal"] for x in rows),
            "accepted_removed_total": sum(x["accepted_removed"] for x in rows),
            "accepted_added_total": sum(x["accepted_added"] for x in rows),
            "NDOF_changed_count": sum(x["baseline_NDOF"] != x["fit_NDOF"] for x in rows),
            "ERRFLAG_nonzero_count": sum(x["ERRFLAG_FIT"] != "0" for x in rows),
        }
    result = {"scope": "Conditional ten-object timing sensitivity after frozen baseline gate; no inferred real correction",
              "pass": not failures, "failures": failures, "protocol_sha256": sha(PROTO),
              "baseline_gate_sha256": sha(OUT / "baseline-gate.json"),
              "per_object_sha256": sha(path), "summary": summary,
              "inputs_sha256": input_hashes}
    (OUT / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"pass": result["pass"], "failures": failures, "summary": summary}, indent=2))


if __name__ == "__main__":
    main()
