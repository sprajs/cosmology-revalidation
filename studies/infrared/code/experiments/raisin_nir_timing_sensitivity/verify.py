"""Independent input-difference and fixed-peak closure after timing runs."""
from pathlib import Path
import csv
import hashlib
import json
import re
import struct

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
REL = ROOT / "sources/repos/djones1040__RAISIN_DataRelease@a383c4b/photometry/RAISIN/DES_RAISIN"
CASES = ("baseline", "baseline_copy", "minus1", "minus0p5", "plus0p5", "plus1", "author_t0")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fitres(path):
    cols = None; out = {}
    for line in path.read_text().splitlines():
        if line.startswith("VARNAMES:"):
            cols = line.split()[1:]
        elif line.startswith("SN:"):
            vals = line.split()[1:]
            assert cols and len(vals) == len(cols)
            row = dict(zip(cols, vals))
            out[row["CID"]] = row
    return out


def main():
    protocol = json.loads((OUT / "protocol-shortvpec.json").read_text())
    cohort = protocol["cohort"]
    with (OUT / "header-ledger.csv").open(newline="") as f:
        ledger = {(x["condition"], x["CID"]): x for x in csv.DictReader(f)}
    failures = []
    peak_differences = []
    for case in CASES:
        fit = fitres(OUT / "fits-shortvpec" / case / "fit.FITRES.TEXT")
        if set(fit) != set(cohort):
            failures.append(case + ": CID set")
        for cid in cohort:
            source = (REL / f"{cid}.snana.dat").read_text()
            variant = (OUT / "data" / case / "DES_RAISIN" / f"{cid}.snana.dat").read_text()
            source_lines, variant_lines = source.splitlines(keepends=True), variant.splitlines(keepends=True)
            diffs = [(a, b) for a, b in zip(source_lines, variant_lines) if a != b]
            if len(source_lines) != len(variant_lines) or any(not a.startswith("PEAKMJD:") or not b.startswith("PEAKMJD:") for a, b in diffs):
                failures.append(f"{case}/{cid}: nonheader source change")
            expected_diffs = 0 if case in ("baseline", "baseline_copy") else 2
            if len(diffs) != expected_diffs:
                failures.append(f"{case}/{cid}: expected {expected_diffs} changed header lines, got {len(diffs)}")
            peaks = re.findall(r"(?m)^PEAKMJD:\s*([^\s]+)", variant)
            if len(peaks) != 2 or peaks[0] != peaks[1] or peaks[0] != ledger[(case, cid)]["new_header"]:
                failures.append(f"{case}/{cid}: header ledger mismatch")
            expected32 = struct.unpack("f", struct.pack("f", float(peaks[0])))[0]
            printed = float(fit[cid]["PKMJD"])
            if abs(printed - expected32) > 5.1e-5 or float(fit[cid]["PKMJDERR"]) != 0:
                failures.append(f"{case}/{cid}: fixed peak != float32 header")
            peak_differences.append(abs(printed - expected32))
    if (OUT / "fits-shortvpec/baseline/fit.FITRES.TEXT").read_bytes() != (OUT / "fits-shortvpec/baseline_copy/fit.FITRES.TEXT").read_bytes():
        failures.append("baseline copy FITRES byte mismatch")
    result = {
        "pass": not failures, "failures": failures, "n_input_checks": len(CASES) * len(cohort),
        "n_shifted_header_only": (len(CASES) - 2) * len(cohort),
        "max_abs_printed_peak_minus_float32_header_day": max(peak_differences),
        "baseline_copy_FITRES_byte_identical": sha(OUT / "fits-shortvpec/baseline/fit.FITRES.TEXT") == sha(OUT / "fits-shortvpec/baseline_copy/fit.FITRES.TEXT"),
        "protocol_sha256": sha(OUT / "protocol-shortvpec.json"),
        "checker_sha256": sha(Path(__file__)),
    }
    (OUT / "independent-verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
