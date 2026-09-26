"""Frozen, gated ten-object NIR header-timing sensitivity (SNANA v11_03c)."""
from __future__ import annotations

from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/research_2026_09_26/raisin_nir_timing_sensitivity"
REL = ROOT / "sources/repos/djones1040__RAISIN_DataRelease@a383c4b"
PHOT = REL / "photometry/RAISIN/DES_RAISIN"
NML = REL / "lcfitting/REFAC_DES_RAISIN_nir_sys.nml"
AUTHOR = ROOT / "runs/research_2026_09_26/raisin_nir_inheritance_review/author"
COHORT = ROOT / "runs/research_2026_09_26/raisin_historical_native/fit-protocol-shortprefix.json"
BUILD = OUT / "SNANA-v11_03c-source"
BIN = BUILD / "bin/snlc_fit.exe"
PROTO = OUT / "protocol-shortvpec.json"
ORIGINAL_PROTO = OUT / "protocol.json"
CONDITIONS = {"baseline": 0.0, "minus1": -1.0, "minus0p5": -0.5,
              "plus0p5": 0.5, "plus1": 1.0, "author_t0": None}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path: Path, obj: object) -> None:
    path.write_text(json.dumps(obj, indent=2) + "\n")


def fitres(path: Path) -> dict[str, dict[str, str]]:
    cols = None
    result = {}
    for line in path.read_text().splitlines():
        if line.startswith("VARNAMES:"):
            cols = line.split()[1:]
        if line.startswith("SN:"):
            assert cols is not None
            vals = line.split()[1:]
            assert len(vals) == len(cols), (path, line)
            row = dict(zip(cols, vals))
            assert row["CID"] not in result
            result[row["CID"]] = row
    return result


def accepted_counts(path: Path) -> dict[str, int]:
    cols = None
    result = {}
    for line in path.read_text().splitlines():
        if line.startswith("VARNAMES:"):
            cols = line.split()[1:]
        elif line.startswith("OBS:"):
            assert cols is not None
            vals = line.split()[1:]
            assert len(vals) == len(cols)
            row = dict(zip(cols, vals))
            if row["DATAFLAG"] == "1":
                result[row["CID"]] = result.get(row["CID"], 0) + 1
    return result


def cohort_and_t0() -> tuple[list[str], dict[str, str]]:
    cohort = json.loads(COHORT.read_text())["cohort"]
    assert len(cohort) == len(set(cohort)) == 10
    t0 = {}
    for line in (AUTHOR / "data/raisin_t0.txt").read_text().splitlines():
        bits = line.split()
        if bits and bits[0] in cohort:
            assert bits[0] not in t0
            t0[bits[0]] = bits[1]
    assert set(t0) == set(cohort)
    return cohort, t0


def set_value(text: str, key: str, value: str) -> str:
    pattern = rf"(?m)^(\s*{re.escape(key)}\s*=).*$"
    new, count = re.subn(pattern, rf"\1 {value}", text)
    assert count == 1, (key, count)
    return new


def prepare() -> None:
    assert not ORIGINAL_PROTO.exists(), "Never overwrite a frozen protocol"
    out_source = OUT / "SNANA-v11_03c-source"
    commit = subprocess.check_output(["git", "-C", str(out_source), "rev-parse", "HEAD"], text=True).strip()
    assert commit == "06f2ccfdbf99d62c23dec66d0bf7b7e9452604d2"
    assert BIN.exists()
    cohort, t0 = cohort_and_t0()
    original = NML.read_text().split("&SNLCINP", 1)[1]
    original = "  &SNLCINP" + original
    base_files = [Path(__file__), ROOT / "scripts/research_2026_09_26/build_raisin_nir_v11_03c.sh",
                  OUT / "build.log", OUT / "build-generated-source-diff.patch", BIN, NML, COHORT, AUTHOR / "data/raisin_t0.txt",
                  AUTHOR / "output/fit_nir/DES_RAISIN.FITRES.TEXT",
                  REL / "kcor/kcor_DES_NIR.fits", REL / "model/snoopy.B18/SALT2.INFO",
                  REL / "vpec/vpec_baseline_raisin.list", PHOT / "DES_RAISIN.README"]
    # Model grid is a directory; hash each file, not just a directory name.
    base_files = [x for x in base_files if x.is_file()]
    base_files.extend(sorted((REL / "model/snoopy.B18").rglob("*")))
    base_files = [x for x in base_files if x.is_file()]
    source_files = [x for x in out_source.rglob("src/*") if x.is_file() and x.suffix in (".car", ".c", ".h")]
    # The git commit and root src tree identify the full source; retain direct source-file hash rollup.
    source_file_hashes = {str(x.relative_to(out_source)): sha(x) for x in sorted(source_files)}
    modified = subprocess.check_output(["git", "-C", str(out_source), "diff", "--name-only"], text=True).splitlines()
    assert modified == ["src/sntools_output.h"], modified
    assert subprocess.check_output(["git", "-C", str(out_source), "diff", "--", "src/sntools_output.h"]) == (OUT / "build-generated-source-diff.patch").read_bytes()
    rows = []
    jobs = []
    for condition, delta in [("baseline", 0.0), ("baseline_copy", 0.0),
                             ("minus1", -1.0), ("minus0p5", -0.5),
                             ("plus0p5", 0.5), ("plus1", 1.0),
                             ("author_t0", None)]:
        data_version = OUT / "data" / condition / "DES_RAISIN"
        work = OUT / "fits" / condition
        data_version.mkdir(parents=True, exist_ok=False)
        work.mkdir(parents=True, exist_ok=False)
        shutil.copyfile(PHOT / "DES_RAISIN.README", data_version / "DES_RAISIN.README")
        (data_version / "DES_RAISIN.LIST").write_text("".join(f"{cid}.snana.dat\n" for cid in cohort))
        for cid in cohort:
            raw = (PHOT / f"{cid}.snana.dat").read_text()
            found = re.findall(r"(?m)^PEAKMJD:\s*([^\s]+)", raw)
            assert len(found) == 2 and found[0] == found[1], (cid, found)
            new_peak = f"{float(found[0]) + delta:.1f}" if delta is not None else t0[cid]
            changed, count = re.subn(r"(?m)^(PEAKMJD:\s*)[^\s]+", lambda m: m.group(1) + new_peak, raw)
            assert count == 2
            (data_version / f"{cid}.snana.dat").write_text(changed)
            rows.append({"condition": condition, "CID": cid, "original_header": found[0],
                         "new_header": new_peak, "delta_day": float(new_peak) - float(found[0]),
                         "phot_source_sha256": sha(PHOT / f"{cid}.snana.dat"),
                         "phot_input_sha256": sha(data_version / f"{cid}.snana.dat")})
        text = original
        substitutions = {
            "PRIVATE_DATA_PATH": f"'{data_version.parent}'",
            "VERSION_PHOTOMETRY": "'DES_RAISIN'",
            "KCOR_FILE": f"'{REL / 'kcor/kcor_DES_NIR.fits'}'",
            "TEXTFILE_PREFIX": "'fit'",
            "SNTABLE_LIST": "'FITRES(text:key) LCPLOT(text:key)'",
            "HEADER_OVERRIDE_FILE": f"'{REL / 'vpec/vpec_baseline_raisin.list'}'",
            "FITMODEL_NAME": f"'{REL / 'model/snoopy.B18'}'",
        }
        for key, val in substitutions.items():
            text = set_value(text, key, val)
        nml = work / "fit.nml"
        nml.write_text(text)
        jobs.append({"condition": condition, "nml": str(nml.relative_to(ROOT)),
                     "private_data_path": str(data_version.parent.relative_to(ROOT)),
                     "nml_sha256": sha(nml)})
    with (OUT / "header-ledger.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader(); writer.writerows(rows)
    p = {
        "state": "frozen before any NIR fit from this source clone",
        "purpose": "Conditional response of v11_03c JH-only distances and accepted masks to fixed photometry-header timing in ten DES16 objects. No actual peak bias or cosmology correction inferred.",
        "source_tag": "v11_03c", "source_commit": commit,
        "source_src_tree": subprocess.check_output(["git", "-C", str(out_source), "rev-parse", "HEAD:src"], text=True).strip(),
        "build_modified_tracked_files": modified,
        "source_file_sha256": source_file_hashes,
        "cohort": cohort,
        "conditions": CONDITIONS,
        "t0_rule": "Use unrounded author data/raisin_t0.txt value in each PEAKMJD header, not the SIMLIB's one-decimal rendering.",
        "baseline": "Released DES_RAISIN JH flux/error rows and photometry PEAKMJD headers, released NIR NML settings (path/output-format substitutions only), v11_03c source, project-local modern compiler/libs.",
        "engineering_gate": "Run baseline and baseline_copy with byte-identical ten-file inputs in isolated directories. Require identical FITRES SN rows and LCPLOT bytes, ten ERRFLAG_FIT=0, fixed PKMJD/stretch/AV/RV, and NDOF closure.",
        "author_gate": "All ten baseline DLMAG differ by at most 0.001 mag from archived raw author NIR FITRES and NDOF match exactly. If any fail, stop sensitivity fits and diagnose without retuning.",
        "sensitivity": "Only after both gates: fixed-header shifts -1,-0.5,+0.5,+1 observer days and separate unrounded author raisin_t0; preserve identical JH epoch flux/error source rows and fit settings. Record DLMAG, PKMJD, FITCHI2, NDOF, ERRFLAG and accepted LCPLOT epoch multisets, including any failures.",
        "start_check": "Not separately varied: although INIVAL_DLMAG exists, source may reinitialize amplitude. With fixed peak/shape/AV/RV, report table and accepted-mask closure; do not rank across evolving C or assume a universal common objective.",
        "limitations": "No inferred true peak offset, regenerated bias correction, selection closure, full optical-header creator, or historical executable/toolchain identity.",
        "jobs": jobs,
        "inputs_sha256": {str(x.relative_to(ROOT)): sha(x) for x in base_files},
        "header_ledger_sha256": sha(OUT / "header-ledger.csv"),
    }
    dump(ORIGINAL_PROTO, p)
    (OUT / "executed-preparation-source.py").write_bytes(Path(__file__).read_bytes())
    print(json.dumps({"protocol_sha256": sha(ORIGINAL_PROTO), "jobs": len(jobs), "cohort": len(cohort)}))


def amend_short_vpec() -> None:
    assert not PROTO.exists(), "Never overwrite amended protocol"
    original = json.loads(ORIGINAL_PROTO.read_text())
    failed_dir = OUT / "fits/baseline"
    failed_exec = json.loads((failed_dir / "execution.json").read_text())
    assert failed_exec["returncode"] == -6 and failed_exec["protocol_sha256"] == sha(ORIGINAL_PROTO)
    assert "check_file_docana" in (failed_dir / "fit.log").read_text().lower()
    vpec = REL / "vpec/vpec_baseline_raisin.list"
    jobs = []
    hashes = dict(original["inputs_sha256"])
    hashes[str(Path(__file__).relative_to(ROOT))] = sha(Path(__file__))
    for path in (ORIGINAL_PROTO, failed_dir / "fit.nml", failed_dir / "fit.log", failed_dir / "execution.json"):
        hashes[str(path.relative_to(ROOT))] = sha(path)
    for job in original["jobs"]:
        old = ROOT / job["nml"]
        work = OUT / "fits-shortvpec" / job["condition"]
        work.mkdir(parents=True, exist_ok=False)
        new_vpec = work / "vpec.list"
        shutil.copyfile(vpec, new_vpec)
        assert new_vpec.read_bytes() == vpec.read_bytes()
        nml = work / "fit.nml"
        nml.write_text(set_value(old.read_text(), "HEADER_OVERRIDE_FILE", "'vpec.list'"))
        jobs.append({**job, "nml": str(nml.relative_to(ROOT)), "nml_sha256": sha(nml),
                     "local_vpec_sha256": sha(new_vpec)})
    amended = dict(original)
    amended["state"] = "frozen technical amendment before any successful NIR fit"
    amended["amendment"] = {
        "original_protocol_sha256": sha(ORIGINAL_PROTO),
        "failed_baseline_returncode": -6,
        "failed_log_sha256": sha(failed_dir / "fit.log"),
        "failure": "v11_03c aborted with stack smash in check_file_docana while reading absolute HEADER_OVERRIDE_FILE; no FITRES output. Exact underlying buffer fault not independently localized.",
        "technical_change": "Only HEADER_OVERRIDE_FILE path changed to work-local 'vpec.list'; each file is byte-identical to released vpec_baseline_raisin.list. Original input, source, baseline NML and failure remain preserved.",
    }
    amended["jobs"] = jobs
    amended["inputs_sha256"] = hashes
    dump(PROTO, amended)
    (OUT / "executed-amendment-source.py").write_bytes(Path(__file__).read_bytes())
    print(json.dumps({"amended_protocol_sha256": sha(PROTO), "jobs": len(jobs)}))


def run(condition: str) -> None:
    p = json.loads(PROTO.read_text())
    assert sha(Path(__file__)) == p["inputs_sha256"][str(Path(__file__).relative_to(ROOT))]
    for name, digest in p["inputs_sha256"].items():
        assert sha(ROOT / name) == digest, name
    assert sha(OUT / "header-ledger.csv") == p["header_ledger_sha256"]
    with (OUT / "header-ledger.csv").open(newline="") as f:
        for row in csv.DictReader(f):
            assert sha(PHOT / f"{row['CID']}.snana.dat") == row["phot_source_sha256"]
            assert sha(OUT / "data" / row["condition"] / "DES_RAISIN" / f"{row['CID']}.snana.dat") == row["phot_input_sha256"]
    job = next(x for x in p["jobs"] if x["condition"] == condition)
    work = (ROOT / job["nml"]).parent
    assert sha(ROOT / job["nml"]) == job["nml_sha256"]
    assert sha(work / "vpec.list") == job["local_vpec_sha256"]
    assert not (work / "fit.log").exists()
    env = os.environ.copy()
    sysroot = ROOT / "phase2/official/build/sysroot/usr"
    env.update(SNANA_DIR=str(BUILD), SNDATA_ROOT=str(ROOT / "phase2/official/inputs/SNDATA_ROOT"),
               LD_LIBRARY_PATH=str(sysroot / "lib"), OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1")
    start = time.monotonic()
    with (work / "fit.log").open("x") as log:
        proc = subprocess.run([str(BIN), str(ROOT / job["nml"])], cwd=work, env=env,
                              stdout=log, stderr=subprocess.STDOUT, timeout=900)
    dump(work / "execution.json", {"condition": condition, "returncode": proc.returncode,
                                   "seconds": time.monotonic() - start,
                                   "protocol_sha256": sha(PROTO), "log_sha256": sha(work / "fit.log"),
                                   "binary_sha256": sha(BIN)})
    print(json.dumps({"condition": condition, "returncode": proc.returncode,
                      "seconds": time.monotonic() - start}))


def gate() -> None:
    p = json.loads(PROTO.read_text())
    cohort = p["cohort"]
    baseline = OUT / "fits-shortvpec/baseline"
    copy = OUT / "fits-shortvpec/baseline_copy"
    author = fitres(AUTHOR / "output/fit_nir/DES_RAISIN.FITRES.TEXT")
    rows = fitres(baseline / "fit.FITRES.TEXT")
    copyrows = fitres(copy / "fit.FITRES.TEXT")
    accepted = accepted_counts(baseline / "fit.LCPLOT.TEXT")
    failures = []
    for condition, work in (("baseline", baseline), ("baseline_copy", copy)):
        execution = json.loads((work / "execution.json").read_text())
        if execution["returncode"] != 0 or execution["protocol_sha256"] != sha(PROTO):
            failures.append(condition + " execution")
    if set(rows) != set(cohort) or set(copyrows) != set(cohort):
        failures.append("ten CID membership")
    for cid in cohort:
        if cid not in rows or cid not in copyrows or cid not in author:
            failures.append(cid + " missing row")
            continue
        if rows[cid] != copyrows[cid]:
            failures.append(cid + " copy FITRES mismatch")
        if rows[cid]["ERRFLAG_FIT"] != "0":
            failures.append(cid + " ERRFLAG")
        if not (float(rows[cid]["PKMJDERR"]) == 0 and float(rows[cid]["STRETCH"]) == 1
                and float(rows[cid]["AV"]) == 0 and float(rows[cid]["RV"]) == 1.518):
            failures.append(cid + " fixed coordinates")
        if abs(float(rows[cid]["DLMAG"]) - float(author[cid]["DLMAG"])) > .001:
            failures.append(cid + " author DLMAG difference >.001")
        if float(rows[cid]["NDOF"]) != float(author[cid]["NDOF"]):
            failures.append(cid + " author NDOF")
        if float(rows[cid]["NDOF"]) != accepted.get(cid, 0) - 1:
            failures.append(cid + " accepted epoch NDOF closure")
    if (baseline / "fit.LCPLOT.TEXT").read_bytes() != (copy / "fit.LCPLOT.TEXT").read_bytes():
        failures.append("copy LCPLOT not byte-identical")
    result = {"pass": not failures, "failures": failures, "protocol_sha256": sha(PROTO),
              "baseline_cids": len(rows), "copy_cids": len(copyrows),
              "max_abs_author_DLMAG_difference": max(abs(float(rows[c]["DLMAG"]) - float(author[c]["DLMAG"])) for c in cohort if c in rows and c in author),
              "NDOF_matches": sum(float(rows[c]["NDOF"]) == float(author[c]["NDOF"]) for c in cohort if c in rows and c in author),
              "accepted_epoch_counts": accepted,
              "Rcopy_LCPLOT_byte_identical": (baseline / "fit.LCPLOT.TEXT").read_bytes() == (copy / "fit.LCPLOT.TEXT").read_bytes(),
              "baseline_fitres_sha256": sha(baseline / "fit.FITRES.TEXT"),
              "copy_fitres_sha256": sha(copy / "fit.FITRES.TEXT")}
    dump(OUT / "baseline-gate.json", result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["prepare", "amend-shortvpec", "run", "gate"])
    parser.add_argument("--condition", choices=["baseline", "baseline_copy", *CONDITIONS])
    args = parser.parse_args()
    if args.action == "prepare":
        prepare()
    elif args.action == "amend-shortvpec":
        amend_short_vpec()
    elif args.action == "run":
        assert args.condition
        run(args.condition)
    else:
        gate()
