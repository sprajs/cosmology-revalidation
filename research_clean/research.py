#!/usr/bin/env python3
"""One entry point for acquisition, verification and reproducible calculations."""

import argparse
import datetime
import importlib
import importlib.metadata
import json
import os
import platform
import sys
import time
from pathlib import Path

# Set before importing NumPy/JAX; one CPU thread avoids uncontrolled fit fan-out.
for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"
os.environ.setdefault("JAX_PLATFORMS", "cpu")
os.environ.setdefault(
    "XLA_FLAGS", "--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1"
)
os.environ.setdefault("MPLBACKEND", "Agg")

from lib.paths import ROOT

WORKFLOWS = {
    "cosmology": (
        "distances",
        "Released-distance likelihood and cosmological fits",
        ["distances", "bao"],
    ),
    "bao-shape": ("bao", "Flat nonaccelerating shape-cone projection", ["bao"]),
    "ages": (
        "ages",
        "Published host-age extraction and conditional regression",
        ["distances", "ages"],
    ),
    "populations": (
        "populations",
        "Star formation and progenitor delay convolution",
        [],
    ),
    "dust": ("dust", "Absorption geometry and latent-colour degeneracy", []),
    "des-flux": ("des_flux", "Independent SALT3 calibrated-flux fitting", ["des_flux"]),
    "des-predictors": (
        "des_predictors",
        "Convergence-gated selected-sample prediction",
        ["des_predictors"],
    ),
    "calibration": (
        "calibration",
        "Shared calibration-mode inference and distance information",
        ["calibration"],
    ),
    "raisin": (
        "raisin",
        "Paired released distances, covariance and signed-flux inventory",
        ["raisin"],
    ),
    "timing": (
        "timing",
        "Coherent archived simulation timing and selection",
        ["timing"],
    ),
    "csp-lineage": (
        "csp_lineage",
        "Raw photometry to released physical-filter lineage",
        ["raisin", "csp"],
    ),
    "csp-passbands": (
        "csp_passbands",
        "Photon integrals and conditional passband contrasts",
        ["raisin", "csp"],
    ),
    "hst-repeat": (
        "hst",
        "Signed HST repeat photometry with frozen operators",
        ["hst"],
    ),
    "hst-geometry": (
        "hst_geometry",
        "Rebuild HST aperture geometry and source masks",
        ["hst"],
    ),
    "hst-flat": ("hst_flat", "Propagate quoted HST flat-reference variance", ["hst"]),
    "hst-dark-raw": (
        "hst_dark_raw",
        "Raw dark read moments and fixed signed aperture slopes",
        ["hst_dark_geometry", "hst_dark_raw"],
    ),
    "hst-dark-calibrated": (
        "hst_dark_calibrated",
        "Native dark slope/variance and influence accounting",
        ["hst_dark_geometry", "hst_dark_calibrated"],
    ),
    "signed-baseline": (
        "signed_baseline",
        "Signed author photometry and pre-explosion noise",
        ["signed"],
    ),
    "sign-selection": (
        "sign_selection",
        "Gaussian sign likelihood and synthetic amplitude recovery",
        [],
    ),
}


def inputs():
    return json.loads((ROOT / "provenance/inputs.json").read_text())["files"]


def verify(groups=None):
    from lib.records import sha256

    selected = [x for x in inputs() if groups is None or x["group"] in groups]
    failures = []
    for item in selected:
        path = ROOT / item["path"]
        if not path.is_file():
            failures.append({"path": item["path"], "error": "missing"})
        elif path.stat().st_size != item["bytes"] or sha256(path) != item["sha256"]:
            failures.append({"path": item["path"], "error": "size or SHA-256 mismatch"})
    return {"checked": len(selected), "pass": not failures, "failures": failures}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="action", required=True)
    sub.add_parser("list", help="List the curated calculations")
    sub.add_parser(
        "verify", help="Verify all bundled input bytes without changing them"
    )
    f = sub.add_parser(
        "fetch", help="Restore missing public inputs using recorded URLs and hashes"
    )
    f.add_argument("--group", choices=sorted({x["group"] for x in inputs()}))
    f.add_argument("--path", help="Restore a single manifest path, e.g. data/bao/...")
    f.add_argument("--dry-run", action="store_true")
    r = sub.add_parser("run", help="Run one workflow in a new output directory")
    r.add_argument("workflow", choices=WORKFLOWS)
    r.add_argument(
        "--name",
        help="Unique result directory name; existing results are never overwritten",
    )
    r.add_argument(
        "--config", type=Path, help="JSON settings, merged with documented defaults"
    )
    a = p.parse_args()
    if a.action == "list":
        for name, (_, description, _) in WORKFLOWS.items():
            print(f"{name:18s} {description}")
        return
    if a.action == "verify":
        report = verify()
        print(json.dumps(report, indent=2))
        if not report["pass"]:
            raise SystemExit(1)
        return
    if a.action == "fetch":
        from workflows.acquire import restore

        restore(inputs(), a.group, a.path, a.dry_run)
        return
    module_name, _, groups = WORKFLOWS[a.workflow]
    report = verify(groups)
    if not report["pass"]:
        print(json.dumps(report, indent=2), file=sys.stderr)
        raise SystemExit("Input verification failed; no analysis was started")
    from lib.records import sha256, write_json

    code = [
        ROOT / "research.py",
        *sorted((ROOT / "lib").glob("*.py")),
        *sorted((ROOT / "workflows").glob("*.py")),
    ]
    source_hashes = {str(x.relative_to(ROOT)): sha256(x) for x in code}
    module = importlib.import_module("workflows." + module_name)
    config = dict(module.DEFAULTS)
    if a.config:
        supplied = json.loads(a.config.read_text())
        unknown = set(supplied) - set(config)
        if unknown:
            raise ValueError(f"Unknown settings for {a.workflow}: {sorted(unknown)}")
        config.update(supplied)
    name = a.name or a.workflow + "-" + datetime.datetime.now().strftime(
        "%Y%m%d-%H%M%S"
    )
    if name in {".", ".."} or Path(name).name != name:
        raise ValueError("Result name must be one directory name")
    out = ROOT / "results" / name
    out.mkdir(parents=True, exist_ok=False)
    packages = {}
    for package in (
        "numpy",
        "scipy",
        "pandas",
        "astropy",
        "emcee",
        "sncosmo",
        "extinction",
        "iminuit",
        "jax",
        "jaxlib",
        "numpyro",
    ):
        try:
            packages[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            pass
    record = {
        "workflow": a.workflow,
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "config": config,
        "python": platform.python_version(),
        "packages": packages,
        "platform": platform.platform(),
        "inputs_sha256": {
            x["path"]: x["sha256"] for x in inputs() if x["group"] in groups
        },
        "code_sha256": source_hashes,
        "lock_sha256": sha256(ROOT / "uv.lock"),
        "status": "running",
    }
    write_json(out / "run.json", record)
    start = time.monotonic()
    try:
        result = module.run(out, config)
        write_json(out / "summary.json", result)
        record["status"] = "completed"
    except BaseException as error:
        record["status"] = "failed"
        record["error"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        record["elapsed_seconds"] = time.monotonic() - start
        record["outputs_sha256"] = {
            str(x.relative_to(out)): sha256(x)
            for x in sorted(out.rglob("*"))
            if x.is_file() and x.name != "run.json"
        }
        record["source_changed_during_run"] = [
            str(x.relative_to(ROOT))
            for x in code
            if sha256(x) != record["code_sha256"][str(x.relative_to(ROOT))]
        ]
        if record["source_changed_during_run"] and record["status"] == "completed":
            record["status"] = "completed_with_source_change"
        write_json(out / "run.json", record)
    print(
        json.dumps(
            {
                "workflow": a.workflow,
                "output": str(out),
                "elapsed_seconds": round(record["elapsed_seconds"], 2),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
