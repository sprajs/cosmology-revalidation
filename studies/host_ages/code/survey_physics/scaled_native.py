#!/usr/bin/env python3
"""Independent, paired native training shards; preserve the earlier campaign."""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import datetime
import json
import os
from pathlib import Path
import subprocess
import time
from common import HERE, WORK, RESULTS, DATA, ARCHIVE, sha
from prepare import replace_key
import run_native

SCALE = WORK / "scaled-bbc"
MANIFEST = RESULTS / "scaled-campaign.json"


def prepare(shards, attempts):
    assert shards > 0 and attempts > 0
    for sub in ["inputs", "fits", "simulations", "bbc"]:
        (SCALE / sub).mkdir(parents=True, exist_ok=True)
    link = SCALE / "SNANA"
    if not link.exists():
        link.symlink_to(WORK / "SNANA", target_is_directory=True)
    previous = json.loads(MANIFEST.read_text()) if MANIFEST.exists() else None
    if previous:
        assert previous["attempts_per_shard"] == attempts
        assert shards >= previous["paired_shards"]
    jobs = []
    for i in range(shards):
        for arm, label in [("nominal", "N"), ("age", "A")]:
            name = f"SPB_{label}{i:03d}"
            source = WORK / "inputs" / ("SPP_TRAIN_" + arm.upper() + ".input")
            text = source.read_text()
            changes = {
                "GENVERSION": name,
                "GENPREFIX": name,
                "NGENTOT_LC": attempts,
                "RANSEED": 271000 + i,
                "CIDOFF": 1000000 * (i + 1),
                "PATH_SNDATA_SIM": SCALE / "simulations",
                "GENMODEL_EXTRAP_LATETIME": WORK
                / "inputs/SNIa_Extrap_LateTime_2expon.TEXT",
            }
            for k, v in changes.items():
                text = replace_key(text, k, v)
            path = SCALE / "inputs" / (name + ".input")
            if path.exists():
                assert (
                    path.read_text() == text
                ), "Frozen shard changed; use a new campaign."
            else:
                path.write_text(text)
            jobs.append(
                {
                    "name": name,
                    "arm": arm,
                    "shard": i,
                    "attempts": attempts,
                    "seed": 271000 + i,
                    "cid_offset": 1000000 * (i + 1),
                    "input": str(path),
                    "input_sha256": sha(path),
                    "physics_parent_sha256": sha(source),
                }
            )
    r = {
        "frozen_utc": (
            previous["frozen_utc"]
            if previous
            else datetime.datetime.now(datetime.timezone.utc).isoformat()
        ),
        "paired_shards": shards,
        "attempts_per_shard": attempts,
        "attempts_per_arm": shards * attempts,
        "attempts_total": 2 * shards * attempts,
        "jobs": jobs,
        "design": {
            "population": "Same declared W22 mocked host ages and width distribution, fixedRV3.1 physically positive screens, inherited intrinsiccolour and EBV. Extra greyshift −.030*(SN_age−3) enters photons before selection.",
            "seeds": "Each shard has a distinct seed; nominal/age counterparts share that seed. All training seeds differ from fixed evaluation seed270922.",
            "evaluation": "Retain existing3000-attempt nominal/age physical evaluation arms; do not select evaluation based on outcome.",
            "identity": "CIDOFF is unique acrossshards, shared within nominal/age counterpart; validate(shard,CID,LIBID) truth and independent evaluation.",
            "classifier": "Unchanged reconstructed SNNV19 and pIa>.5; pureIa signal acceptance only.",
            "native_BBC": "opt_biascor4336, released field groups, binning, SNR60 intrinsic scatter gate and interpolation minimum counts unchanged. Native capacity bounds fix only. No CCprior or BEAMS purity claim.",
            "estimands": "Compare nominal evaluation, injected evaluation with frozen nominal nuisance/correction, injected evaluation with nuisance refit but nominal training, and independently age-retrained correction. Distinguish native conditional covariance from training MonteCarlo.",
            "adaptation": "Inspect first independent shards for throughput, then accumulate training until support stabilizes or a specific independent barrier is established; never relax scientific support gates.",
            "not_claimed": [
                "Actual progenitor-age measurement",
                "Full survey likelihood or CC population",
                "SALT surface or classifier retraining",
                "Empirical cosmological correction",
            ],
        },
        "code_sha256": sha(__file__),
        "original_driver_sha256": sha(HERE / "run_native.py"),
        "classifier_wrapper_sha256": sha(HERE / "scaled_classifier.py"),
    }
    if previous:
        r["first_executed_driver_sha256"] = previous.get(
            "first_executed_driver_sha256", previous["code_sha256"]
        )
        r["volume_expansions"] = previous.get("volume_expansions", [])
        if r["code_sha256"] != previous["code_sha256"]:
            r["coordinator_update"] = (
                "Added explicit classifier Python/PYTHONPATH overrides and atomic per-job progress writes to prevent concurrent partial-JSON reads; frozen native inputs, helper and numerical outputs unchanged and cache hashes reverified."
            )
    if previous and shards > previous["paired_shards"]:
        r["volume_expansions"] = previous.get("volume_expansions", []) + [
            {
                "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                "shards": shards,
                "attempts_total": r["attempts_total"],
            }
        ]
    MANIFEST.write_text(json.dumps(r, indent=2) + "\n")
    print(
        json.dumps(
            {k: r[k] for k in ["paired_shards", "attempts_per_arm", "attempts_total"]},
            indent=2,
        )
    )


def one(job, classifier_python=None, classifier_pythonpath=None):
    run_native.WORK = SCALE
    run = run_native.run(job)
    if not run.get("fit_graceful"):
        return run
    folder = SCALE / "fits" / job["name"]
    dest = folder / "classifier.json"
    if dest.exists():
        assert json.loads(dest.read_text())["native_fit_sha256"] == sha(
            folder / "fit.FITRES.TEXT"
        )
    else:
        env = os.environ.copy()
        env.update(
            OPENBLAS_NUM_THREADS="1",
            OMP_NUM_THREADS="1",
            PYTHONPATH=classifier_pythonpath
            or str(ARCHIVE / "phase2/classification/env/lib/python3.10/site-packages")
            + ":"
            + str(ARCHIVE / "phase2/classification/sources/SuperNNova"),
        )
        python = classifier_python or (
            ARCHIVE
            / "phase2/classification/python/cpython-3.10.21-linux-x86_64-gnu/bin/python3.10"
        )
        start = time.monotonic()
        with (folder / "classify.log").open("w") as out:
            subprocess.run(
                [
                    str(python),
                    str(HERE / "scaled_classifier.py"),
                    "--name",
                    job["name"],
                ],
                env=env,
                stdout=out,
                stderr=subprocess.STDOUT,
                check=True,
            )
        run["classifier_process_seconds"] = time.monotonic() - start
    record = {
        "name": job["name"],
        "arm": job["arm"],
        "shard": job["shard"],
        "attempts": job["attempts"],
        "generation_seconds": run["generation_seconds"],
        "fit_seconds": run["fit_seconds"],
        "fitres_rows": run["fit.FITRES.TEXT"]["rows"],
        "native_run_sha256": sha(folder / "run.json"),
        "classifier": json.loads(dest.read_text()),
    }
    temporary = folder / "scaled-run.json.tmp"
    temporary.write_text(json.dumps(record, indent=2) + "\n")
    temporary.replace(folder / "scaled-run.json")
    print(
        json.dumps(
            {
                "completed": job["name"],
                "attempts": job["attempts"],
                "fitres": record["fitres_rows"],
                "seconds": run["generation_seconds"] + run["fit_seconds"],
            }
        ),
        flush=True,
    )
    return record


def execute(workers, classifier_python=None, classifier_pythonpath=None):
    assert 1 <= workers <= 12
    jobs = json.loads(MANIFEST.read_text())["jobs"]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [
            pool.submit(one, j, classifier_python, classifier_pythonpath) for j in jobs
        ]
        for future in as_completed(futures):
            future.result()
            records = [
                json.loads(p.read_text())
                for p in sorted((SCALE / "fits").glob("SPB_*/scaled-run.json"))
            ]
            (RESULTS / "scaled-native-execution.json").write_text(
                json.dumps(
                    {
                        "recorded_utc": datetime.datetime.now(
                            datetime.timezone.utc
                        ).isoformat(),
                        "workers": workers,
                        "completed_shards": len(records),
                        "records": records,
                    },
                    indent=2,
                )
                + "\n"
            )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("stage", choices=["prepare", "run"])
    p.add_argument("--shards", type=int, default=8)
    p.add_argument("--attempts", type=int, default=15000)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument(
        "--classifier-python",
        help="Python interpreter in the separately restored pinned classifier environment.",
    )
    p.add_argument(
        "--classifier-pythonpath",
        help="Pinned Torch site-packages and SuperNNova checkout, separated by a colon.",
    )
    a = p.parse_args()
    if a.stage == "prepare":
        prepare(a.shards, a.attempts)
    else:
        execute(a.workers, a.classifier_python, a.classifier_pythonpath)


if __name__ == "__main__":
    main()
