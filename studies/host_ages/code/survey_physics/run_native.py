#!/usr/bin/env python3
"""Execute frozen generated-attempt simulations and native SALT fits, at most four jobs."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime
import json
import re
import subprocess
import time
from pathlib import Path
from common import HERE, WORK, RESULTS, DATA, sha, native_env


def run(job):
    name = job["name"]
    folder = WORK / "fits" / name
    folder.mkdir(parents=True, exist_ok=True)
    assert sha(job["input"]) == job["input_sha256"]
    manifest = folder / "run.json"
    if manifest.exists():
        previous = json.loads(manifest.read_text())
        if previous.get("fit_graceful"):
            recorded = previous.get("hashes", {}).get(job["input"])
            assert (
                recorded == job["input_sha256"]
            ), "Changed frozen input requires a new version name, not cached outputs."
            return previous
    env = native_env()
    simexe = WORK / "SNANA/bin/snlc_sim.exe"
    fitexe = WORK / "SNANA/bin/snlc_fit.exe"
    start = time.monotonic()
    with (folder / "generation.log").open("w") as output:
        p = subprocess.run(
            [str(simexe), job["input"]],
            env=env,
            cwd=folder,
            stdout=output,
            stderr=subprocess.STDOUT,
        )
    simseconds = time.monotonic() - start
    generated = (folder / "generation.log").read_text()
    if p.returncode or "DONE with snlc_sim." not in generated:
        result = {
            "name": name,
            "generation_returncode": p.returncode,
            "generation_seconds": simseconds,
            "fit_graceful": False,
        }
        manifest.write_text(json.dumps(result, indent=2) + "\n")
        return result
    nmlsrc = DATA / "sample_input_files/DES-SN5YR/base_files/lcfit/lcfit_desSMP_5yr.nml"
    text = nmlsrc.read_text()
    changes = {
        "PRIVATE_DATA_PATH": str(WORK / "simulations"),
        "VERSION_PHOTOMETRY": name,
        "KCOR_FILE": "$SNDATA_ROOT/kcor/Dovekie/calib_DES-SN5YR_DES.fits.gz",
        "FITMODEL_NAME": "$SNDATA_ROOT/models/SALT3/SALT3.DOVEKIE",
        # Upstream detects dots anywhere in prefix; basename avoids .work collision.
        "TEXTFILE_PREFIX": "fit",
        "SNTABLE_LIST": "FITRES(text:host) SNANA(text:host)",
    }
    for key, value in changes.items():
        pattern = r"(?m)^\s*" + key + r"\s*=.*$"
        assert len(re.findall(pattern, text)) == 1, key
        text = re.sub(pattern, "    " + key + " = '" + value + "'", text)
    # Published Dovekie override, modern exact F99 (default false, also explicit).
    text = text.replace("&SNLCINP", "&SNLCINP\n    RESTORE_DES5YR = F")
    nml = folder / "fit.nml"
    nml.write_text(text)
    start = time.monotonic()
    with (folder / "fit.log").open("w") as output:
        p = subprocess.run(
            [str(fitexe), str(nml)],
            env=env,
            cwd=folder,
            stdout=output,
            stderr=subprocess.STDOUT,
        )
    fitseconds = time.monotonic() - start
    log = (folder / "fit.log").read_text()
    result = {
        "name": name,
        "finished_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "code_sha256": sha(__file__),
        "generation_seconds": simseconds,
        "fit_seconds": fitseconds,
        "generation_returncode": 0,
        "fit_returncode": p.returncode,
        "fit_graceful": p.returncode == 0 and "ENDING PROGRAM GRACEFULLY" in log,
        "minuit_warning_occurrences": len(re.findall("MINUIT WARNING", log)),
        "generated_attempts_requested": job["attempts"],
        "hashes": {
            str(path): sha(path)
            for path in [
                Path(job["input"]),
                nml,
                nmlsrc,
                simexe,
                fitexe,
                folder / "generation.log",
                folder / "fit.log",
            ]
        },
    }
    for filename in ["fit.FITRES.TEXT", "fit.SNANA.TEXT"]:
        f = folder / filename
        if f.exists():
            result[filename] = {
                "sha256": sha(f),
                "rows": sum(line.startswith("SN:") for line in f.open()),
            }
    manifest.write_text(json.dumps(result, indent=2) + "\n")
    print(name, result["fit_graceful"], round(simseconds + fitseconds, 1), flush=True)
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--name")
    p.add_argument(
        "--campaign",
        choices=["campaign.json", "positive-dust-campaign.json"],
        default="campaign.json",
    )
    args = p.parse_args()
    assert 1 <= args.workers <= 4
    jobs = json.loads((RESULTS / args.campaign).read_text())["jobs"]
    if args.name:
        jobs = [j for j in jobs if j["name"] == args.name]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        outcomes = list(pool.map(run, jobs))
    # Rebuild from completed run manifests so a filtered invocation cannot erase other jobs.
    paths = sorted((WORK / "fits").glob("SP[HP]_*/run.json"))
    outcomes = [json.loads(path.read_text()) for path in paths]
    (RESULTS / "native-execution.json").write_text(
        json.dumps(outcomes, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
