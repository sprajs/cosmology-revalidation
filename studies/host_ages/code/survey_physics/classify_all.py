#!/usr/bin/env python3
"""Run the recovered classifier environment over a completed frozen campaign."""
import argparse
import json
import os
import subprocess
from common import HERE, WORK, RESULTS, ARCHIVE, sha


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--campaign",
        default="campaign.json",
        choices=["campaign.json", "positive-dust-campaign.json"],
    )
    parser.add_argument(
        "--python",
        default=str(
            ARCHIVE
            / "phase2/classification/python/cpython-3.10.21-linux-x86_64-gnu/bin/python3.10"
        ),
    )
    parser.add_argument(
        "--pythonpath",
        default=str(ARCHIVE / "phase2/classification/env/lib/python3.10/site-packages")
        + ":"
        + str(ARCHIVE / "phase2/classification/sources/SuperNNova"),
    )
    args = parser.parse_args()
    env = os.environ.copy()
    env.update(
        PYTHONPATH=args.pythonpath, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1"
    )
    for job in json.loads((RESULTS / args.campaign).read_text())["jobs"]:
        name = job["name"]
        folder = WORK / "fits" / name
        if (folder / "classifier.json").exists():
            previous = json.loads((folder / "classifier.json").read_text())
            assert previous["native_fit_sha256"] == sha(
                folder / "fit.FITRES.TEXT"
            ), "Changed native fit requires classifier rerun."
            continue
        with (folder / "classify.log").open("w") as out:
            subprocess.run(
                [args.python, str(HERE / "classify.py"), "--name", name],
                env=env,
                stdout=out,
                stderr=subprocess.STDOUT,
                check=True,
            )


if __name__ == "__main__":
    main()
