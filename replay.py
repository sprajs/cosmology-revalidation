#!/usr/bin/env python3
"""Replay each curated workflow once, with explicit failure retention."""

import argparse
import subprocess
import sys
from research import WORKFLOWS, ROOT


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--output", required=True, help="Directory below results for this set of workflow outputs"
    )
    p.add_argument(
        "--skip",
        nargs="*",
        choices=WORKFLOWS,
        default=[],
        help="Explicitly omit workflows; omissions are printed",
    )
    a = p.parse_args()
    for name in WORKFLOWS:
        if name in a.skip:
            print("SKIPPED " + name, flush=True)
            continue
        print("RUN " + name, flush=True)
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "research.py"),
                "run",
                name,
                "--name",
                a.output + "/" + name,
            ],
            check=True,
            cwd=ROOT,
        )


if __name__ == "__main__":
    main()
