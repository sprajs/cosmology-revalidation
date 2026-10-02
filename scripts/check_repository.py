#!/usr/bin/env python3
"""Check packet contracts, local documentation links and the public storage budget."""
from pathlib import Path
import re
import subprocess

from packet import ROOT, validate_all

LOCAL = {"data", "downloads", "results", "runs", "simulation", "simulations", "notebooks", ".work", ".venv"}
GENERATED = (".ipynb", ".npy", ".npz", ".fits", ".pdf", ".png", ".svg", ".tar.gz", ".zip")
PUBLIC_SOURCE_BYTES = 2 * 1024 * 1024


def check():
    count = validate_all()
    files = subprocess.check_output(["git", "-C", str(ROOT), "ls-files", "-z"]).decode().split("\0")[:-1]
    total = 0
    for name in files:
        path = ROOT / name
        if not path.is_file():
            raise ValueError(f"Tracked path is missing; stage deletions before checking: {name}")
        if Path(name).parts[0] in LOCAL or name.endswith(GENERATED):
            raise ValueError(f"Generated/local artifact is tracked: {name}")
        total += path.stat().st_size
        if path.suffix == ".md":
            for target in re.findall(r"\]\(([^)]+)\)", path.read_text()):
                target = target.split("#", 1)[0]
                if not target or "://" in target or target.startswith("mailto:"):
                    continue
                if not (path.parent / target).exists():
                    raise ValueError(f"Broken local link in {name}: {target}")
    if total > PUBLIC_SOURCE_BYTES:
        raise ValueError("Public tree exceeds 2 MiB; archive bulk evidence and review the storage design")
    packets = {}
    for name in files:
        parts = Path(name).parts
        if len(parts) >= 3 and parts[0] == "experiments":
            packets[parts[1]] = packets.get(parts[1], 0) + 1
    for name, size in packets.items():
        if size > 8:
            raise ValueError(f"Packet exceeds eight public source files: {name}")
    print(f"Checked {count} executable/blocked packets, {len(files)} public files, {total:,} bytes.")


if __name__ == "__main__":
    check()
