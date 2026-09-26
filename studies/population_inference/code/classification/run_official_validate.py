#!/usr/bin/env python3
"""Run pinned SuperNNova validate_rnn on isolated DES test HDF5."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reconstruct import MODEL, OUT, RAW, SOURCE


def main():
    task = json.loads((OUT / "test_database_command.json").read_text())
    argv = [sys.executable, str(SOURCE / "run.py"), "--validate_rnn",
            "--model_files", str(MODEL), "--dump_dir", str(OUT / "test_database"),
            "--raw_dir", str(RAW), "--sntypes", json.dumps(task["sntypes"]),
            "--list_filters", "g", "i", "r", "z", "--redshift", "zspe",
            "--norm", "cosmo", "--data_testing"]
    (OUT / "official_validate_command.json").write_text(json.dumps(argv, indent=2) + "\n")
    with (OUT / "official_validate.log").open("w") as log:
        subprocess.run(argv, cwd=SOURCE, stdout=log, stderr=subprocess.STDOUT, check=True)
    print("Official validate_rnn completed")


if __name__ == "__main__":
    main()
