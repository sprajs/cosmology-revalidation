#!/usr/bin/env python3
"""Build isolated DES SuperNNova prediction database to recover its normalization."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import h5py
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reconstruct import OUT, RAW, SOURCE


def main():
    with fits.open(str(RAW / "DES-SN5YR_DES_HEAD.FITS.gz"), memmap=False) as f:
        codes = sorted(set(map(int, f[1].data["SNTYPE"])))
    types = {str(code): ("Ia" if code in (1, 101) else "II") for code in codes}
    dump = OUT / "test_database"
    dump.mkdir(parents=True, exist_ok=True)
    clump = OUT / "clump_run/DES-SN5YR_DES.SNANA.TEXT"
    cmd = [sys.executable, str(SOURCE / "run.py"), "--data", "--data_testing",
           "--sntypes", json.dumps(types), "--list_filters", "g", "i", "r", "z",
           "--phot_reject", "PHOTFLAG", "--phot_reject_list", "8", "16", "32", "64", "128", "256", "512",
           "--redshift_label", "REDSHIFT_FINAL", "--photo_window_files", str(clump),
           "--norm", "cosmo", "--dump_dir", str(dump), "--raw_dir", str(RAW)]
    (OUT / "test_database_command.json").write_text(json.dumps({"argv": cmd, "sntypes": types}, indent=2) + "\n")
    with (OUT / "test_database.log").open("w") as log:
        subprocess.run(cmd, cwd=SOURCE, stdout=log, stderr=subprocess.STDOUT, check=True)
    h5 = dump / "processed/database.h5"
    with h5py.File(h5, "r") as f:
        norm = {name: {key: float(f[f"normalizations/{name}/{key}"][()])
                       for key in ("min", "mean", "std")}
                for name in ("delta_time", "FLUXCAL_g", "FLUXCALERR_g")}
        norm["features"] = f["features"][:].astype(str).tolist()
        norm["n_objects"] = len(f["SNID"])
    (OUT / "test_database_norm.json").write_text(json.dumps(norm, indent=2) + "\n")
    print(json.dumps({"delta_time": norm["delta_time"], "n_objects": norm["n_objects"]}, indent=2))


if __name__ == "__main__":
    main()
