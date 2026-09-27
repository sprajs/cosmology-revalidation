#!/usr/bin/env python3
"""Acquire public input bytes and verify checksums; never install system packages."""
import argparse
import gzip
import hashlib
import json
import shutil
import tarfile
import urllib.request
from common import WORK, DATA, ARCHIVE, sha


def fetch(url, dest):
    if not dest.exists():
        temporary = dest.with_name(dest.name + ".part")
        with (
            urllib.request.urlopen(url, timeout=60) as source,
            temporary.open("wb") as target,
        ):
            shutil.copyfileobj(source, target, 2**20)
        temporary.replace(dest)
    return dest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--extract", action="store_true")
    args = parser.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)
    urls = {
        "release-era-sntools_trigger.c": "https://raw.githubusercontent.com/RickKessler/SNANA/8ccf23a17cd98ac667dbc9e9666970c1e36246af/src/sntools_trigger.c",
        "sndataroot-zenodo.json": "https://zenodo.org/api/records/17591282",
        "des-science__DES-SN5YR-tree.json": "https://api.github.com/repos/des-science/DES-SN5YR/git/trees/c9a4fcafc4cbd19bd750dee47fc76194a45c181f?recursive=1",
        "des-science__DES-SN5YR-releases.json": "https://api.github.com/repos/des-science/DES-SN5YR/releases",
    }
    for name, url in urls.items():
        fetch(url, WORK / name)
    z = json.loads((WORK / "sndataroot-zenodo.json").read_text())
    entry = next(x for x in z["files"] if x["key"] == "SNDATA_ROOT_2025-11-12.tar.gz")
    archive = fetch(
        "https://zenodo.org/api/records/17591282/files/SNDATA_ROOT_2025-11-12.tar.gz/content",
        WORK / entry["key"],
    )
    h = hashlib.md5()
    with archive.open("rb") as stream:
        for chunk in iter(lambda: stream.read(2**20), b""):
            h.update(chunk)
    assert (
        archive.stat().st_size == entry["size"]
        and "md5:" + h.hexdigest() == entry["checksum"]
    )
    if args.extract and not (DATA / "models/SALT3/SALT3.DOVEKIE").exists():
        DATA.mkdir(exist_ok=True)
        with tarfile.open(archive) as source:
            source.extractall(DATA, filter="data")
    inputs = WORK / "inputs"
    inputs.mkdir(exist_ok=True)
    p = fetch(
        "https://zenodo.org/api/records/6672739/files/SNIa_Extrap_LateTime_2expon.TEXT.gz/content",
        inputs / "SNIa_Extrap_LateTime_2expon.TEXT.gz",
    )
    assert hashlib.md5(p.read_bytes()).hexdigest() == "da08780daf6386b1758352e930bbc775"
    unpacked = inputs / "SNIa_Extrap_LateTime_2expon.TEXT"
    unpacked.write_bytes(gzip.decompress(p.read_bytes()))
    previous = ARCHIVE / "phase2/literature/sources/SNIa_Extrap_LateTime_2expon.TEXT"
    if previous.exists():
        assert sha(unpacked) == sha(previous)
    (DATA / "SIM").mkdir(exist_ok=True)
    (WORK / "simulations").mkdir(exist_ok=True)
    (DATA / "SIM/PATH_SNDATA_SIM.LIST").write_text(str(WORK / "simulations") + "\n")
    print(
        json.dumps(
            {
                "bundle_sha256": sha(archive),
                "late_time_model_sha256": sha(unpacked),
                "native_build": "See README; pinned public source and output-only patch, no downloaded code committed.",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
