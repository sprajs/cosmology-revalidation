#!/usr/bin/env python3
"""Restore only this experiment's inherited public inputs, with byte checks.

A local old-layout archive is an optional cache, never modified. Without it,
ZTF restoration downloads the public 1.37 GB release ZIP, extracting only three
explicit CSV members. The existing main acquisition supplies C25 age tables
and Pantheon+ distances (see docs/methods/data.md).
"""

from pathlib import Path
import argparse, json, shutil, urllib.request, zipfile
from acquire import ROOT, WORK, OUT, sha
from match import ARCHIVE


def download(url, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".partial")
    with urllib.request.urlopen(url, timeout=120) as response, tmp.open("wb") as stream:
        shutil.copyfileobj(response, stream)
    tmp.replace(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--destination", type=Path, default=WORK / "prerequisites")
    ap.add_argument("--archive-cache", type=Path, default=ARCHIVE)
    args = ap.parse_args()
    lock = json.loads(Path(__file__).with_name("prerequisite-lock.json").read_text())[
        "files"
    ]
    records = []
    for rel, expected in lock.items():
        target = args.destination / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        record = {"path": rel, "sha256": expected}
        if target.exists():
            record["source"] = "existing verified destination"
        elif (args.archive_cache / rel).exists():
            shutil.copyfile(args.archive_cache / rel, target)
            record["source"] = "optional local frozen input cache"
        elif "gupta2011" in rel:
            url = "https://cdsarc.cds.unistra.fr/ftp/J/ApJ/740/92/" + Path(rel).name
            download(url, target)
            record["source"] = url
        elif "chung2025" in rel:
            source = ROOT / "data/ages" / Path(rel).name
            if not source.exists():
                raise FileNotFoundError(
                    "Run the main ages acquisition first: missing " + str(source)
                )
            shutil.copyfile(source, target)
            record["source"] = str(source.relative_to(ROOT))
        elif "Ginolin" in rel:
            url = "https://raw.githubusercontent.com/mginolin/standax/903d65df4f52c7096ff44e9f2dc4acd0bcbd36a8/notebooks/Ginolin25ab_masterlist.csv"
            download(url, target)
            record["source"] = url
        elif "ztfsniadr2_lite" in rel:
            url = "https://ztfcosmo.in2p3.fr/download/data"
            p = WORK / "ztfsniadr2_lite.zip"
            if not p.exists():
                download(url, p)
            with zipfile.ZipFile(p) as archive:
                member = "ztfsniadr2_lite/tables/" + Path(rel).name
                if member not in archive.namelist():
                    raise KeyError("Explicit expected ZIP member absent: " + member)
                target.write_bytes(archive.read(member))
            record["source"] = url + "#" + member
            record["zip_sha256"] = sha(p)
        else:
            raise ValueError(rel)
        if sha(target) != expected:
            raise ValueError("Input byte identity differs: " + str(target))
        records.append(record)
    (OUT / "prerequisite-restoration.json").write_text(
        json.dumps(
            {
                "code_sha256": sha(Path(__file__)),
                "lock_sha256": sha(Path(__file__).with_name("prerequisite-lock.json")),
                "records": records,
            },
            indent=2,
        )
        + "\n"
    )
    print("--archive", args.destination)


if __name__ == "__main__":
    main()
