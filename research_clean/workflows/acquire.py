"""Restore public input files, refusing unverified substitutions."""

import json
import os
import urllib.request
from lib.paths import ROOT
from lib.records import sha256


def restore(items, group=None, path=None, dry_run=False):
    rows = [
        x
        for x in items
        if (group is None or x["group"] == group)
        and (path is None or x["path"] == path)
    ]
    if not rows:
        raise ValueError("No input matched the request")
    failed = []
    for item in rows:
        dest = ROOT / item["path"]
        if dest.is_file():
            if sha256(dest) != item["sha256"]:
                failed.append(
                    {
                        "path": item["path"],
                        "error": "Existing file differs; preserve it and resolve explicitly",
                    }
                )
            continue
        if dry_run:
            print(
                json.dumps(
                    {
                        "path": item["path"],
                        "urls": item["urls"],
                        "available_from_public_url": bool(item["urls"]),
                    }
                )
            )
            continue
        if not item["urls"]:
            failed.append(
                {
                    "path": item["path"],
                    "error": "No verified public URL: restore this frozen derived input from the distributed bundle",
                }
            )
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        partial = dest.with_name(dest.name + ".partial")
        errors = []
        for url in item["urls"]:
            try:
                request = urllib.request.Request(
                    url, headers={"User-Agent": "supernova-clean-research/0.1"}
                )
                with (
                    urllib.request.urlopen(request, timeout=90) as response,
                    partial.open("xb") as stream,
                ):
                    count = 0
                    while block := response.read(1024 * 1024):
                        count += len(block)
                        if count > item["bytes"]:
                            raise ValueError("Response exceeds expected byte count")
                        stream.write(block)
                if (
                    partial.stat().st_size != item["bytes"]
                    or sha256(partial) != item["sha256"]
                ):
                    raise ValueError("Downloaded bytes differ from the frozen input")
                os.replace(partial, dest)
                print("Restored " + item["path"])
                break
            except Exception as error:
                errors.append(str(error))
                partial.unlink(missing_ok=True)
        else:
            failed.append({"path": item["path"], "errors": errors})
    if failed:
        print(json.dumps(failed, indent=2))
        raise SystemExit(1)
