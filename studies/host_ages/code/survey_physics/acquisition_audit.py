#!/usr/bin/env python3
"""Verify actual Git trees/LFS identities and the official alternate SNANA bundle."""
from concurrent.futures import ThreadPoolExecutor
import datetime
import hashlib
import json
import urllib.error
import urllib.request
from common import ROOT, WORK, RESULTS, sha


def get(url, data=None, headers=None):
    req = urllib.request.Request(url, data=data, headers=headers or {})
    try:
        with urllib.request.urlopen(req, timeout=60) as response:
            return response.status, response.read()
    except urllib.error.HTTPError as exc:
        return exc.code, exc.read()


def main():
    prior = ROOT / "studies/host_ages/results/age_recovery/survey-availability.json"
    objects = json.loads(prior.read_text())["objects"]
    tree = json.loads((WORK / "des-science__DES-SN5YR-tree.json").read_text())
    entries = {item["path"]: item for item in tree["tree"]}

    def pointer(obj):
        url = (
            "https://raw.githubusercontent.com/des-science/DES-SN5YR/"
            + tree["sha"]
            + "/"
            + obj["path"]
        )
        code, data = get(url)
        text = data.decode()
        blob = hashlib.sha1(
            b"blob " + str(len(data)).encode() + b"\0" + data
        ).hexdigest()
        valid = (
            code == 200
            and blob == entries[obj["path"]]["sha"]
            and ("oid sha256:" + obj["sha256"]) in text
            and ("size " + str(obj["bytes"])) in text
        )
        return {
            "path": obj["path"],
            "oid_sha256": obj["sha256"],
            "bytes": obj["bytes"],
            "pointer_http_status": code,
            "git_blob_sha1": blob,
            "tree_pointer_oid_size_verified": valid,
        }

    with ThreadPoolExecutor(max_workers=4) as pool:
        verified = list(pool.map(pointer, objects))
    payload = {
        "operation": "download",
        "transfers": ["basic"],
        "objects": [{"oid": o["sha256"], "size": o["bytes"]} for o in objects],
    }
    request = json.dumps(payload).encode()
    status, body = get(
        "https://github.com/des-science/DES-SN5YR.git/info/lfs/objects/batch",
        request,
        {
            "Content-Type": "application/vnd.git-lfs+json",
            "Accept": "application/vnd.git-lfs+json",
        },
    )
    (WORK / "lfs-batch-response.json").write_bytes(body)
    zenodo = json.loads((WORK / "sndataroot-zenodo.json").read_text())
    archive = WORK / "SNDATA_ROOT_2025-11-12.tar.gz"
    md5 = hashlib.md5()
    with archive.open("rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            md5.update(block)
    assert "md5:" + md5.hexdigest() == zenodo["files"][0]["checksum"]
    releases = json.loads((WORK / "des-science__DES-SN5YR-releases.json").read_text())
    result = {
        "checked_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "code_sha256": sha(__file__),
        "prior_manifest_sha256": sha(prior),
        "git_tree_commit": tree["sha"],
        "git_tree_truncated": tree.get("truncated"),
        "objects_checked": len(verified),
        "pointers_verified": sum(x["tree_pointer_oid_size_verified"] for x in verified),
        "lfs_batch_http_status": status,
        "lfs_batch_response": json.loads(body),
        "lfs_batch_request_sha256": hashlib.sha256(request).hexdigest(),
        "releases": [
            {"tag": r["tag_name"], "asset_count": len(r["assets"])} for r in releases
        ],
        "official_alternative": {
            "url": "https://zenodo.org/records/17591282",
            "file": archive.name,
            "bytes": archive.stat().st_size,
            "md5_verified": md5.hexdigest(),
            "sha256": sha(archive),
            "recovered": [
                "Dovekie SALT3 surfaces and calibration",
                "DES cadence and observed/mocked host libraries",
                "Pipeline and host redshift efficiency",
                "P23/W22 population PDFs",
                "Classifier weights",
                "Nominal simulation, fitting and BBC inputs",
            ],
            "not_recovered": "The original 250 Dovekie mock byte streams. A fresh simulation is not those realizations.",
        },
        "other_official_archive": {
            "url": "https://zenodo.org/records/12720778",
            "scope": "July2024 original DES-SN5YR release precedes Dovekie; not an interchangeable source of revised mock objects.",
        },
        "objects": verified,
    }
    assert all(x["tree_pointer_oid_size_verified"] for x in verified)
    (RESULTS / "acquisition.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "objects"}, indent=2))


if __name__ == "__main__":
    main()
