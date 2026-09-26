#!/usr/bin/env python3
"""Fetch MAST CAOM cone metadata only; never fetch dataURI products."""
import hashlib
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
URL = "https://mast.stsci.edu/api/v0/invoke"
LIMIT = 20_000_000


def query(request, remaining):
    body = urllib.parse.urlencode({"request": json.dumps(request, separators=(",", ":"))}).encode()
    req = urllib.request.Request(URL, data=body, headers={"Content-Type": "application/x-www-form-urlencoded", "User-Agent": "CodexRAISINMetadataAudit/1.0"})
    with urllib.request.urlopen(req, timeout=60) as response:
        raw = response.read(remaining + 1)
    if len(raw) > remaining:
        raise ValueError("metadata acquisition cap exceeded")
    return raw


def main():
    start = time.monotonic()
    req = {"service": "Mast.Caom.Cone", "params": {"ra": 8.9685, "dec": -43.35812, "radius": 30 / 3600}, "format": "json", "pagesize": 2000, "page": 1}
    raw = query(req, LIMIT)
    result = json.loads(raw)
    if result.get("status") != "COMPLETE":
        raise ValueError(("MAST query incomplete", result.get("status"), result.get("msg")))
    HERE.mkdir(parents=True, exist_ok=True)
    (HERE / "caom-cone.json").write_bytes(raw)
    (HERE / "caom-request.json").write_text(json.dumps(req, indent=2) + "\n")
    summary = {"status": result.get("status"), "rows": len(result.get("data", [])), "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest(), "elapsed_seconds": time.monotonic() - start, "paging": result.get("paging")}
    (HERE / "caom-summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
