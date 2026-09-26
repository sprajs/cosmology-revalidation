#!/usr/bin/env python3
"""Fetch only CAOM product metadata for the frozen RAISIN2 target observations."""
import concurrent.futures
import hashlib
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
URL = "https://mast.stsci.edu/api/v0/invoke"
TOTAL_CAP = 20_000_000


def query(obsid):
    request = {"service": "Mast.Caom.Products", "params": {"obsid": obsid}, "format": "json", "pagesize": 2000, "page": 1}
    body = urllib.parse.urlencode({"request": json.dumps(request, separators=(",", ":"))}).encode()
    req = urllib.request.Request(URL, data=body, headers={"Content-Type": "application/x-www-form-urlencoded", "User-Agent": "CodexRAISINMetadataAudit/1.0"})
    with urllib.request.urlopen(req, timeout=90) as response:
        raw = response.read(1_000_001)
    if len(raw) > 1_000_000:
        raise ValueError(("single observation metadata cap", obsid))
    result = json.loads(raw)
    if result.get("status") != "COMPLETE":
        raise ValueError(("MAST product query incomplete", obsid, result.get("status")))
    return obsid, raw, result


def main():
    start = time.monotonic()
    cone = json.loads((HERE / "caom-cone.json").read_text())
    selected = [row for row in cone["data"] if row["obs_collection"] == "HST" and row["proposal_id"] == "14216" and row["instrument_name"] == "WFC3/IR" and row["target_name"].upper() == "DES16E1DCX" and row["filters"] in ("F125W", "F160W", "detection")]
    obsids = sorted({int(row["obsid"]) for row in selected})
    if len(obsids) != 20:
        raise ValueError(("unexpected target observation count", len(obsids)))
    product_dir = HERE / "products"
    product_dir.mkdir(parents=True, exist_ok=True)
    total = (HERE / "caom-cone.json").stat().st_size
    summary = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(query, obsid): obsid for obsid in obsids}
        for future in concurrent.futures.as_completed(futures):
            obsid, raw, result = future.result()
            total += len(raw)
            if total > TOTAL_CAP:
                raise ValueError("total metadata acquisition cap")
            (product_dir / f"{obsid}.json").write_bytes(raw)
            summary.append({"obsid": obsid, "rows": len(result.get("data", [])), "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest(), "paging": result.get("paging")})
            if time.monotonic() - start > 600:
                raise TimeoutError("metadata time cap")
    summary.sort(key=lambda row: row["obsid"])
    (HERE / "observations.json").write_text(json.dumps(sorted(selected, key=lambda row: (row["t_min"], row["obsid"])), indent=2) + "\n")
    (HERE / "products-summary.json").write_text(json.dumps({"obsids": len(obsids), "product_rows": sum(row["rows"] for row in summary), "bytes_all_metadata": total, "elapsed_seconds": time.monotonic() - start, "entries": summary}, indent=2) + "\n")


if __name__ == "__main__":
    main()
