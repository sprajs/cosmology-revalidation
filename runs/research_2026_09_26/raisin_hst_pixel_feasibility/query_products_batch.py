#!/usr/bin/env python3
"""One bounded MAST product-metadata query for frozen 14216 observation IDs."""
import hashlib
import json
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
cone = json.loads((HERE / "caom-cone.json").read_text())
selected = [r for r in cone["data"] if r["obs_collection"] == "HST" and r["proposal_id"] == "14216" and r["instrument_name"] == "WFC3/IR" and r["target_name"].upper() == "DES16E1DCX" and r["filters"] in ("F125W", "F160W", "detection")]
ids = sorted({int(r["obsid"]) for r in selected})
assert len(ids) == 20
request = {"service": "Mast.Caom.Products", "params": {"obsid": ",".join(map(str, ids))}, "format": "json", "pagesize": 2000, "page": 1}
body = urllib.parse.urlencode({"request": json.dumps(request, separators=(",", ":"))}).encode()
req = urllib.request.Request("https://mast.stsci.edu/api/v0/invoke", data=body, headers={"Content-Type": "application/x-www-form-urlencoded", "User-Agent": "CodexRAISINMetadataAudit/1.0"})
with urllib.request.urlopen(req, timeout=90) as response:
    raw = response.read(19_900_000)
if len(raw) >= 19_900_000:
    raise ValueError("metadata cap")
d = json.loads(raw)
if d.get("status") != "COMPLETE":
    raise ValueError((d.get("status"), d.get("msg")))
(HERE / "products-batch.json").write_bytes(raw)
(HERE / "products-batch-request.json").write_text(json.dumps(request, indent=2) + "\n")
(HERE / "products-batch-summary.json").write_text(json.dumps({"obsids": ids, "rows": len(d.get("data", [])), "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest(), "paging": d.get("paging")}, indent=2) + "\n")
