#!/usr/bin/env python3
"""One-attempt exact MAST RAW acquisition; no FITS arrays opened."""
from __future__ import annotations

import argparse
import hashlib
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
PLAN = json.loads((HERE / "raw-acquisition-protocol.json").read_text())
STATE = HERE / "raw-acquisition.json"
FILES = HERE / "raw"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def save(state: dict) -> None:
    STATE.write_text(json.dumps(state, indent=2, sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("phase", choices=("engineering", "remaining"))
    args = parser.parse_args()
    FILES.mkdir(exist_ok=True)
    state = json.loads(STATE.read_text()) if STATE.exists() else {
        "protocol_sha256": sha(HERE / "raw-acquisition-protocol.json"),
        "downloaded_bytes": 0,
        "active_seconds": 0.0,
        "records": [],
    }
    if state["protocol_sha256"] != sha(HERE / "raw-acquisition-protocol.json"):
        raise SystemExit("frozen protocol hash changed")
    done = {r["root"] for r in state["records"]}
    phase_records = [x for x in PLAN["records"] if x["phase"] == args.phase]
    started = time.monotonic()
    for item in phase_records:
        root = item["root"]
        if root in done:
            continue  # No second attempt, regardless of prior status.
        row = {k: item[k] for k in ("root", "uri", "filename", "size", "phase")}
        row["attempt"] = 1
        row["status"] = "started"
        state["records"].append(row)
        save(state)
        if state["downloaded_bytes"] + item["size"] > PLAN["cap_bytes"]:
            row["status"] = "resource_stop_bytes"
            break
        if state["active_seconds"] + (time.monotonic() - started) > PLAN["cap_seconds"]:
            row["status"] = "resource_stop_time"
            break
        url = "https://mast.stsci.edu/api/v0.1/Download/file?" + urllib.parse.urlencode({"uri": item["uri"]})
        row["url"] = url
        path = FILES / item["filename"]
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "Codex-dark-ramp-fixed-cohort/1"}), timeout=30) as res:
                row["http_status"] = res.status
                if res.status != 200:
                    raise RuntimeError(f"unexpected HTTP {res.status}")
                row["response_length"] = res.headers.get("Content-Length")
                with path.open("xb") as f:
                    while True:
                        if state["active_seconds"] + time.monotonic() - started >= PLAN["cap_seconds"]:
                            raise RuntimeError("cumulative active-time cap")
                        if state["downloaded_bytes"] >= PLAN["cap_bytes"]:
                            raise RuntimeError("cumulative byte cap")
                        chunk = res.read(min(1 << 20, PLAN["cap_bytes"] - state["downloaded_bytes"]))
                        if not chunk:
                            break
                        f.write(chunk)
                        state["downloaded_bytes"] += len(chunk)
                        if state["downloaded_bytes"] % (8 << 20) < len(chunk):
                            save(state)
            row["actual_size"] = path.stat().st_size
            if row["actual_size"] != item["size"]:
                raise RuntimeError(f"size {row['actual_size']} != frozen {item['size']}")
            row["sha256"] = sha(path)
            row["status"] = "downloaded"
        except Exception as exc:
            row["status"] = "failed"
            row["error"] = repr(exc)
            row["partial_size"] = path.stat().st_size if path.exists() else 0
            break
        finally:
            state["active_seconds"] += time.monotonic() - started
            started = time.monotonic()
            save(state)
    state["phase_executed"] = args.phase
    save(state)
    print(json.dumps({"phase": args.phase, "active_seconds": state["active_seconds"],
                      "downloaded_bytes": state["downloaded_bytes"],
                      "statuses": {r["root"]: r["status"] for r in state["records"]}}, indent=2))


if __name__ == "__main__":
    main()
