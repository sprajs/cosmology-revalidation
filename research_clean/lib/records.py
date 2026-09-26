"""Small, strict data readers and result writers shared by workflows."""

import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def jsonable(value):
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, np.generic):
        return jsonable(value.item())
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path, value):
    Path(path).write_text(json.dumps(jsonable(value), indent=2, allow_nan=False) + "\n")


def write_rows(path, rows, fields=None):
    rows = list(rows)
    if not rows and fields is None:
        raise ValueError("An empty table needs explicit column names")
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fitres(path):
    """Read FITRES without guessing identifiers, column order or missing rows."""
    opener = gzip.open if str(path).endswith(".gz") else open
    columns, rows = None, []
    with opener(path, "rt") as stream:
        for line in stream:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == "VARNAMES:":
                columns = parts[1:]
            elif parts[0] == "SN:":
                if columns is None or len(parts) != len(columns) + 1:
                    raise ValueError(f"Invalid FITRES row in {path}")
                rows.append(parts[1:])
    if columns is None or not rows:
        raise ValueError(f"No FITRES data in {path}")
    frame = pd.DataFrame(rows, columns=columns)
    if "CID" not in frame or frame.CID.duplicated().any():
        raise ValueError(f"Missing or duplicate physical identifiers in {path}")
    for name in frame.columns:
        if name == "CID":
            continue
        try:
            frame[name] = pd.to_numeric(frame[name])
        except ValueError:
            pass
    return frame.set_index("CID")
