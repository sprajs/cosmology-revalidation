"""Exact cached FITRES row cloning for repeated Monte Carlo inputs.

Only the occurrence CID changes. Reusing immutable numeric tokens avoids slow
reformatting and removes an otherwise unnecessary floating-point round trip.
"""

from pathlib import Path
from common import sha


def cache_table(record):
    path = Path(record["path"])
    assert sha(path) == record["sha256"]
    lines = path.read_text().splitlines()
    header = next(line for line in lines if line.startswith("VARNAMES:"))
    names = header.split()[1:]
    assert names[0] == "CID"
    rows = {}
    for line in lines:
        if not line.startswith("SN:"):
            continue
        _, cid, suffix = line.split(maxsplit=2)
        assert cid not in rows
        rows[cid] = suffix
    assert len(rows) == record["rows"]
    return {"header": header, "columns": names[1:], "rows": rows, "source": record}


def write_clones(frame, path, cache, target="literal"):
    assert list(frame.columns) == cache["columns"]
    assert frame.index.is_unique
    source = cache["source"]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as stream:
        stream.write(
            "# Paired attempt bootstrap: only CID is replaced; every other token comes from the hashed parent table.\n"
        )
        stream.write(cache["header"] + "\n")
        for cid, original in zip(frame.index, frame.ORIGINAL_CID):
            stream.write(
                "SN: " + str(cid) + " " + cache["rows"][str(int(original))] + "\n"
            )
    return {
        "path": str(path),
        "rows": len(frame),
        "sha256": sha(path),
        "target": target,
        "changed_columns": ["CID"],
        "parent_table": source,
        "mapping": "Every cloned occurrence preserves its original row tokens except CID; training-target transformation is already present in the separately validated hashed parent.",
    }
