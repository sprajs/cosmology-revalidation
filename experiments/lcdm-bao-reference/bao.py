"""Ordered DESI DR2 Gaussian compression scored against one CLASS state.

This adapter uses CLASS distances and its printed, predicted drag horizon. It
does not compute an expansion history, supply a drag redshift, or combine probes.
Production Gaussian scoring belongs to the existing Irreducible SDK. These
functions only admit inputs, map emitted CLASS distances, and form residuals.
"""
import hashlib
import json
import math
import os
from pathlib import Path
import re
import stat

import theory

RELEASE = {"repository": "CobayaSampler/bao_data",
           "revision": "bb0c1c9009dc76d1391300e169e8df38fd1096db",
           "directory": "desi_bao_dr2"}
INPUTS = (
    {"role": "ordered_mean", "relative_path": "data/bao/desi_gaussian_bao_ALL_GCcomb_mean.txt",
     "bytes": 472, "sha256": "9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585"},
    {"role": "full_covariance", "relative_path": "data/bao/desi_gaussian_bao_ALL_GCcomb_cov.txt",
     "bytes": 2547, "sha256": "252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509"},
)
AXES = ((0.295, "DV_over_rs"), (0.510, "DM_over_rs"), (0.510, "DH_over_rs"),
        (0.706, "DM_over_rs"), (0.706, "DH_over_rs"),
        (0.934, "DM_over_rs"), (0.934, "DH_over_rs"),
        (1.321, "DM_over_rs"), (1.321, "DH_over_rs"),
        (1.484, "DM_over_rs"), (1.484, "DH_over_rs"),
        (2.33, "DH_over_rs"), (2.33, "DM_over_rs"))
IDS = tuple(f"desi-dr2-{i:02d}-{kind}" for i, (_, kind) in enumerate(AXES))
CASE_IDS = ("anchor", "precision", "ns-minus", "ns-plus")
MEAN_HEADER = b"# [z] [value at z] [quantity]\n"
NUMBER = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z")


def number(token):
    if type(token) is not str or len(token) > 96 or not NUMBER.fullmatch(token):
        raise ValueError("invalid numeric token")
    value = float(token)
    if not math.isfinite(value):
        raise ValueError("nonfinite numeric token")
    return value


def finite(value, name):
    if type(value) not in (float, int) or not math.isfinite(value):
        raise ValueError("nonfinite/non-numeric " + name)
    return float(value)


def read_pinned(path, expected):
    """Hash the exact bounded regular bytes parsed, with a linked-file check.

    The two concrete dataset readers supply their fixed release pins. This
    function neither acquires inputs nor accepts a replacement on hash failure.
    """
    path = Path(path)
    if not path.is_absolute() or ".." in path.parts or len(str(path).encode()) > 4096:
        raise ValueError("input path must be absolute and bounded")
    for parent in reversed(path.parents):
        if not parent.is_dir() or parent.is_symlink():
            raise ValueError("input ancestor is not a regular directory")
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    try:
        before = os.fstat(fd)
        if not stat.S_ISREG(before.st_mode) or before.st_size != expected["bytes"]:
            raise ValueError("input regular-file byte identity differs")
        chunks, count = [], 0
        while count <= expected["bytes"]:
            chunk = os.read(fd, min(65536, expected["bytes"] + 1 - count))
            if not chunk:
                break
            count += len(chunk)
            if count > expected["bytes"]:
                raise ValueError("consumed input byte limit")
            chunks.append(chunk)
        after, linked = os.fstat(fd), path.lstat()
        facts = lambda s: (s.st_dev, s.st_ino, s.st_size, s.st_mode, s.st_mtime_ns, s.st_ctime_ns)
        if facts(before) != facts(after) or facts(after) != facts(linked):
            raise ValueError("input changed during its read")
        raw = b"".join(chunks)
        digest = hashlib.sha256(raw).hexdigest()
        if len(raw) != expected["bytes"] or digest != expected["sha256"]:
            raise ValueError("consumed input hash identity differs")
        return raw, {"role": expected["role"], "path": str(path), "bytes": len(raw),
                     "sha256": digest, "device": after.st_dev, "inode": after.st_ino,
                     "mode": oct(stat.S_IMODE(after.st_mode)), "mtime_ns": after.st_mtime_ns}
    finally:
        os.close(fd)


def parse_data(mean_raw, covariance_raw, source_identities=None):
    """Pure structural parser; production read_data admits the byte pins first."""
    if type(mean_raw) is not bytes or type(covariance_raw) is not bytes:
        raise ValueError("BAO parser requires immutable bytes")
    if len(mean_raw) > 4096 or len(covariance_raw) > 16384:
        raise ValueError("BAO source byte limit")
    if mean_raw.startswith(MEAN_HEADER):
        mean_raw = mean_raw[len(MEAN_HEADER):]
    mean = [line.split() for line in mean_raw.decode("ascii").splitlines() if line.strip()]
    if len(mean) != len(AXES) or any(len(row) != 3 for row in mean):
        raise ValueError("BAO mean requires exactly 13 ordered rows")
    rows = []
    for index, (tokens, axis) in enumerate(zip(mean, AXES)):
        z, observed, kind = number(tokens[0]), number(tokens[1]), tokens[2]
        if (z, kind) != axis or observed <= 0:
            raise ValueError("BAO redshift/observable order or observation domain differs")
        rows.append({"index": index, "id": IDS[index],
                     "z": z, "kind": kind, "observed": observed})
    fields = [line.split() for line in covariance_raw.decode("ascii").splitlines() if line.strip()]
    if len(fields) != len(rows) or any(len(row) != len(rows) for row in fields):
        raise ValueError("BAO covariance must be full 13 by 13")
    covariance = [[number(x) for x in row] for row in fields]
    if any(covariance[i][i] <= 0 for i in range(len(rows))) or any(
            covariance[i][j] != covariance[j][i] for i in range(len(rows)) for j in range(i)):
        raise ValueError("BAO covariance must have positive diagonal and exact symmetry")
    return {"schema": "desi-dr2-ordered-gaussian-data/v1", "rows": rows,
            "covariance": covariance, "source_identities": source_identities}


def read_data(input_root):
    root = Path(input_root)
    if not root.is_absolute() or ".." in root.parts:
        raise ValueError("BAO input root must be absolute")
    buffers, sources = [], []
    for expected in INPUTS:
        raw, identity = read_pinned(root / expected["relative_path"], expected)
        buffers.append(raw)
        sources.append({**identity, "release": RELEASE,
                        "member": RELEASE["directory"] + "/" + Path(expected["relative_path"]).name})
    return parse_data(*buffers, source_identities=sources)


def predict(state, data, case_id="anchor"):
    """All 13 external means and rounded binary64 y-mu residuals, without a score."""
    if case_id not in CASE_IDS:
        raise ValueError("unsupported CLASS case identity")
    if set(data) != {"schema", "rows", "covariance", "source_identities"} or \
            data["schema"] != "desi-dr2-ordered-gaussian-data/v1":
        raise ValueError("BAO dataset schema")
    if type(state["parameters"].get("Omega_k")) is not str or \
            number(state["parameters"]["Omega_k"]) != 0:
        raise ValueError("this BAO distance adapter requires explicit flat CLASS Omega_k=0")
    drag = state["drag"]
    if drag.get("source") != "CLASS-thermodynamics-verbose":
        raise ValueError("BAO requires the same CLASS state's predicted drag")
    ruler = finite(drag["r_drag_Mpc"], "CLASS predicted ruler")
    if ruler <= 0:
        raise ValueError("nonpositive CLASS predicted ruler")
    if len(data["rows"]) != len(AXES):
        raise ValueError("BAO dataset row count")
    rows = []
    points = theory.background_at_many(state["background"], [x[0] for x in AXES])
    if len(points) != len(AXES):
        raise ValueError("CLASS background query count differs")
    for index, (row, axis, point) in enumerate(zip(data["rows"], AXES, points)):
        if set(row) != {"index", "id", "z", "kind", "observed"} or type(row["index"]) is not int or \
                row["index"] != index or row["id"] != IDS[index] or \
                (row["z"], row["kind"]) != axis:
            raise ValueError("BAO row identity/order differs")
        if finite(row["observed"], "BAO observation") <= 0:
            raise ValueError("nonpositive BAO observation")
        h = finite(point["H_1_Mpc"], "CLASS H/c")
        dm = finite(point["DM_Mpc"], "CLASS transverse distance")
        if point["z"] != row["z"] or h <= 0 or dm <= 0:
            raise ValueError("BAO background positivity")
        dh = 1.0 / h  # CLASS emits H/c in 1/Mpc, so this is c/H in Mpc.
        distance = (dh if row["kind"] == "DH_over_rs" else dm if row["kind"] == "DM_over_rs"
                    else (row["z"] * dm * dm * dh) ** (1.0 / 3.0))
        prediction = finite(distance / ruler, "BAO ratio")
        if prediction <= 0:
            raise ValueError("nonpositive BAO ratio")
        rows.append({**row, "predicted": prediction, "distance_Mpc": distance,
                     "residual": finite(row["observed"] - prediction, "BAO residual"),
                     "background": point})
    common_state = {"parameters": state["parameters"], "drag": drag,
                    "source_identities": state.get("source_identities")}
    state_raw = json.dumps(common_state, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    if len(state_raw) > 65536:
        raise ValueError("CLASS state identity metadata bound")
    return {"schema": "desi-dr2-class-predictions/v1", "case_id": case_id, "rows": rows,
            "common_state_sha256": hashlib.sha256(state_raw).hexdigest(),
            "r_drag_Mpc": ruler, "drag": dict(drag), "parameters": dict(state["parameters"]),
            "coordinate_measure": "product-of-13-dimensionless-DM-DH-DV-over-rdrag-coordinates",
            "input_source_identities": data["source_identities"],
            "class_source_identities": state.get("source_identities"),
            "residual_convention": "binary64-observed-minus-predicted",
            "qualification": {"input_byte_admission": "owned-by-read_data-caller",
                              "background_interpolation_error_bound": None,
                              "drag_source": "CLASS-printed-six-decimal-predicted-ruler",
                              "inference": "not-assessed", "joint_probe_target": None,
                              "score": "unavailable-until-separate-native-Gaussian-batch",
                              "scope": "fixed-parameter-means-for-released-Gaussian-BAO-compression"}}


def evaluate(state, config, case_id="anchor"):
    """Input-admitted prediction adapter; this does not evaluate a likelihood."""
    if type(config) is not dict or set(config) != {"input_root"} or type(config["input_root"]) is not str:
        raise ValueError("closed BAO runtime configuration")
    data = read_data(config["input_root"])
    result = predict(state, data, case_id)
    after = read_data(config["input_root"])
    if data["source_identities"] != after["source_identities"]:
        raise ValueError("terminal BAO input identity drift")
    result["qualification"]["input_byte_admission"] = "exact-release-bytes-before-and-after"
    return result
