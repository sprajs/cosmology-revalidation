"""Bounded file/JSON transport for this SN working target only."""
import base64
import hashlib
import json
import math
import os
from pathlib import Path
import re
import stat
import struct
import sys
import zlib

TARGET = "released-sn-symmetric-working-profile/v1"
DERIVATION = "selected-binary64-rational-pair-projection-RN-even/v1"
N = 1657
CASES = ("anchor", "precision", "ns-minus", "ns-plus")
NUMBER = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z")
MAX_FILE = 128 * 1024**2
MAX_STORE = 256 * 1024**2
PORTS = {"assemble.py", "controller.py", "consumer.cpp", "reference.py", "transport.py", "controls.py"}
REQUEST_KEYS = {"schema","target_id","derivation_id","selected_n","execution_admitted","ancestry_authority","selection",
                "source_ports","reference_runtime","reference_python","structural_python","sdk","native_loader_cache","criteria","resources","history"}

def need(ok, message):
    if not ok:
        raise ValueError(message)

def pairs(items):
    result = {}
    for k, v in items:
        need(k not in result, "duplicate JSON key")
        result[k] = v
    return result

def document(raw):
    return json.loads(raw, object_pairs_hook=pairs, parse_constant=lambda _: (_ for _ in ()).throw(ValueError("nonfinite JSON")))

def encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()

def runtime_json_bytes(value):
    raw=bytearray()
    encoder=json.JSONEncoder(sort_keys=True,separators=(",",":"),allow_nan=False)
    for token in encoder.iterencode(value):
        block=token.encode("utf-8")
        need(len(raw)+len(block)+1<=4*1024**2,"runtime decoded JSON allocation4MiB")
        raw.extend(block)
    raw.extend(b"\n")
    return bytes(raw)

def unpack_runtime(envelope):
    need(type(envelope) is dict and set(envelope)=={"schema","encoding","decoded_bytes","decoded_sha256","blob"},"closed SN runtime envelope")
    need(envelope["schema"]=="SN-spectral-runtime-packed/v1" and envelope["encoding"]=="canonical-JSON-UTF8-newline/zlib9/base64","exact runtime encoding")
    count=envelope["decoded_bytes"];digest=envelope["decoded_sha256"];blob=envelope["blob"]
    need(type(count) is int and 0<count<=4*1024**2 and type(digest) is str and re.fullmatch("[0-9a-f]{64}",digest),"typed decoded identity")
    need(type(blob) is str and blob.isascii() and len(blob)<=256*1024,"bounded ASCII runtime blob")
    need(len(encoded(envelope))<=256*1024,"stored runtime fingerprint256KiB")
    compressed=base64.b64decode(blob,validate=True)
    need(base64.b64encode(compressed).decode("ascii")==blob,"canonical base64")
    decoder=zlib.decompressobj()
    raw=decoder.decompress(compressed,count+1)
    need(len(raw)==count and decoder.eof and not decoder.unused_data and not decoder.unconsumed_tail,"exact bounded zlib stream/no trailing or unconsumed data")
    need(hashlib.sha256(raw).hexdigest()==digest,"decoded runtime SHA")
    need(zlib.compress(raw,9)==compressed,"canonical zlib9 representation")
    value=document(raw)
    need(runtime_json_bytes(value)==raw,"canonical duplicate-free finite runtime JSON")
    return value

def pack_runtime(value):
    raw=runtime_json_bytes(value)
    envelope={"schema":"SN-spectral-runtime-packed/v1","encoding":"canonical-JSON-UTF8-newline/zlib9/base64",
              "decoded_bytes":len(raw),"decoded_sha256":hashlib.sha256(raw).hexdigest(),
              "blob":base64.b64encode(zlib.compress(raw,9)).decode("ascii")}
    need(len(encoded(envelope))<=256*1024,"stored runtime fingerprint256KiB")
    need(unpack_runtime(envelope)==value,"lossless runtime pack roundtrip")
    return envelope

def facts(s):
    return [s.st_dev, s.st_ino, s.st_size, s.st_mode, s.st_mtime_ns, s.st_ctime_ns]

class SourceIdentityError(ValueError):
    def __init__(self,message,evidence):
        super().__init__(message)
        self.evidence=evidence

class RequestAdmissionError(ValueError):
    def __init__(self,message,pin,expected_sha):
        super().__init__(message)
        self.request_identity={"expected_sha256":expected_sha,"consumed_complete_request":pin,"request_admitted":False}

def failure(exc):
    return {"type":type(exc).__name__,"message":str(exc)[:512],
            "source_identity":exc.evidence if isinstance(exc,SourceIdentityError) else None,"request_identity":getattr(exc,"request_identity",None)}

def read(pin, keep=True, limit=256*1024**2, allow_unpinned=False, unknown_size=False):
    authority={"path":pin.get("path"),"bytes":pin.get("bytes"),"sha256":pin.get("sha256")}
    if type(authority["path"]) is str:authority["path"]=authority["path"][:4096]
    else:authority["path"]=None
    if type(authority["sha256"]) is str:authority["sha256"]=authority["sha256"][:64]
    else:authority["sha256"]=None
    if type(authority["bytes"]) is not int or not 0<=authority["bytes"]<=256*1024**2:authority["bytes"]=None
    evidence={"schema":"SN-consumed-source-read/v1","expected_authority":authority,
              "descriptor_opened":False,"bytes_consumed":None,"consumed_sha256":None,"eof_observed":None,
              "descriptor_before":None,"descriptor_after":None,"linked_before":None,"linked_after":None,
              "observation_errors":[],"accepted_complete_identity":False}
    fd=None;h=hashlib.sha256();size=0;read_failed=False
    try:
        need(type(pin["path"]) is str and len(pin["path"].encode())<=4096 and ((unknown_size and pin["bytes"] is None) or (type(pin["bytes"]) is int and 0<=pin["bytes"]<=256*1024**2)) and ((allow_unpinned and pin["sha256"] is None) or (type(pin["sha256"]) is str and re.fullmatch("[0-9a-f]{64}",pin["sha256"]))),"bounded typed source authority")
        p = Path(pin["path"])
        need(p.is_absolute() and ".." not in p.parts, "absolute bounded source path")
        need(all(x.is_dir() and not x.is_symlink() for x in p.parents), "source parent")
        evidence["linked_before"]=facts(p.lstat())
        fd = os.open(p, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
        evidence.update(descriptor_opened=True,bytes_consumed=0,consumed_sha256=h.hexdigest(),eof_observed=False)
        before = os.fstat(fd)
        evidence["descriptor_before"]=facts(before)
        need(stat.S_ISREG(before.st_mode) and (pin["bytes"] is None or before.st_size == pin["bytes"]), "source regular/bytes")
        need(before.st_size<=limit,"source role bound before allocation")
        expected_size=before.st_size if pin["bytes"] is None else pin["bytes"]
        chunks = []
        while True:
            b = os.read(fd, min(65536, expected_size + 1 - size,limit+1-size))
            if not b:
                evidence["eof_observed"]=True
                break
            size += len(b)
            h.update(b)
            evidence.update(bytes_consumed=size,consumed_sha256=h.hexdigest())
            need(size <= expected_size and size<=limit, "source bound")
            if keep:
                chunks.append(b)
        after = os.fstat(fd)
        evidence["descriptor_after"]=facts(after)
        evidence["linked_after"]=facts(p.lstat())
        need(facts(before) == facts(after) == evidence["linked_before"] == evidence["linked_after"], "source drift")
        need(size == expected_size and (pin["sha256"] is None or h.hexdigest() == pin["sha256"]), "source digest")
        evidence["accepted_complete_identity"]=True
        return b"".join(chunks) if keep else None, {"path": str(p), "bytes": size, "sha256": h.hexdigest(), "stat": facts(after)}
    except BaseException as exc:
        read_failed=True
        evidence["predicate_error"]={"type":type(exc).__name__,"message":str(exc)[:512]}
        if fd is not None:
            try:evidence["descriptor_failure_final"]=facts(os.fstat(fd))
            except BaseException as err:evidence["observation_errors"].append("descriptor-after:"+type(err).__name__)
        if "p" in locals():
            try:evidence["linked_failure_final"]=facts(p.lstat())
            except BaseException as err:evidence["observation_errors"].append("linked-after:"+type(err).__name__)
        raise SourceIdentityError(str(exc)[:512],evidence) from exc
    finally:
        if fd is not None:
            try:os.close(fd)
            except BaseException as exc:
                evidence["observation_errors"].append("descriptor-close:"+type(exc).__name__)
                if not read_failed:raise SourceIdentityError("descriptor close failed",evidence) from exc

def owned_json_bytes(path,limit):
    path=Path(path)
    return read({"path":str(path),"bytes":None,"sha256":None},limit=limit,allow_unpinned=True,unknown_size=True)

def identity(path):
    path = Path(path)
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(65536), b""):
            h.update(b)
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": h.hexdigest()}

def store_bytes(root):
    total = sum(p.stat().st_size for p in root.rglob("*") if p.is_file())
    need(total <= MAX_STORE, "SN owned256MiB store")
    return total

def emit(path, raw, limit=MAX_FILE):
    need(len(raw) <= limit, "individual output bound")
    with path.open("xb") as f:
        need(f.write(raw) == len(raw), "short output")
        f.flush()
        os.fsync(f.fileno())
    path.chmod(0o444)
    return identity(path)

def number(token):
    need(type(token) is str and len(token) <= 96 and NUMBER.fullmatch(token), "finite numeric token grammar")
    x = float(token)
    need(math.isfinite(x), "nonfinite source number")
    return x

def f64_payload(path, values):
    with path.open("xb") as f:
        for start in range(0, len(values), 4096):
            raw = struct.pack(">" + "d" * min(4096, len(values)-start), *values[start:start+4096])
            need(f.write(raw) == len(raw), "short binary64 output")
        f.flush()
        os.fsync(f.fileno())
    path.chmod(0o444)
    return identity(path)

def request_fields(raw):
    request = document(raw)
    need(type(request) is dict and set(request)==REQUEST_KEYS and request["schema"]=="SN-S_hat-request-proposal/v1", "closed request schema")
    need(request["target_id"] == TARGET and request["derivation_id"] == DERIVATION and request["selected_n"] == N,
         "exact new SN target/derivation/dimension")
    need(request["execution_admitted"] is True, "proposal is not execution admission")
    need(type(request["source_ports"]) is dict and set(request["source_ports"]) == PORTS, "exact six new source ports")
    need(request["criteria"] == {"native_sensitivity":"1e-8","spectral_scaled_gate":"1e-8","stationarity":"1e-8","same_target_relative_profile_absolute":"1e-6","synthetic_scalar_and_lane_absolute":"1e-6"}, "fixed engineering criteria")
    need(request["resources"] == {"family_wall_seconds":1800,"assembly_wall_seconds":120,"compile_wall_seconds":120,"runtime_wall_seconds":30,"controls_family_wall_seconds":840,"major_control_wall_seconds":180,"minor_control_wall_seconds":20,"actual_native_wall_seconds":300,"actual_reference_wall_seconds":300,"terminal_reservation_seconds":60,"address_bytes":2147483648,"file_bytes":134217728,"store_bytes":268435456,"logs_family_bytes":16777216,"json_bytes":4194304,"jobs":1,"threads":1,"direct_children_maximum":16}, "exact distinct proposed phase bounds")
    need(all(type(v) is int for v in request["resources"].values()) and type(request["selected_n"]) is int,"integers exclude bool/float")
    sdk=request["sdk"]
    need(set(sdk)=={"revision","build_id","include_directory","archive","compiler","admission","build_manifest","headers","libraries"} and sdk["revision"]=="f3b19b6539c72fda5a8b489a0546ca4d0266ace3" and sdk["build_id"]=="6d495efb166006c6ce651359a366af5d87686ce516ecd7f0f49eed4e8d8be07e" and sdk["archive"]["sha256"]=="be08158b14c77f7692fe2a55caf9aaac4abe97eef1edd75e4966e57f27516d77" and sdk["admission"]["sha256"]=="406675b96e741b8bd1aee6d1b9008c581d47f413442e945357334c38d9dc3ba1" and len(sdk["headers"])==3 and len(sdk["libraries"])==5,"closed unchanged historical SDK roles")
    cache=request["native_loader_cache"]
    need(type(cache) is dict and set(cache)=={"path","bytes","sha256"} and type(cache["path"]) is str and cache["path"]=="/etc/ld.so.cache" and type(cache["bytes"]) is int and cache["bytes"]==175035 and type(cache["sha256"]) is str and cache["sha256"]=="5c863dcc3df28b7ecb622e4d7e8b494eacce37bc8109bbfcb525b0c69a700184","exact separately pinned native loader cache")
    for name, source in request["source_ports"].items():
        need(Path(source["path"]).name==name and Path(source["path"]).parent==Path(__file__).resolve().parent, "one actual source packet")
        read(source, False)
    need(Path(request["source_ports"]["transport.py"]["path"]).resolve()==Path(__file__).resolve() and Path(sys.argv[0]).resolve()==Path(request["source_ports"][Path(sys.argv[0]).name]["path"]).resolve(), "executed source identity")
    return request

def admitted_request(path,expected_sha):
    raw,pin=read({"path":str(Path(path)),"bytes":None,"sha256":expected_sha},limit=65536,unknown_size=True)
    try:
        return request_fields(raw),pin
    except SourceIdentityError as exc:
        exc.request_identity={"expected_sha256":expected_sha,"consumed_complete_request":pin,"request_admitted":False}
        raise
    except BaseException as exc:
        raise RequestAdmissionError(str(exc)[:512],pin,expected_sha) from exc
