"""Minimal official PLC3.01 C API transport, loaded only after caller admission.

Signatures and ownership come from captured clik.h 873b4c..., clik.c b81cc7...,
errorlist.h 8a133a... and errorlist.c fe6c0e.... No likelihood is implemented here.
Actual ABI/library/dependency qualification requires the separately built runtime.
"""
import ctypes as C
import hashlib
import math
from numbers import Real
import os
from pathlib import Path
import stat
import struct


class Error(C.Structure):
    pass


ErrorPointer = C.POINTER(Error)
Error._fields_ = [("errWhere", C.c_char * 2048), ("errText", C.c_char * 4192),
                  ("errValue", C.c_int), ("next", ErrorPointer)]
Parname = C.c_char * 256
_backend = None


class ClikCError(ValueError):
    def __init__(self, record):
        self.record = record
        super().__init__("official clik C error at " + record["phase"] +
                         ": " + str(record["native_code"]))


def file_identity(pin):
    """Streaming exact file admission; no dlopen before its positive result."""
    if (type(pin) is not dict or set(pin) != {"path", "bytes", "sha256"}
            or type(pin["bytes"]) is not int or pin["bytes"] < 1):
        raise ValueError("require an exact library file pin")
    path = Path(pin["path"])
    if not path.is_absolute() or ".." in path.parts or len(str(path).encode()) > 4096:
        raise ValueError("library absolute path")
    if any(p.is_symlink() or not p.is_dir() for p in path.parents):
        raise ValueError("library path ancestor")
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    try:
        before = os.fstat(fd)
        if not stat.S_ISREG(before.st_mode) or before.st_size != pin["bytes"]:
            raise ValueError("library type/size differs")
        digest, count = hashlib.sha256(), 0
        while count <= pin["bytes"]:
            chunk = os.read(fd, min(65536, pin["bytes"] + 1 - count))
            if not chunk:
                break
            digest.update(chunk); count += len(chunk)
        after, linked = os.fstat(fd), path.lstat()
        facts = lambda x: (x.st_dev, x.st_ino, x.st_size, x.st_mtime_ns, x.st_ctime_ns)
        observed = {"path": str(path), "bytes": count, "sha256": digest.hexdigest()}
        if observed != pin or facts(before) != facts(after) or facts(after) != facts(linked):
            raise ValueError("library consumed bytes/stat identity differs")
        return observed
    finally:
        os.close(fd)


class Backend:
    def __init__(self, pin):
        self.identity = file_identity(pin)
        # clik_helper.c uses dlsym(RTLD_DEFAULT) for likelihood factories.
        self.lib = C.CDLL(pin["path"], mode=os.RTLD_NOW | os.RTLD_GLOBAL)
        self.objects = []
        pp = C.POINTER(ErrorPointer)
        signatures = {
            "clik_init": ([C.c_char_p, pp], C.c_void_p),
            "clik_get_lmax": ([C.c_void_p, C.POINTER(C.c_int), pp], None),
            "clik_get_extra_parameter_names":
                ([C.c_void_p, C.POINTER(C.POINTER(Parname)), pp], C.c_int),
            "clik_compute": ([C.c_void_p, C.POINTER(C.c_double), pp], C.c_double),
            "clik_cleanup": ([C.POINTER(C.c_void_p)], None),
            "clik_get_version": ([C.c_void_p, pp], C.c_void_p),
            "_isError": ([ErrorPointer], C.c_int),
            "getErrorValue": ([ErrorPointer], C.c_int),
            "purgeError": ([pp], None),
            "free": ([C.c_void_p], None),
            "fflush": ([C.c_void_p], C.c_int),
        }
        for name, (arguments, returned) in signatures.items():
            function = getattr(self.lib, name)
            function.argtypes, function.restype = arguments, returned

    def call(self, phase, function, *arguments):
        error = ErrorPointer()  # NULL error*, but a REAL non-null error** argument.
        try:
            returned = function(*arguments, C.byref(error))
            if self.lib._isError(error):
                rows, current, seen = [], error, set()
                while current and len(rows) < 64:
                    address = C.addressof(current.contents)
                    if address in seen:
                        break
                    seen.add(address)
                    item = current.contents
                    rows.append({"code": item.errValue,
                                 "where": bytes(item.errWhere).decode("utf-8", "replace"),
                                 "text": bytes(item.errText).decode("utf-8", "replace")})
                    current = item.next
                record = {"phase": phase, "native_code": self.lib.getErrorValue(error),
                          "errors": rows, "error_chain_truncated": bool(current)}
                if phase == "clik_compute":
                    record["native_return"] = {"value": returned if math.isfinite(returned) else None,
                                               "binary64_le_hex": struct.pack("<d", returned).hex()}
                raise ClikCError(record)
            return returned
        finally:
            if error:
                self.lib.purgeError(C.byref(error))

    def version(self):
        pointer = self.call("clik_get_version", self.lib.clik_get_version, None)
        if not pointer:
            raise ValueError("null clik version")
        try:
            # clik_get_version(NULL) allocates exactly 500 bytes in pinned clik.c.
            raw = C.string_at(pointer, 500)
            end = raw.find(b"\0")
            if end < 0:
                raise ValueError("unterminated clik base version")
            return raw[:end].decode("utf-8", "strict")
        finally:
            self.lib.free(pointer)


def configure(library_pin):
    global _backend
    if _backend is not None:
        raise ValueError("clik C backend already configured")
    _backend = Backend(library_pin)
    return {"library": _backend.identity, "transport": "official-PLC3.01-ctypes-C-API",
            "ctypes_sizes": {"int": C.sizeof(C.c_int), "double": C.sizeof(C.c_double),
                             "pointer": C.sizeof(C.c_void_p), "error": C.sizeof(Error)},
            "actual_ABI_independent_qualification": None}


def version():
    if _backend is None:
        raise ValueError("clik C backend not admitted")
    return _backend.version()


def flush():
    if _backend is not None and _backend.lib.fflush(None) != 0:
        raise OSError("native stdio flush refused")


class clik:
    """Single-vector subset used by the reviewed retained PlanckPrimary owner."""
    def __init__(self, filename):
        if _backend is None:
            raise ValueError("clik C backend not admitted")
        self.backend, self.pointer, self.closed = _backend, C.c_void_p(), False
        path = Path(filename)
        if not path.is_absolute() or "\0" in str(path) or len(os.fsencode(path)) > 4096:
            raise ValueError("clik product path")
        value = self.backend.call("clik_init", self.backend.lib.clik_init, os.fsencode(path))
        if not value:
            raise ValueError("clik_init returned null without an error")
        self.pointer = C.c_void_p(value)
        self.backend.objects.append(self)
        try:
            self.maxima = self.get_lmax()
            self.names = self.get_extra_parameter_names()
            self.ndim = sum(n + 1 for n in self.maxima) + len(self.names)
        except BaseException:
            self.close()
            raise

    def _open(self):
        if self.closed or not self.pointer.value:
            raise ValueError("clik owner closed")

    def get_lmax(self):
        self._open()
        values = (C.c_int * 6)()
        self.backend.call("clik_get_lmax", self.backend.lib.clik_get_lmax, self.pointer, values)
        result = tuple(int(v) for v in values)
        if any(not -1 <= v <= 10000 for v in result):
            raise ValueError("clik multipole maximum bound")
        return result

    def get_extra_parameter_names(self):
        self._open()
        names = C.POINTER(Parname)()
        try:
            count = self.backend.call("clik_get_extra_parameter_names",
                                      self.backend.lib.clik_get_extra_parameter_names,
                                      self.pointer, C.byref(names))
            if not 0 <= count <= 256 or (count and not names):
                raise ValueError("clik nuisance count/pointer")
            result = []
            for index in range(count):
                raw = bytes(names[index])
                end = raw.find(b"\0")
                if not 0 < end <= 128:
                    raise ValueError("clik nuisance name bound")
                result.append(raw[:end].decode("utf-8", "strict"))
            if len(set(result)) != len(result):
                raise ValueError("clik duplicate nuisance name")
            return tuple(result)
        finally:
            if names:
                # Header and official C/Python callers assign caller-owned malloc.
                self.backend.lib.free(C.cast(names, C.c_void_p))

    def __call__(self, values):
        self._open()
        if type(values) not in (list, tuple) or len(values) != self.ndim:
            raise ValueError("clik exact vector length differs")
        if any(isinstance(v, bool) or not isinstance(v, Real) or not math.isfinite(v) for v in values):
            raise ValueError("clik vector must be finite real scalars")
        buffer = (C.c_double * self.ndim)(*values)
        value = self.backend.call("clik_compute", self.backend.lib.clik_compute,
                                  self.pointer, buffer)
        if not math.isfinite(value):
            raise ClikCError({"phase": "clik_compute_return_validation", "native_code": None,
                              "errors": [], "error_chain_truncated": False,
                              "native_return": {"value": None, "binary64_le_hex": struct.pack("<d", value).hex()}})
        # No sentinel is changed. The likelihood adapter refuses <= -1e30.
        return [value]

    def close(self):
        if self.closed:
            return {"status": "already_closed"}
        self.closed = True
        pointer, self.pointer = self.pointer, C.c_void_p()
        self.backend.lib.clik_cleanup(C.byref(pointer))
        return {"status": "returned", "pointer_cleared_by_C": not bool(pointer.value)}


def close_all():
    rows = []
    if _backend is not None:
        for obj in reversed(_backend.objects):
            try:
                rows.append(obj.close())
            except BaseException as exc:
                rows.append({"status": "unavailable", "kind": type(exc).__name__, "message": str(exc)[:1024]})
    return rows
