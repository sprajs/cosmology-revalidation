#!/usr/bin/env python3
"""Bounded transport fork of the original 60-digit reference; no native calls.

Original scientific ancestry SHA256:
f912a6f3b0b351b91ece9d26af959200618776908d1842bb64caa63ed87ebe0e.
The transport/runtime path is separate from scientific initialization.
"""
import argparse
import codecs
import decimal
import encodings.ascii
import encodings.latin_1
import encodings.utf_8
import hashlib
import importlib.machinery
import importlib.metadata
import importlib.util
import json
import math
import os
import pathlib
import platform
import re
import resource
import stat
import sys
import time
import mpmath as mp
import mpmath.libmp.backend as mp_backend

mp.mp.dps = 60
BASE = pathlib.Path(__file__).resolve().parent
evals = {"time_frequency": 0, "characteristic_function": 0}
_PARTIAL = None


def _initialize_scientific_constants():
    """Original initializers, in original order, reached only after admission."""
    global mpf, ARITH, ZERO, ONE, weights, T, C, h, c, r, gain, qe
    mpf = mp.mpf
    ARITH = mpf("1e-50")
    ZERO = mpf(0)
    ONE = mpf(1)
    weights = (mpf(1)/4, mpf(3)/4)
    T = ((mpf(1)/4, mpf(3)/4), (mpf(3)/4, mpf(1)/4))
    C = mpf(2)**-80
    h = mpf("6.62607015e-34")
    c = mpf(299792458)
    r = mpf(2)
    gain = mpf(2)
    qe = (mpf(3)/4, mpf(1)/2)


def reserve(x):
    return ARITH*(1+abs(x))


def text(x):
    return mp.nstr(x, 60)


def rule(n):
    nodes, weights_ = mp.gauss_quadrature(n, "legendre")
    return list(zip(nodes, weights_))


rules = {}


def integrate(f, lo, hi, n):
    middle = (hi+lo)/2
    half = (hi-lo)/2
    return half*mp.fsum(w*f(middle+half*x) for x, w in rules[n])


def source_luminosity(nu, t):
    # Direct original analytic source expression, no native interpolator.
    tau = (t-10)/r
    wavelength_rest = c/(r*nu)
    x = wavelength_rest-1
    return C*(1+tau/4+x/2+tau*x/16)


def frequency_photons(e, b, n):
    lower, upper = ((mpf(9), mpf(11)), (mpf(12), mpf(15)))[e]
    # Mapped source knots are8,10,14; declared zero outside8..14.
    lo, hi = max(lower, mpf(8)), min(upper, mpf(14))
    splits = sorted(set([lo, hi]+[mpf(v) for v in (8, 10, 14) if lo<v<hi]))
    band_lower, band_upper = ((mpf(2), mpf(3)), (mpf(3), mpf(4)))[b]
    nu_lower, nu_upper = c/band_upper, c/band_lower
    def integrand(nu, t):
        evals["time_frequency"] += 1
        return source_luminosity(nu, t)*c/nu**3
    total = mp.fsum(integrate(lambda t: integrate(lambda nu: integrand(nu, t), nu_lower, nu_upper, n), a, z, n)
                    for a, z in zip(splits, splits[1:]))
    return total/(4*mp.pi*r*h)


def poisson(k, lam):
    return mp.exp(-lam)*lam**k/mp.factorial(k)


def panel_cache(n):
    sigma = mpf(3)/4
    panels = []
    for i in range(96):
        lo, hi = mpf(i)/4, mpf(i+1)/4
        middle, half = (hi+lo)/2, (hi-lo)/2
        for node, weight in rules[n]:
            t = middle+half*node
            panels.append((t, half*weight, mp.exp(1j*t)-1,
                           mp.exp(-sigma**2*t**2/2),
                           {z: mp.exp(-1j*t*z) for z in (1, 4, 5)}))
    return panels


def fourier(lam, z, panels):
    density_terms = []
    cdf_terms = []
    for t, weight, ep, envelope, phase in panels:
        # Shared CF evaluation supplies density at z and thresholdCDF at1.
        phi = mp.exp(lam*ep)*envelope
        evals["characteristic_function"] += 1
        density_terms.append(weight*mp.re(phase[z]*phi))
        cdf_terms.append(weight*mp.im(phase[1]*phi)/t)
    density_adu = gain/mp.pi*mp.fsum(density_terms)
    below_threshold = mpf(1)/2-mp.fsum(cdf_terms)/mp.pi
    return density_adu, below_threshold


def mixture(values, errors):
    value = mp.fsum(weights[s]*values[s] for s in range(2))
    error = mp.fsum(weights[s]*errors[s] for s in range(2))+reserve(value)
    if not value>error>=0:
        raise ArithmeticError(json.dumps({"stage":"positive-mixture-error-gate","value":text(value),"error":text(error),
                                         "conditional_values":[text(v) for v in values],"conditional_errors":[text(e) for e in errors]}))
    return value, error, -mp.log1p(-error/value)


def product(values, errors):
    value = mp.fprod(values)
    error = mp.fprod(v+e for v, e in zip(values, errors))-value+reserve(value)
    return value, error


def run():
    started = time.perf_counter()
    for n in (32, 64):
        rules[n] = rule(n)
    fractions = ((mpf(17)/12, mpf(29)/12), (mpf(23)/12, mpf(307)/96))
    photons = []
    photon_errors = []
    frequency_rows = []
    base_frequency = {}
    for e in range(2):
        for b in range(2):
            q32, q64 = (frequency_photons(e, b, n) for n in (32, 64))
            base_frequency[e, b] = (q32, q64)
    for s in range(2):
        for e in range(2):
            for b in range(2):
                analytic = r*C*T[s][b]*fractions[e][b]/(4*mp.pi*h*c)
                q32, q64 = (q*T[s][b] for q in base_frequency[e, b])
                analytic_error = reserve(analytic)
                q32_error = reserve(q32)
                q64_error = abs(q64-q32)+reserve(q64)+q32_error
                photons.append(analytic)
                photon_errors.append(analytic_error)
                frequency_rows.append({"polynomial": text(analytic), "polynomial_error": text(analytic_error),
                                       "frequency32": text(q32), "frequency64": text(q64),
                                       "frequency32_error": text(q32_error), "frequency64_error": text(q64_error)})
                _observe_partial("photons", locals())
    frequency_seconds = time.perf_counter()-started
    _observe_partial("frequency", locals())
    sigma = mpf(3)/4
    cutoff = mpf(24)
    density_tail = gain*mp.exp(-sigma**2*cutoff**2/2)/(mp.pi*sigma**2*cutoff)
    cdf_tail = mp.exp(-sigma**2*cutoff**2/2)/(mp.pi*sigma**2*cutoff**2)
    caches = {n: panel_cache(n) for n in (32, 64)}
    results = []
    detector_rows = []
    lam = [qe[i%2]*photons[i]+mpf(1)/4+mpf(1)/8*(2 if i%4<2 else 3) for i in range(8)]
    _observe_partial("lambda", locals())
    for noise in range(2):
        densities = []
        below = []
        density_errors = []
        below_errors = []
        sigma_rows = []
        _observe_partial("sigma-start", locals())
        for i in range(8):
            k = (2, 3, 3, 2)[i%4]
            if noise:
                z = (4, 5, 5, 4)[i%4]
                d32, b32 = fourier(lam[i], z, caches[32])
                d64, b64 = fourier(lam[i], z, caches[64])
                de = abs(d64-d32)+density_tail+reserve(d64)+reserve(d32)
                be = abs(b64-b32)+cdf_tail+reserve(b64)+reserve(b32)
                density, nondetection = d64, b64
                sigma_rows.append({"density32_per_adu": text(d32), "density64_per_adu": text(d64),
                                   "density_error_per_adu": text(de), "below32": text(b32),
                                   "below64": text(b64), "below_error": text(be),
                                   "density_tail_per_adu": text(density_tail), "cdf_tail": text(cdf_tail)})
            else:
                density, nondetection = poisson(k, lam[i]), mp.exp(-lam[i])
                de, be = reserve(density), reserve(nondetection)
                sigma_rows.append({"count_probability": text(density), "count_probability_error": text(de),
                                   "below": text(nondetection), "below_error": text(be)})
            if not (density>de and nondetection>be and 1-nondetection>be):
                raise ArithmeticError(json.dumps({"stage":"conditional-positive-reference-gate","sigma_index":noise,"channel_index":i,
                                                 "lambda":text(lam[i]),"density":text(density),"density_error":text(de),"nondetection":text(nondetection),
                                                 "nondetection_error":text(be),"coarse_refined_attempt":sigma_rows[-1]}))
            densities.append(density);below.append(nondetection)
            density_errors.append(de);below_errors.append(be)
            _observe_partial("sigma-row", locals())
        detector_rows.append({"sigma_electrons": noise*.75, "rows": sigma_rows})
        _observe_partial("sigma-complete", locals())
        for record in range(3):
            conditional = []
            event_conditional = []
            conditional_errors = []
            event_errors = []
            per_exposure = [[], []]
            per_exposure_event = [[], []]
            per_exposure_errors = [[], []]
            per_exposure_event_errors = [[], []]
            for s in range(2):
                start = 4*s
                data = densities[start:start+4].copy()
                data_errors = density_errors[start:start+4].copy()
                if record==2:
                    data[3]=below[start+3];data_errors[3]=below_errors[start+3]
                v, ve = product(data, data_errors)
                q, qe_ = product([1-below[start+j] for j in range(4)], below_errors[start:start+4])
                conditional.append(v);conditional_errors.append(ve)
                event_conditional.append(q);event_errors.append(qe_)
                for e in range(2):
                    pv, pe = product(data[2*e:2*e+2], data_errors[2*e:2*e+2])
                    pq, pqe = product([1-below[start+j] for j in range(2*e,2*e+2)], below_errors[start+2*e:start+2*e+2])
                    per_exposure[e].append(pv);per_exposure_errors[e].append(pe)
                    per_exposure_event[e].append(pq);per_exposure_event_errors[e].append(pqe)
            v, ve, vlog_error = mixture(conditional, conditional_errors)
            q, qe_, qlog_error = mixture(event_conditional, event_errors)
            value_log = mp.log(v)-(mp.log(q) if record==1 else 0)
            joint_log_error=vlog_error+reserve(mp.log(v))
            event_log_error=qlog_error+reserve(mp.log(q))
            value_error=joint_log_error+(event_log_error if record==1 else 0)+reserve(value_log)
            independent_v_parts=[mixture(per_exposure[e], per_exposure_errors[e]) for e in range(2)]
            independent_q_parts=[mixture(per_exposure_event[e], per_exposure_event_errors[e]) for e in range(2)]
            independent_v, independent_ve=product([a[0] for a in independent_v_parts],[a[1] for a in independent_v_parts])
            independent_q, independent_qe=product([a[0] for a in independent_q_parts],[a[1] for a in independent_q_parts])
            independent_vlog_error=-mp.log1p(-independent_ve/independent_v)+reserve(mp.log(independent_v))
            independent_qlog_error=-mp.log1p(-independent_qe/independent_q)+reserve(mp.log(independent_q))
            def relative_error(a, ae, b, be):
                ratio=a/b
                return max((a+ae)/(b-be)-ratio,ratio-(a-ae)/(b+be))+reserve(ratio-1)
            record_relative_error=relative_error(v,ve,independent_v,independent_ve)
            event_relative_error=relative_error(q,qe_,independent_q,independent_qe)
            results.append({"sigma_index": noise, "record_index": record,
                            "conditional_record": [text(v_) for v_ in conditional],
                            "conditional_event": [text(v_) for v_ in event_conditional],
                            "log_joint": text(mp.log(v)), "joint_log_error": text(joint_log_error),
                            "log_event": text(mp.log(q)), "event_log_error": text(event_log_error),
                            "log_value": text(value_log), "value_log_error": text(value_error),
                            "product_per_exposure_log_record": text(mp.log(independent_v)),
                            "product_per_exposure_log_event": text(mp.log(independent_q)),
                            "product_per_exposure_record_log_error": text(independent_vlog_error),
                            "product_per_exposure_event_log_error": text(independent_qlog_error),
                            "record_relative_difference": text(v/independent_v-1),
                            "event_relative_difference": text(q/independent_q-1),
                            "record_signed_difference":text(v-independent_v),"record_signed_difference_error":text(ve+independent_ve+reserve(v-independent_v)),
                            "event_signed_difference":text(q-independent_q),"event_signed_difference_error":text(qe_+independent_qe+reserve(q-independent_q)),
                            "record_relative_difference_error":text(record_relative_error),"event_relative_difference_error":text(event_relative_error),
                            "positive_joint": text(v), "positive_joint_error": text(ve),
                            "positive_event": text(q), "positive_event_error": text(qe_)})
            _observe_partial("record", locals())
    return {"schema": "original-high-precision-temporal-optical-reference/v1",
            "mpmath_version": mp.__version__, "decimal_digits": mp.mp.dps,
            "shared_ancestry": ["declared equations and exact SI constants", "mpmath arithmetic/GL node generator; no native kernels"],
            "frequency_seconds": frequency_seconds, "elapsed_seconds": time.perf_counter()-started,
            "integrand_evaluations": evals, "photons": [text(v) for v in photons],
            "photon_errors": [text(v) for v in photon_errors], "frequency_rows": frequency_rows,
            "lambda_electrons": [text(v) for v in lam], "detector_rows": detector_rows, "records": results}


# All helpers below implement the frozen transport contract, not reference physics.
JSON_LIMIT = 1048576
INVENTORY_LIMIT = 4096
INVENTORY_BYTES = 67108864
CONTRACT_SHA = "e3a1ceb57d96f599cf82d938462d03feef5d424ad9d2c70f82316d63c34ba6d7"
ADMISSION_SHA = "1ee77abd937be8342f3be946c17453a162f75e9cd315a9a53b05f8fe2db401cc"
ORIGINAL_REFERENCE_SHA = "f912a6f3b0b351b91ece9d26af959200618776908d1842bb64caa63ed87ebe0e"
PYTHON_SHA = "815b1275bf87e7595fe1cdf3f2efc6ae58e27b842262a7b530de8610ac2760ba"
CFG_SHA = "11a0241468a7d4cb44419cd7c7a87cd6149c3e90e38b75be6e1a201279621e43"
MPMATH_INIT_SHA = "b241584d2c1fc0304b0a1015ea923749d7b0800411dd406dcab7c82bf25d9fe8"
LEGACY_SHA = "982dc09f2ddf1206cc9cd4467c7be1ce4238a612489b43bf5f5d196a8b2495a8"
PYTHON_PATH = "/home/szymon/Projects/irreducible/evidence/project-review/science/recovery-abundance-growth-20261002/reference-venv/bin/python"
RESOLVED_PYTHON = "/home/szymon/.local/share/mise/installs/python/3.14.8/bin/python3.14"
PYTHON_VERSION = "3.14.8 (main, Oct  1 2026, 20:55:04) [GCC 16.2.1 20260810]"
VENV_PREFIX = "/home/szymon/Projects/irreducible/evidence/project-review/science/recovery-abundance-growth-20261002/reference-venv"
BASE_PREFIX = "/home/szymon/.local/share/mise/installs/python/latest"
ENV_KEYS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "PYTHONDONTWRITEBYTECODE")
SHA_PATTERN = re.compile(r"[0-9a-f]{64}\Z", re.ASCII)
DECIMAL_PATTERN = re.compile(r"-?(0|[1-9][0-9]*)(\.[0-9]+)?([eE][+-]?[0-9]+)?\Z", re.ASCII)


class TransportRefusal(ValueError):
    def __init__(self, stage, reason, kind="admission", path=None):
        super().__init__(reason)
        self.stage, self.kind, self.path = stage, kind, path


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf8")


def _sha(raw):
    return hashlib.sha256(raw).hexdigest()


def _finite_decimal(token):
    if len(token) > 256 or DECIMAL_PATTERN.fullmatch(token) is None:
        raise TransportRefusal("json-number", "Invalid or excessive decimal token")
    value = decimal.Decimal(token)
    if not value.is_finite() or abs(value.adjusted()) > 1000:
        raise TransportRefusal("json-number", "Nonfinite or excessive decimal exponent")
    return value


def _integer_token(token):
    if len(token) > 256:
        raise TransportRefusal("json-number", "Excessive integer token")
    return int(token)


def _pairs(items):
    out = {}
    for key, value in items:
        if key in out:
            raise TransportRefusal("json-keys", "Duplicate JSON key: " + key[:256])
        out[key] = value
    return out


def _tree_bounds(value, depth=0, count=None):
    if count is None:
        count = [0]
    count[0] += 1
    if depth > 32 or count[0] > 65536:
        raise TransportRefusal("json-bounds", "JSON depth/node limit exceeded")
    if type(value) is str:
        if len(value.encode("utf8")) > 8192:
            raise TransportRefusal("json-bounds", "JSON string limit exceeded")
    elif type(value) is dict:
        for key, child in value.items():
            _tree_bounds(key, depth + 1, count)
            _tree_bounds(child, depth + 1, count)
    elif type(value) is list:
        for child in value:
            _tree_bounds(child, depth + 1, count)
    elif type(value) is decimal.Decimal:
        if not value.is_finite() or abs(value.adjusted()) > 1000:
            raise TransportRefusal("json-bounds", "JSON decimal limit exceeded")
    elif value is not None and type(value) not in (bool, int):
        raise TransportRefusal("json-type", "Unsupported parsed JSON type")


def strict_json(raw):
    if type(raw) is not bytes or not raw or len(raw) > JSON_LIMIT:
        raise TransportRefusal("json-bounds", "JSON byte limit or empty input")
    decoded = raw.decode("utf8", errors="strict")
    # Reject deep syntax before the decoder can allocate a deeply nested object.
    depth = 0
    quoted = escaped = False
    for char in decoded:
        if quoted:
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == '"':
                quoted = False
        elif char == '"':
            quoted = True
        elif char in "[{":
            depth += 1
            if depth > 32:
                raise TransportRefusal("json-bounds", "JSON depth limit exceeded")
        elif char in "]}":
            depth -= 1
    def bad_constant(token):
        raise TransportRefusal("json-number", "Nonfinite JSON constant: " + token)
    result = json.loads(decoded, object_pairs_hook=_pairs, parse_float=_finite_decimal,
                        parse_int=_integer_token, parse_constant=bad_constant)
    _tree_bounds(result)
    return result


def _json_bytes(value):
    """Bounded exact Decimal tokens; floating timing metadata never enters math."""
    chunks, byte_count, nodes = [], 0, 0
    def add(part):
        nonlocal byte_count
        byte_count += len(part.encode("utf8"))
        if byte_count > JSON_LIMIT:
            raise TransportRefusal("output-size", "Reference JSON exceeds 1 MiB", "io")
        chunks.append(part)
    def emit(item, depth=0):
        nonlocal nodes
        nodes += 1
        if depth > 32 or nodes > 65536:
            raise TransportRefusal("output-size", "Reference JSON structure exceeds bounds", "io")
        if item is None:
            add("null")
        elif type(item) is bool:
            add("true" if item else "false")
        elif type(item) is str:
            if len(item.encode("utf8")) > 8192:
                raise TransportRefusal("output-size", "Reference JSON string exceeds bounds", "io")
            add(json.dumps(item, ensure_ascii=False))
        elif type(item) is int:
            token = str(item)
            if len(token) > 256:
                raise TransportRefusal("output-number", "Excessive output integer", "io")
            add(token)
        elif type(item) is decimal.Decimal:
            _finite_decimal(str(item))
            add(str(item))
        elif type(item) is float:
            if not math.isfinite(item):
                raise TransportRefusal("output-number", "Nonfinite timing/resource metadata", "io")
            add(json.dumps(item, allow_nan=False))
        elif type(item) is dict:
            add("{")
            for i, key in enumerate(sorted(item)):
                if type(key) is not str:
                    raise TransportRefusal("output-type", "Nonstring JSON key", "io")
                if i:
                    add(",")
                emit(key, depth + 1); add(":"); emit(item[key], depth + 1)
            add("}")
        elif type(item) is list:
            add("[")
            for i, child in enumerate(item):
                if i:
                    add(",")
                emit(child, depth + 1)
            add("]")
        else:
            raise TransportRefusal("output-type", "Unsupported output type", "io")
    emit(value)
    add("\n")
    return "".join(chunks).encode("utf8")


def _path(raw):
    if type(raw) is not str or not raw.startswith("/") or "\0" in raw or len(raw.encode("utf8")) > 4096:
        raise TransportRefusal("path", "Expected bounded absolute path")
    p = pathlib.Path(raw)
    if str(p) != raw or ".." in p.parts:
        raise TransportRefusal("path", "Path is not canonical", path=raw)
    return p


def _directory(raw):
    p = _path(raw)
    flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC
    fd = os.open("/", flags)
    try:
        for part in p.parts[1:]:
            following = os.open(part, flags, dir_fd=fd)
            os.close(fd)
            fd = following
        return fd
    except BaseException:
        os.close(fd)
        raise


def _read_fd(fd, maximum):
    before = os.fstat(fd)
    if not stat.S_ISREG(before.st_mode) or before.st_size > maximum:
        raise TransportRefusal("regular-file", "Nonregular or excessive input file")
    data = bytearray()
    while True:
        part = os.read(fd, min(65536, maximum + 1 - len(data)))
        if not part:
            break
        data.extend(part)
        if len(data) > maximum:
            raise TransportRefusal("regular-file", "Input grew past bound")
    after = os.fstat(fd)
    if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns) or len(data) != after.st_size:
        raise TransportRefusal("file-drift", "Input file changed while read")
    return bytes(data), after


def _read_leaf(directory_fd, name, maximum=JSON_LIMIT):
    if type(name) is not str or name in ("", ".", "..") or "/" in name:
        raise TransportRefusal("path", "Expected one regular source leaf")
    fd = os.open(name, os.O_RDONLY | os.O_NONBLOCK | os.O_NOFOLLOW | os.O_CLOEXEC, dir_fd=directory_fd)
    try:
        return _read_fd(fd, maximum)[0]
    finally:
        os.close(fd)


def _runtime_file(raw):
    resolved = pathlib.Path(raw).resolve(strict=True)
    parent_fd = _directory(str(resolved.parent))
    try:
        fd = os.open(resolved.name, os.O_RDONLY | os.O_NONBLOCK | os.O_NOFOLLOW | os.O_CLOEXEC, dir_fd=parent_fd)
        try:
            before = os.fstat(fd)
            if not stat.S_ISREG(before.st_mode) or before.st_size > INVENTORY_BYTES:
                raise TransportRefusal("runtime-file", "Invalid runtime file", "runtime", str(resolved))
            digest, total = hashlib.sha256(), 0
            while True:
                chunk = os.read(fd, 65536)
                if not chunk:
                    break
                total += len(chunk)
                if total > INVENTORY_BYTES:
                    raise TransportRefusal("runtime-file", "Runtime file exceeds scan bound", "runtime", str(resolved))
                digest.update(chunk)
            after = os.fstat(fd)
            current = resolved.stat()
            key = lambda s: (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns)
            if key(before) != key(after) or key(after) != key(current) or total != after.st_size:
                raise TransportRefusal("runtime-file-drift", "Runtime file changed while read", "runtime", str(resolved))
            return {"path": str(resolved), "bytes": total, "sha256": digest.hexdigest()}
        finally:
            os.close(fd)
    finally:
        os.close(parent_fd)


def _metadata_discovery():
    """Complete helper/backend metadata imports; no scientific function calls."""
    # Codec lookup/metadata parsing must finish before module/mapping inventories.
    for encoding in ("utf8", "ascii", "latin1"):
        codecs.lookup(encoding)
    distribution = importlib.metadata.distribution("mpmath")
    distribution_files = list(distribution.files or [])
    package = pathlib.Path(mp.__file__).parent
    selected = set(package.rglob("*.py"))
    expanded = set(selected)
    for entry in distribution_files:
        located = pathlib.Path(distribution.locate_file(entry))
        name = str(entry)
        if name.endswith(("LICENSE", "METADATA", "RECORD", "WHEEL", "top_level.txt")):
            selected.add(located)
        if ".dist-info/" in name or name.endswith((".py", ".so")):
            expanded.add(located)
    return selected, expanded


def _proc_bytes(path):
    if path not in ("/proc/self/status", "/proc/self/maps"):
        raise TransportRefusal("runtime-proc", "Unsupported process metadata path", "runtime")
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC)
    try:
        data = bytearray()
        while True:
            chunk = os.read(fd, min(65536, JSON_LIMIT + 1 - len(data)))
            if not chunk:
                break
            data.extend(chunk)
            if len(data) > JSON_LIMIT:
                raise TransportRefusal("runtime-proc", "Process metadata exceeds bound", "runtime", path)
        return bytes(data)
    finally:
        os.close(fd)


def _runtime_python_identity():
    if (sys.executable != PYTHON_PATH or str(pathlib.Path(sys.executable).resolve(strict=True)) != RESOLVED_PYTHON or sys.version != PYTHON_VERSION or platform.python_implementation() != "CPython" or sys.prefix != VENV_PREFIX or sys.base_prefix != BASE_PREFIX or pathlib.Path(sys.base_prefix).resolve(strict=True) != pathlib.Path(RESOLVED_PYTHON).parent.parent):
        raise TransportRefusal("runtime-python", "Reference Python identity differs from original pin", "runtime")


def _collect_imported_modules(expanded):
    module_snapshot = tuple(sorted(sys.modules.items()))
    if not 1 <= len(module_snapshot) <= 4096:
        raise TransportRefusal("runtime-module", "Imported module-count bound exceeded", "runtime")
    modules = []
    for name, module in module_snapshot:
        if len(name.encode("ascii")) > 256 or module is None:
            raise TransportRefusal("runtime-module", "Unsupported imported module identity", "runtime")
        spec = getattr(module, "__spec__", None)
        origin = getattr(spec, "origin", None)
        filename = getattr(module, "__file__", None)
        namespace_paths = []
        if origin in ("built-in", "frozen"):
            kind, recorded_origin = origin, None
        elif origin is None and filename is None and getattr(spec, "submodule_search_locations", None) is not None:
            kind, recorded_origin = "namespace", None
            namespace_paths = sorted({str(pathlib.Path(p).resolve(strict=True)) for p in spec.submodule_search_locations})
            if len(namespace_paths) > 32:
                raise TransportRefusal("runtime-module", "Excessive namespace search paths", "runtime")
        elif filename is not None:
            recorded_origin = str(pathlib.Path(filename).resolve(strict=True))
            kind = "extension" if any(recorded_origin.endswith(suffix) for suffix in importlib.machinery.EXTENSION_SUFFIXES) else "source"
        else:
            raise TransportRefusal("runtime-module", "Unbound module origin: " + name, "runtime")
        if filename is not None:
            expanded.add(pathlib.Path(filename))
        cached = getattr(module, "__cached__", None)
        if cached and pathlib.Path(cached).exists():
            expanded.add(pathlib.Path(cached))
        modules.append({"name": name, "kind": kind, "origin": recorded_origin, "namespace_paths": namespace_paths})
    return modules, module_snapshot


def _collect_mapped_libraries():
    mapping_paths = set()
    map_text = _proc_bytes("/proc/self/maps")
    for line in map_text.decode("utf8", "strict").splitlines():
        parts = line.split(maxsplit=5)
        if len(parts) < 6 or not parts[5].startswith("/"):
            continue
        raw = re.sub(r"\\([0-7]{3})", lambda m: chr(int(m.group(1), 8)), parts[5])
        if re.search(r"\.so(?:\.|$| )", raw):
            if raw.endswith(" (deleted)"):
                raise TransportRefusal("runtime-maps", "Deleted mapped library", "runtime", raw)
            mapping_paths.add(pathlib.Path(raw).resolve(strict=True))
            if len(mapping_paths) > 256:
                raise TransportRefusal("runtime-maps", "Mapped library-count bound exceeded", "runtime")
    return mapping_paths


def runtime_fingerprint():
    selected, expanded = _metadata_discovery()
    _runtime_python_identity()
    if mp.__version__ != "1.3.0" or mp_backend.BACKEND != "python" or mp.mp.dps != 60 or "gmpy2" in sys.modules or not sys.dont_write_bytecode:
        raise TransportRefusal("runtime-backend", "Reference backend/context/bytecode policy differs", "runtime")
    environment = {name: os.environ.get(name) for name in ENV_KEYS}
    if any(value != "1" for value in environment.values()):
        raise TransportRefusal("runtime-environment", "Thread/bytecode environment is not pinned", "runtime")
    # Native metadata/parsing extensions are inventoried; scientific backends refuse.
    if any(name.split(".")[0] in ("numpy", "scipy", "gmpy2", "gmpy", "torch", "jax") for name in sys.modules):
        raise TransportRefusal("runtime-backend", "Unexpected native numerical backend loaded", "runtime")
    status_text = _proc_bytes("/proc/self/status").decode("ascii", "strict")
    thread_lines = [line for line in status_text.splitlines() if line.startswith("Threads:")]
    if len(thread_lines) != 1 or thread_lines[0].split() != ["Threads:", "1"]:
        raise TransportRefusal("runtime-thread", "Actual reference thread count is not one", "runtime")
    configuration = {"arithmetic_backend": "mpmath-libmp-python", "decimal_digits": mp.mp.dps,
                     "binary_precision_bits": mp.mp.prec, "active_native_numeric_backends": [],
                     "gmpy2_loaded": False, "byteorder": sys.byteorder}
    modules, module_snapshot = _collect_imported_modules(expanded)
    mapping_paths = _collect_mapped_libraries()
    expanded.update(mapping_paths)
    executable = pathlib.Path(sys.executable).resolve(strict=True)
    config = pathlib.Path(sys.prefix) / "pyvenv.cfg"
    expanded.update((executable, config, pathlib.Path(__file__)))
    records = {}
    total = 0
    for path in sorted({p.resolve(strict=True) for p in expanded}):
        if len(records) >= INVENTORY_LIMIT:
            raise TransportRefusal("runtime-inventory", "Runtime file-count bound exceeded", "runtime")
        identity = _runtime_file(path)
        total += identity["bytes"]
        if total > INVENTORY_BYTES:
            raise TransportRefusal("runtime-inventory", "Runtime total scan-byte bound exceeded", "runtime")
        records[identity["path"]] = identity
    legacy = {}
    for path in sorted(selected):
        identity = records[str(path.resolve(strict=True))]
        legacy[str(path)] = identity["sha256"]
    legacy_sha = _sha(_canonical(legacy))
    cfg_identity = records[str(config.resolve(strict=True))]
    executable_identity = records[str(executable)]
    init_identity = records[str(pathlib.Path(mp.__file__).resolve(strict=True))]
    if legacy_sha != LEGACY_SHA or cfg_identity["sha256"] != CFG_SHA or executable_identity["sha256"] != PYTHON_SHA or init_identity["sha256"] != MPMATH_INIT_SHA:
        raise TransportRefusal("runtime-pin", "Original runtime artifact/complete package subset changed", "runtime")
    # Metadata collector itself must not change the imported-module set mid-seal.
    if tuple(sorted(sys.modules)) != tuple(name for name, _ in module_snapshot):
        raise TransportRefusal("runtime-module-drift", "Metadata collection imported an unsealed module", "runtime")
    inventory = [records[path] for path in sorted(records)]
    mapped = [records[str(path)] for path in sorted(mapping_paths, key=str)]
    return {"schema_version": 1, "interface_id": "temporal-shared-optical-reference-runtime/v1",
            "python": {"version": sys.version, "implementation": "CPython", "requested_executable": sys.executable,
                       "resolved_executable": str(executable), "executable_bytes": executable_identity["bytes"],
                       "executable_sha256": executable_identity["sha256"], "prefix": sys.prefix, "base_prefix": sys.base_prefix,
                       "pyvenv_cfg": cfg_identity, "dont_write_bytecode": sys.dont_write_bytecode},
            "mpmath": {"version": mp.__version__, "module_origin": str(pathlib.Path(mp.__file__).resolve(strict=True)),
                       "module_root": str(pathlib.Path(mp.__file__).parent.resolve(strict=True)), "backend": mp_backend.BACKEND,
                       "decimal_digits": mp.mp.dps, "binary_precision_bits": mp.mp.prec, "gmpy2_loaded": False,
                       "original_subset_map_sha256": legacy_sha},
            "thread_environment": environment, "actual_thread_count": 1,
            "imported_modules": modules, "mapped_libraries": mapped, "inventory": inventory,
            "inventory_sha256": _sha(_canonical(inventory)), "numeric_configuration": configuration}


def _same_constants(expected, actual, path="request"):
    if type(actual) is not type(expected):
        raise TransportRefusal("request-type", "Exact type mismatch at " + path)
    if type(expected) is dict:
        if set(expected) != set(actual):
            raise TransportRefusal("request-keys", "Closed key mismatch at " + path)
        for key in expected:
            _same_constants(expected[key], actual[key], path + "/" + key)
    elif type(expected) is list:
        if len(expected) != len(actual):
            raise TransportRefusal("request-order", "Array length mismatch at " + path)
        for i, (left, right) in enumerate(zip(expected, actual)):
            _same_constants(left, right, path + "/" + str(i))
    elif expected != actual:
        raise TransportRefusal("request-value", "Frozen constant mismatch at " + path)


def _request_expected(request):
    expected = strict_json(_EXPECTED_REQUEST_TEXT.encode("utf8"))
    expected.update(request_id="temporal-shared-optical-reproducible-request/v1",
                    interface_id="temporal-shared-optical-bounded-controller/v1",
                    status="reviewed-native-synthetic-control", execution=None)
    if type(request) is not dict or type(request.get("source_ports")) is not dict:
        raise TransportRefusal("request-type", "Expected closed request/source_ports object")
    for name in ("reference", "controller", "sdk_verifier"):
        port = request["source_ports"].get(name)
        value = port.get("sha256") if type(port) is dict else None
        if type(value) is not str or SHA_PATTERN.fullmatch(value) is None:
            raise TransportRefusal("source-pin", "Unresolved or invalid implemented source hash: " + name)
        expected["source_ports"][name]["sha256"] = value
    _same_constants(expected, request)


def _source_seal(attempt_fd, packet_fd):
    sources = {name: _sha(_read_leaf(packet_fd, name)) for name in ("consumer.cpp", "reference.py", "controller.py", "verify_sdk.py")}
    sources["contract.snapshot.json"] = _sha(_read_leaf(attempt_fd, "contract.snapshot.json"))
    return {"request_sha256": _sha(_read_leaf(attempt_fd, "request.json")), "current_source_sha256": sources}


def _admit(request_path, output_path):
    request_p, output_p = _path(request_path), _path(output_path)
    if request_p.name != "request.json" or output_p.name != "reference.json" or request_p.parent != output_p.parent:
        raise TransportRefusal("output-path", "Request/output must be approved immediate attempt leaves")
    attempt = request_p.parent
    if attempt.parent.name != "temporal-shared-optical-control" or attempt.parent.parent.name != "results" or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", attempt.name, re.ASCII) is None:
        raise TransportRefusal("attempt-path", "Attempt is outside the fixed packet output layout")
    attempt_fd = _directory(str(attempt))
    packet_fd = None
    try:
        raw = _read_leaf(attempt_fd, "request.json")
        request = strict_json(raw)
        _request_expected(request)
        contract_raw = _read_leaf(attempt_fd, "contract.snapshot.json")
        if len(contract_raw) != 26350 or _sha(contract_raw) != CONTRACT_SHA:
            raise TransportRefusal("contract-pin", "Original contract snapshot differs")
        strict_json(contract_raw)
        packet_path = attempt / "source" / "experiments" / "temporal-shared-optical-control"
        if _path(str(pathlib.Path(__file__).absolute())) != packet_path / "reference.py":
            raise TransportRefusal("source-layout", "Reference is not the approved committed snapshot leaf")
        packet_fd = _directory(str(packet_path))
        seal = _source_seal(attempt_fd, packet_fd)
        if seal["current_source_sha256"]["contract.snapshot.json"] != _sha(contract_raw):
            raise TransportRefusal("contract-drift", "Verified contract changed during admission")
        for port_name in ("consumer", "reference", "controller", "sdk_verifier"):
            port = request["source_ports"][port_name]
            if seal["current_source_sha256"][port["path"]] != port["sha256"]:
                raise TransportRefusal("source-pin", "Source bytes differ from request: " + port_name)
        if seal["request_sha256"] != _sha(raw):
            raise TransportRefusal("request-drift", "Request changed during admission")
        identities = {**seal, "original_contract_sha256": CONTRACT_SHA, "admission_sha256": ADMISSION_SHA,
                      "original_reference_sha256": ORIGINAL_REFERENCE_SHA,
                      "reference_script_sha256": seal["current_source_sha256"]["reference.py"]}
        return request, identities, seal, attempt_fd, packet_fd
    except BaseException:
        os.close(attempt_fd)
        if packet_fd is not None:
            os.close(packet_fd)
        raise


def _empty_partial():
    return {"completed_groups": [], "photon_prefix": [], "photon_error_prefix": [], "frequency_prefix": [],
            "lambda_prefix": [], "detector_prefix": [], "active_sigma_index": None,
            "active_sigma_rows": [], "record_prefix": [], "control_failure": None}


def _observe_partial(stage, local):
    if _PARTIAL is None:
        return
    p = _PARTIAL
    if stage == "photons":
        p["photon_prefix"] = [text(x) for x in local["photons"]]
        p["photon_error_prefix"] = [text(x) for x in local["photon_errors"]]
        p["frequency_prefix"] = list(local["frequency_rows"])
    elif stage == "frequency":
        p["completed_groups"].extend(("frequency", "photons"))
    elif stage == "lambda":
        p["lambda_prefix"] = [text(x) for x in local["lam"]]
        p["completed_groups"].append("lambda")
    elif stage == "sigma-start":
        p["active_sigma_index"] = local["noise"]
        p["active_sigma_rows"] = []
    elif stage == "sigma-row":
        p["active_sigma_rows"] = list(local["sigma_rows"])
    elif stage == "sigma-complete":
        p["detector_prefix"] = list(local["detector_rows"])
        p["completed_groups"].append("detector-sigma0" if local["noise"] == 0 else "detector-sigma075")
        p["active_sigma_index"] = None
        p["active_sigma_rows"] = []
    elif stage == "record":
        p["record_prefix"] = list(local["results"])
        if len(p["record_prefix"]) == 6:
            p["completed_groups"].append("records")


def _clip(value, maximum):
    raw = value.encode("utf8", errors="backslashreplace")
    return raw[:maximum].decode("utf8", errors="ignore"), len(raw) > maximum, _sha(raw)


def _failure_details(candidate):
    """Closed original failure shapes; failed values are diagnostics, not zeros."""
    if type(candidate) is not dict:
        return False
    stage = candidate.get("stage")
    if stage == "positive-mixture-error-gate":
        if set(candidate) != {"stage", "value", "error", "conditional_values", "conditional_errors"}:
            return False
        for key in ("conditional_values", "conditional_errors"):
            if type(candidate[key]) is not list or len(candidate[key]) != 2:
                return False
        values = [("value", candidate["value"]), ("error", candidate["error"])]
        values += [("conditional_value", v) for v in candidate["conditional_values"]]
        values += [("conditional_error", v) for v in candidate["conditional_errors"]]
    elif stage == "conditional-positive-reference-gate":
        if set(candidate) != {"stage", "sigma_index", "channel_index", "lambda", "density", "density_error", "nondetection", "nondetection_error", "coarse_refined_attempt"}:
            return False
        if type(candidate["sigma_index"]) is not int or candidate["sigma_index"] not in (0, 1) or type(candidate["channel_index"]) is not int or not 0 <= candidate["channel_index"] < 8:
            return False
        coarse = candidate["coarse_refined_attempt"]
        row_keys = {"count_probability", "count_probability_error", "below", "below_error"} if candidate["sigma_index"] == 0 else {"density32_per_adu", "density64_per_adu", "density_error_per_adu", "below32", "below64", "below_error", "density_tail_per_adu", "cdf_tail"}
        if type(coarse) is not dict or set(coarse) != row_keys:
            return False
        values = [(key, candidate[key]) for key in ("lambda", "density", "density_error", "nondetection", "nondetection_error")]
        values += list(coarse.items())
    else:
        return False
    for key, value in values:
        if type(value) is not str:
            return False
        number = _finite_decimal(value)
        if ("error" in key or "tail" in key) and number < 0:
            return False
    return True


def _error(error, default_kind="reference"):
    message = str(error)
    clipped, truncated, digest = _clip(message, 4096)
    kind = error.kind if isinstance(error, TransportRefusal) else default_kind
    stage = error.stage if isinstance(error, TransportRefusal) else type(error).__name__
    details = None
    if type(error) is ArithmeticError:
        try:
            candidate = strict_json(message.encode("utf8"))
            if _failure_details(candidate):
                details, stage = candidate, candidate["stage"]
        except (ValueError, UnicodeError, TypeError, RecursionError):
            pass
    if details is None:
        details = {"operation": stage[:256], "path": getattr(error, "path", None), "reason": clipped}
    return {"kind": kind, "stage": stage[:256], "message": clipped,
            "message_truncated": truncated, "message_sha256": digest, "details": details}


def _work(started):
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return {"elapsed_seconds": time.perf_counter() - started,
            "time_frequency_evaluations": evals["time_frequency"], "characteristic_function_evaluations": evals["characteristic_function"],
            "cpu_user_seconds": usage.ru_utime, "cpu_system_seconds": usage.ru_stime, "maximum_rss_kib": usage.ru_maxrss}


def _summary(status, request_sha=None, output_sha=None, output_bytes=None):
    return {"interface_id": "temporal-shared-optical-reference-summary/v1", "status": status,
            "request_sha256": request_sha, "reference_json_sha256": output_sha, "reference_json_bytes": output_bytes}


def _write_all(fd, raw):
    if len(raw) > JSON_LIMIT:
        raise TransportRefusal("output-size", "Reference JSON exceeds bound", "io")
    offset = 0
    while offset < len(raw):
        written = os.write(fd, raw[offset:])
        if written <= 0:
            raise TransportRefusal("output-write", "Reference write made no progress", "io")
        offset += written
    os.fsync(fd)


def _bounded_limits():
    for limit, desired in ((resource.RLIMIT_CPU, 180), (resource.RLIMIT_AS, 1073741824), (resource.RLIMIT_FSIZE, JSON_LIMIT)):
        soft, hard = resource.getrlimit(limit)
        if (soft != resource.RLIM_INFINITY and soft < desired) or (hard != resource.RLIM_INFINITY and hard < desired):
            raise TransportRefusal("resource-policy", "Child inherited tighter unsupported resource bound", "resource")
        resource.setrlimit(limit, (desired, desired))


def _normal(request_path, output_path):
    global _PARTIAL
    started = time.perf_counter()
    request = identities = before = after = payload = problem = None
    attempt_fd = packet_fd = output_fd = None
    baseline = None
    phase = "admission"
    _PARTIAL = _empty_partial()
    limits = ["Synthetic fixed source/shared-optical numerical control only; physical and observational qualification remain blocked.",
              "Shared declared equations, exact SI constants and mpmath GL ancestry are retained; empirical errors are not an original-input certificate."]
    try:
        request, identities, baseline, attempt_fd, packet_fd = _admit(request_path, output_path)
        # Reserve only the new descriptor; never overwrite a receipt or old output.
        phase = "io"
        output_fd = os.open("reference.json", os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW | os.O_CLOEXEC, 0o600, dir_fd=attempt_fd)
        if not stat.S_ISREG(os.fstat(output_fd).st_mode):
            raise TransportRefusal("output-file", "Output descriptor is not regular", "io")
        os.fchmod(output_fd, 0o600)
        phase = "runtime"
        before = runtime_fingerprint()
        phase = "reference"
        _initialize_scientific_constants()
        payload = run()
    except Exception as error:
        problem = _error(error, phase)
        if type(error) is ArithmeticError and type(problem["details"]) is dict and problem["details"].get("stage") in ("positive-mixture-error-gate", "conditional-positive-reference-gate"):
            _PARTIAL["control_failure"] = problem["details"]
    finally:
        try:
            after = runtime_fingerprint()
            if before is not None and before != after:
                raise TransportRefusal("runtime-terminal-drift", "Complete runtime/module/mapping seal changed", "runtime")
        except Exception as error:
            if problem is None:
                problem = _error(error, "runtime")
            else:
                limits.append("Additional terminal refusal: " + _clip(str(error), 2048)[0])
        try:
            if baseline is not None:
                current_directory = _directory(str(_path(request_path).parent))
                try:
                    held, current = os.fstat(attempt_fd), os.fstat(current_directory)
                    if (held.st_dev, held.st_ino) != (current.st_dev, current.st_ino):
                        raise TransportRefusal("directory-terminal-drift", "Admitted attempt directory changed", "admission")
                finally:
                    os.close(current_directory)
                if _source_seal(attempt_fd, packet_fd) != baseline:
                    raise TransportRefusal("source-terminal-drift", "Request/contract/source bytes changed", "admission")
        except Exception as error:
            if problem is None:
                problem = _error(error, "admission")
            else:
                limits.append("Additional source terminal refusal: " + _clip(str(error), 2048)[0])
    status = "completed" if problem is None and payload is not None and before is not None and after is not None else "refused"
    if request is not None:
        axes = dict(request["shape_contract"]["axes"])
        axes.pop("unit_policy", None)
        settings = {"scientific": request["shape_contract"]["unchanged_scientific_settings"], "resources": request["resources"]}
    else:
        # Frozen design metadata is available even when actual input admission fails.
        expected = strict_json(_EXPECTED_REQUEST_TEXT.encode("utf8"))
        axes = dict(expected["shape_contract"]["axes"])
        axes.pop("unit_policy", None)
        settings = {"scientific": expected["shape_contract"]["unchanged_scientific_settings"], "resources": expected["resources"]}
    envelope = {"schema_version": 1, "interface_id": "temporal-shared-optical-reference-transport/v1", "status": status,
                "error": problem, "identities": identities, "axes": axes, "settings": settings,
                "runtime_before": before, "runtime_after": after, "payload": payload,
                "partial": _PARTIAL, "work": _work(started),
                "qualification": {"scope": "empirical-synthetic-fixed-temporal-shared-optical-control/v1", "native_outputs_consumed": False,
                                  "physical_qualification": False, "observational_qualification": False,
                                  "original_input_certificate": False, "limits": limits}}
    output_sha = output_bytes = None
    try:
        if output_fd is not None:
            encoded = _json_bytes(envelope)
            _write_all(output_fd, encoded)
            os.fsync(attempt_fd)
            output_sha, output_bytes = _sha(encoded), len(encoded)
        else:
            print(_clip(json.dumps(problem, sort_keys=True), 8192)[0], file=sys.stderr)
    except Exception as error:
        status = "refused"
        # Preserve already-written bytes. No second write, retry, replacement or unlink.
        print(_clip(str(error), 4096)[0], file=sys.stderr)
    finally:
        for fd in (output_fd, packet_fd, attempt_fd):
            if fd is not None:
                os.close(fd)
    print(_json_bytes(_summary(status, identities["request_sha256"] if identities else None, output_sha, output_bytes)).decode("utf8"), end="")
    return 0 if status == "completed" and output_sha is not None else 1


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime-fingerprint", action="store_true")
    parser.add_argument("--request")
    parser.add_argument("--output")
    args = parser.parse_args(argv)
    if args.runtime_fingerprint:
        if args.request is not None or args.output is not None:
            parser.error("--runtime-fingerprint is exclusive")
    elif args.request is None or args.output is None:
        parser.error("use --runtime-fingerprint or both --request and --output")
    try:
        _bounded_limits()
        if args.runtime_fingerprint:
            print(_json_bytes(runtime_fingerprint()).decode("utf8"), end="")
            return 0
        return _normal(args.request, args.output)
    except Exception as error:
        print(_clip(str(error), 4096)[0], file=sys.stderr)
        print(_json_bytes(_summary("refused")).decode("utf8"), end="")
        return 1


# Literal frozen request design follows. Three implemented source hashes are
# admitted only after checking their exact actual sibling/source bytes. The four
# executable identifiers are exact parent-approved values in _request_expected.
# This template is source code, not an ignored runtime file dependency.
_EXPECTED_REQUEST_TEXT = r'''
{
  "schema_version": 1,
  "request_id": "temporal-shared-optical-reproducible-request/v2-proposal",
  "status": "source-only-proposal-unexecutable-until-parent-freeze",
  "execution": null,
  "packet": "temporal-shared-optical-control",
  "interface_id": "temporal-shared-optical-bounded-controller/v2-proposal",
  "origin": {
    "kind": "development_smoke",
    "component_repository": "Irreducible",
    "component_revision": "bd2dff2d7d84fec45e384326d8e0e9d2aa1e9d09",
    "component_root": "/home/szymon/.codex/worktrees/all19-optical/irreducible",
    "source_directory": "evidence/next-wave-temporal-optical-sdk",
    "source_tracking": "ignored component files; each exact byte identity is independent of tracked HEAD",
    "admission_path": ".work/temporal-optical-admission-20261003/admission.json",
    "admission_sha256": "1ee77abd937be8342f3be946c17453a162f75e9cd315a9a53b05f8fe2db401cc",
    "preserved_directory": ".work/temporal-optical-admission-20261003/preserved",
    "original_sources": [
      {
        "path": "contract-v2.json",
        "bytes": 26350,
        "sha256": "e3a1ceb57d96f599cf82d938462d03feef5d424ad9d2c70f82316d63c34ba6d7",
        "role": "exact current original design copied as contract.snapshot.json"
      },
      {
        "path": "consumer.cpp",
        "bytes": 17548,
        "sha256": "46e17d76d0b03390048fdc66ba97ca7012d8eda215beb9ab2aeea96a92905ab8",
        "role": "unchanged native SDK consumer and external reducer"
      },
      {
        "path": "reference.py",
        "bytes": 12857,
        "sha256": "f912a6f3b0b351b91ece9d26af959200618776908d1842bb64caa63ed87ebe0e",
        "role": "original scientific arithmetic for I/O-only fork review"
      },
      {
        "path": "controller.py",
        "bytes": 16090,
        "sha256": "4f76eca4f9ec4f01ead57978a7031f6dbf28e07eaf14c7e7e286424a1412e535",
        "role": "original169 comparison logic and budgets for semantic preservation"
      },
      {
        "path": "verify_sdk.py",
        "bytes": 4566,
        "sha256": "3354c6fbb7b8b2ec201351b859a3b3e3774cead6db17bb27a015e7bebcd70a68",
        "role": "original full311 SDK fingerprint review ancestry"
      }
    ],
    "historical_result": {
      "attempt": "attempt03",
      "result_sha256": "cc2d51012ee7167205719adcd99d6eea040e00fec696590b4b39c7dcbc26a740",
      "comparison_sha256": "35ddcde6d6cc6b94901db69c5f4d77a283d411a34e218b02e69594fdef3bcc2d",
      "reported_checks": 169,
      "new_execution_evidence": false
    }
  },
  "original_contract": {
    "path": "contract.snapshot.json",
    "bytes": 26350,
    "sha256": "e3a1ceb57d96f599cf82d938462d03feef5d424ad9d2c70f82316d63c34ba6d7",
    "all_numerical_model_endpoint_budgets": "unchanged exact original contract bytes"
  },
  "source_ports": {
    "consumer": {
      "path": "consumer.cpp",
      "bytes": 17548,
      "sha256": "46e17d76d0b03390048fdc66ba97ca7012d8eda215beb9ab2aeea96a92905ab8",
      "input": "no arguments or input JSON; exact source declares synthetic design",
      "output_schema": "external-temporal-optical-sdk-control/v1"
    },
    "reference": {
      "path": "reference.py",
      "sha256": null,
      "interface_id": "temporal-shared-optical-reference-transport/v1",
      "scientific_payload_schema": "original-high-precision-temporal-optical-reference/v1",
      "change_scope": "transport/runtime/failure accounting only; original arithmetic unchanged"
    },
    "controller": {
      "path": "controller.py",
      "sha256": null,
      "arguments": [
        "--attempt",
        "SLUG"
      ],
      "comparison_source_ancestry_sha256": "4f76eca4f9ec4f01ead57978a7031f6dbf28e07eaf14c7e7e286424a1412e535",
      "comparison_semantics": "original169 checks and all budgets unchanged"
    },
    "sdk_verifier": {
      "path": "verify_sdk.py",
      "sha256": null,
      "change_scope": "bounded full inventory/raw CLI/runtime/source admission only"
    },
    "closed_output_contract": {
      "source_owner": "/root/reference_owner",
      "documentary_proposal_sha256": "791ea3981b8eef42eec3a5622e0193cf8ecfe06a5f28472dda032c50f1f1c266",
      "embedded_schema_key": "shape_contract",
      "ignored_proposal_file_required_at_runtime": false
    },
    "pin_rule": "final request pins all new scripts; experiment metadata pins final request bytes; source admission requires clean committed blobs, avoiding cyclic script/request hash literals"
  },
  "sdk": {
    "source_root": "/home/szymon/.codex/worktrees/box-target-integrated/irreducible",
    "artifacts_root": "/home/szymon/.codex/worktrees/box-target-integrated/irreducible",
    "revision": "6f869532c1951ed1afd9f2506b5d05c6cfd03c82",
    "tracked_status": "clean",
    "build_id": "f22c25423cfb9cbac3c2b91a4e514b13ce604e92e7010f55a9aa0bdd42f40f59",
    "manifest": {
      "path": "build/build-manifest-release.json",
      "sha256": "f2e4d6a22257c13f8ab46de51cdfc92654d9857a82a7bc3ea5b8c0edf234f062"
    },
    "archive": {
      "path": "build/native-release/libirred_core.a",
      "sha256": "1cb2b85ad292334f3b5041d669187b04a6fca0916bca8c41070b0ff878499dff"
    },
    "cli": {
      "path": "target/release/irred",
      "sha256": "09b5bfd05ce6c057423b1d41f6f81f8e4d0db44dfd35424fdccd0b34fd9b22d9",
      "discovery": "fixed describe --json; no temporal operation inferred"
    },
    "compiler": {
      "path": "/usr/bin/c++",
      "sha256": "f04191f6a7b2cd7d9a62e1745872b8a6088791e5af6955488c69c9b2c4668bc9"
    },
    "standard_library": {
      "path": "/usr/lib/gcc/x86_64-pc-linux-gnu/16/libstdc++.so",
      "sha256": "f5fc7380f2ae46fa4053a64be04e7b98109f1066a4bbfff3c37042488aa0be0e"
    },
    "source_inventory_count": 311,
    "source_patterns": [
      "src/**/*.rs",
      "tests/**/*.rs",
      "cpp/**/*.cpp",
      "cpp/**/*.h",
      "cpp/**/*.hpp",
      "cpp/**/*.inc",
      "cpp/**/*.cmake",
      "schema/*.json",
      "tools/*.py"
    ],
    "source_explicit": [
      "Cargo.toml",
      "Cargo.lock",
      "build.rs",
      "cpp/CMakeLists.txt"
    ],
    "headers": [
      {
        "path": "cpp/include/irred/temporal_photometry.hpp",
        "sha256": "78df6251066b9272482f79a43acd3f4a48a7128e2594e8d5e4004a7d198ab151"
      },
      {
        "path": "cpp/include/irred/detector_selection.hpp",
        "sha256": "2973496bcf85f51ef42a1f5ab2be319019aa34076d7603b8e2d6eb19ea5a1dbd"
      },
      {
        "path": "cpp/include/irred/photometry_calibration.hpp",
        "sha256": "0a825b6fd7bb659491801dec667e08e3c52b739ea2c36dff1c4a55113a1980fd"
      },
      {
        "path": "cpp/include/irred/photometry.hpp",
        "sha256": "0361c9e43696bbc49791e4b9b1331da6fcae430273cdb2b71dad3a74da64dec4"
      }
    ],
    "guides": [
      {
        "path": "docs/temporal-photometry.md",
        "sha256": "68b93d0278bd148970a9656fb49f84aa040dfadd3fc656c04540c416c8a6ef67"
      },
      {
        "path": "docs/detector-selection.md",
        "sha256": "b02d34e08d0cdc8c96be7f1543011caf90adc9a5841fe42edcc4513226484262"
      },
      {
        "path": "docs/optical-detector.md",
        "sha256": "af5fb3ca260222146ae669decbecb567ff9e63fae6fe2e59ab64479504af9017"
      },
      {
        "path": "docs/photometry-calibration.md",
        "sha256": "98bc0a8ea98ef1e04468b86dde217954bc75997026e30bbcc164267bfbb6125e"
      }
    ],
    "compiler_flags": [
      "-std=c++20",
      "-O3",
      "-DNDEBUG",
      "-Wall",
      "-Wextra",
      "-Wpedantic",
      "-fno-fast-math",
      "-ffp-contract=off"
    ],
    "verification": "complete source/header/build/artifact maps at initial admission, precompile, prenative, prereference and finally; exact original before-each-numerical boundaries; required artifact/source/tool hash checks at consumption boundaries; actual compiled consumer identity before/after"
  },
  "reference_runtime": {
    "requested_executable": "/home/szymon/Projects/irreducible/evidence/project-review/science/recovery-abundance-growth-20261002/reference-venv/bin/python",
    "resolved_executable": "/home/szymon/.local/share/mise/installs/python/3.14.8/bin/python3.14",
    "executable_sha256": "815b1275bf87e7595fe1cdf3f2efc6ae58e27b842262a7b530de8610ac2760ba",
    "implementation": "CPython",
    "version": "3.14.8 (main, Oct  1 2026, 20:55:04) [GCC 16.2.1 20260810]",
    "prefix": "/home/szymon/Projects/irreducible/evidence/project-review/science/recovery-abundance-growth-20261002/reference-venv",
    "base_prefix": "/home/szymon/.local/share/mise/installs/python/latest",
    "pyvenv_cfg_sha256": "11a0241468a7d4cb44419cd7c7a87cd6149c3e90e38b75be6e1a201279621e43",
    "mpmath_version": "1.3.0",
    "mpmath_backend": "python",
    "mpmath_init_sha256": "b241584d2c1fc0304b0a1015ea923749d7b0800411dd406dcab7c82bf25d9fe8",
    "gmpy2_loaded": false,
    "original_seal_sha256": "ae4e395559ac113c6d7ecbbe165128509d350e4b5fa5d9656251e19b0f096745",
    "legacy_package_inventory_sha256": "982dc09f2ddf1206cc9cd4467c7be1ce4238a612489b43bf5f5d196a8b2495a8",
    "legacy_inventory_hash_encoding": "original seal reference_runtime.source_inventory absolute-path-to-sha map; sorted keys compact JSON UTF-8, no trailing LF",
    "expanded_inventory_rule": "complete metadata/scientific-code module imports before first actual map, no scientific calls; independently equal complete imported_modules and mapped_libraries sets plus inventory/config before/after; no numerical prewarm or installation"
  },
  "native_dependencies_proposed": {
    "enabled": true,
    "scope": "ELF loader resolved dependencies under scrubbed environment, not actual native process mapped inventory",
    "loader_path": "/usr/lib/ld-linux-x86-64.so.2",
    "loader_sha256": "f5e11cc62f8f2c24dff532982559f4303dba3b17e7450551372a88c8e7ea8757",
    "boundaries": [
      "before native",
      "terminal finally"
    ],
    "freeze_decision": "parent approved in principle; final literal invocation/environment/schema reviewed before executable request freeze"
  },
  "resources": {
    "jobs": 1,
    "threads": 1,
    "compile_cpu_seconds": 180,
    "compile_wall_seconds": 180,
    "native_cpu_seconds": 180,
    "native_wall_seconds": 180,
    "reference_cpu_seconds": 180,
    "reference_wall_seconds": 180,
    "structural_child_cpu_seconds": 180,
    "structural_child_wall_seconds": 180,
    "git_metadata_cpu_seconds": 30,
    "git_metadata_wall_seconds": 30,
    "git_metadata_regular_file_limit_bytes": 1048576,
    "git_archive_regular_file_limit_bytes": 4194304,
    "git_archive_stdout_bytes": 1048576,
    "git_archive_independent_output_bytes": 4194304,
    "address_space_bytes_each_child": 1073741824,
    "stdout_bytes_each_child": 1048576,
    "stderr_bytes_each_child": 1048576,
    "reference_json_bytes": 1048576,
    "compiled_executable_bytes": 67108864,
    "compile_regular_file_limit_bytes": 67108864,
    "native_reference_regular_file_limit_bytes": 1048576,
    "source_tree_bytes": 2097152,
    "source_archive_bytes": 4194304,
    "source_file_bytes": 1048576,
    "attempt_store_bytes": 268435456,
    "receipt_file_bytes": 1048576,
    "terminal_fallback_bytes": 65536,
    "terminal_fallback_artifact_count": 128,
    "terminal_artifact_relative_path_bytes": 256,
    "terminal_fallback_message_bytes": 2048,
    "runtime_inventory_files": 4096,
    "runtime_inventory_total_bytes_scanned": 67108864,
    "maximum_children": 28,
    "child_count_scope": "DIRECT controller launches; compiler descendants are one serial build job under process-group wall/store kill; RLIMIT_CPU is per-process, not aggregate descendant CPU",
    "git_metadata_children": 19,
    "non_git_children": 9,
    "maximum_aggregate_child_wall_seconds": 2190,
    "aggregate_wall_scope": "sum9*180+19*30 direct-child watchdog budgets, excluding bounded in-process hashing/serialization; not a whole-controller wall assertion",
    "cpu_measurement_scope": "actual wait4 rusage of reaped child as attributed by OS; do not claim an enforced aggregate compiler-tree CPU180 limit",
    "retry_count": 0,
    "schedule": [
      "repro-initial-head",
      "repro-initial-branch",
      "repro-initial-status",
      "repro-initial-tree",
      "repro-source-archive",
      "sdk-initial-head",
      "sdk-initial-status",
      "compiler-version",
      "discovery",
      "runtime-before",
      "sdk-precompile-head",
      "sdk-precompile-status",
      "compile",
      "loader-before",
      "sdk-prenative-head",
      "sdk-prenative-status",
      "native",
      "sdk-prereference-head",
      "sdk-prereference-status",
      "reference",
      "loader-after-finally",
      "runtime-after-finally",
      "sdk-final-head",
      "sdk-final-status",
      "repro-final-head",
      "repro-final-branch",
      "repro-final-status",
      "repro-final-tree"
    ],
    "git_fixed_global_options": [
      "-c",
      "core.fsmonitor=false",
      "-c",
      "core.preloadIndex=false",
      "-c",
      "index.threads=1",
      "--no-pager"
    ],
    "git_fixed_environment": {
      "GIT_OPTIONAL_LOCKS": "0",
      "GIT_CONFIG_NOSYSTEM": "1"
    },
    "metadata_ledger_rule": "literal controller operations described in PROPOSAL.md, never execute commands from JSON; one archive and no per-file git show/other helper launches",
    "environment": {
      "MKL_NUM_THREADS": "1",
      "NUMEXPR_NUM_THREADS": "1",
      "OMP_NUM_THREADS": "1",
      "OPENBLAS_NUM_THREADS": "1",
      "BLIS_NUM_THREADS": "1",
      "VECLIB_MAXIMUM_THREADS": "1",
      "PYTHONDONTWRITEBYTECODE": "1"
    },
    "environment_policy": "fixed PATH/tool paths, scrub LD_PRELOAD/LD_LIBRARY_PATH/PYTHONPATH/PYTHONHOME and unrelated numerical overrides; private compile TMPDIR inside attempt; actual effective environment recorded",
    "native_payload_policy_bytes": 1048576,
    "native_and_reference_math_policy": "unchanged original exact contract"
  },
  "transport": {
    "output_root": "results/temporal-shared-optical-control",
    "attempt_slug_pattern": "[A-Za-z0-9][A-Za-z0-9_-]{0,63}",
    "fresh_only": true,
    "reference_output_name": "reference.json",
    "source_layout": {
      "committed_source_directory": "source",
      "packet_directory": "source/experiments/temporal-shared-optical-control",
      "immediate_request": "request.json",
      "immediate_contract": "contract.snapshot.json",
      "reference_script": "source/experiments/temporal-shared-optical-control/reference.py",
      "contract_relative_to": "request parent",
      "source_ports_relative_to": "approved source/packet directory"
    },
    "maximum_json_depth": 32,
    "maximum_json_nodes": 65536,
    "maximum_string_bytes": 8192,
    "maximum_numeric_token_bytes": 256,
    "maximum_decimal_adjusted_exponent_absolute": 1000,
    "native_float_parser": "finite Decimal from exact emitted JSON number tokens; no binary64 intermediate",
    "reference_scalar_parser": "strict finite decimal strings retaining original60 digits",
    "decimal_string_pattern": "-?(0|[1-9][0-9]*)(\\.[0-9]+)?([eE][+-]?[0-9]+)?",
    "decimal_regex_encoding_rule": "JSON source contains two backslashes before dot, decoded regex exactly one backslash then dot (UTF8hex5c2e); decoded grammar is authoritative",
    "duplicate_keys": "reject",
    "nonfinite_constants": "reject",
    "source_admission": "clean committed repository HEAD and exact current/snapshot source blobs",
    "capture": "read raw stdout/stderr/reference JSON once bounded, hash and parse same bytes; preserve all outcomes; truncated overflow prefixes never impersonate complete streams",
    "terminal": "always rehash captured streams/files and all source/SDK/runtime/input/executable identities; raw drift revokes acceptance without replacing original parsed objects/checks",
    "terminal_size_strategy": "sealed artifacts retain full runtime/source/SDK maps, scientific objects and169checks; compact terminal binds artifact identities; oversize/serialization refuses via distinct closed64KiBfallback before creating record once",
    "earned_payload_rule": "retain a completed scientific payload if only later identity/transport fails; statusrefused and comparisonacceptance unavailable; null only never-earned fullpayload with retained partialprefixes",
    "seal": "exclusive new terminal record; regular output files0444/executable0555, no prior receipt edits"
  },
  "qualification": {
    "execution": "unperformed",
    "numerical": "new committed attempt and original169 comparisons required",
    "inference": "blocked",
    "interpretation": "synthetic-conditional-only",
    "physical_qualification": false,
    "observational_qualification": false,
    "unmet_gates": "all original contract unmet_gates retained; no measured template/calibration/rawdetector/population/covariance qualification"
  },
  "freeze_sequence": {
    "before_implementation": [
      "parent accepts embedded closed output/runtime/interface/model/domain/allocation contract",
      "parent accepts layout/path/resource/direct-child ledger and loader policy"
    ],
    "after_code_before_execution": [
      "I/O-only reference fork and dedicated controller/verifier source completed and deliberately reviewed",
      "all null port hashes replaced",
      "final executable request hash bound by experiment metadata",
      "clean committed checkpoint",
      "parent explicit single-job grant before tests/compile/science"
    ]
  },
  "shape_contract": {
    "schema_version": 1,
    "id": "temporal-shared-optical-closed-shape-proposal/v2",
    "status": "source-only proposal; requires parent contract/request approval before implementation",
    "scope": "transport/type/domain/order proposal only; original scientific arithmetic and budgets unchanged",
    "provenance": {
      "admission_sha256": "1ee77abd937be8342f3be946c17453a162f75e9cd315a9a53b05f8fe2db401cc",
      "original_contract_sha256": "e3a1ceb57d96f599cf82d938462d03feef5d424ad9d2c70f82316d63c34ba6d7",
      "original_consumer_sha256": "46e17d76d0b03390048fdc66ba97ca7012d8eda215beb9ab2aeea96a92905ab8",
      "original_reference_sha256": "f912a6f3b0b351b91ece9d26af959200618776908d1842bb64caa63ed87ebe0e",
      "original_controller_sha256": "4f76eca4f9ec4f01ead57978a7031f6dbf28e07eaf14c7e7e286424a1412e535",
      "original_runtime_subset_map_sha256": "982dc09f2ddf1206cc9cd4467c7be1ce4238a612489b43bf5f5d196a8b2495a8",
      "derivation": "Actual preserved consumer emitters and reference return/exception definitions; original sealed records used only for shape/identity/failure semantics. No computed acceptance values reused as constants or expected answers.",
      "previous_proposal_sha256": "30f3ba5add7849a1ebc2d6d223a0ddbc001258ad0b8bdd8be645a40e39a9df36"
    },
    "notation": {
      "object_rule": "Every named object below is closed; every listed field is required unless explicitly optional. Nullable means exactly JSON null or the named type. Unknown, duplicate or missing keys refuse transport admission.",
      "array_rule": "array(T,min,max) is ordered, bounded and element-typed. Complete routes require exact stated axis order; partial arrays are prefixes, never silently compressed or renumbered.",
      "bool": "type(value) is bool; never accept integers0/1",
      "uint": "type(value) is int and not bool, nonnegative, bounded as stated; integral Decimal or floating tokens do not become integer identity/work fields",
      "number": "finite JSON integer or decimal number token, decoded directly to int/Decimal; bool and decoded float objects refused. Preserve exact emitted decimal digits, including21-digit native longdouble serialization; no binary64 cast for comparisons.",
      "nonnegative": "number>=0, finite",
      "probability": "number in[0,1], finite",
      "decimal_string": "ASCII finite decimal string length<=256; grammar -?(0|[1-9][0-9]*)(\\.[0-9]+)?([eE][+-]?[0-9]+)?; adjusted exponent absolute value<=1000; parse directly as Decimal without float intermediary",
      "decimal_error": "decimal_string>=0; positivity/range relations checked by context",
      "sha256": "exact64 lowercase hexadecimal ASCII characters",
      "path": "absolute path string, UTF8<=4096bytes, no NUL; input/output/source admission refuses symlinks/escaping paths. Pinned readonly runtime requested aliases may be symlinks only where explicitly identified as requested paths; record and check their exact resolved target. Inventory filepaths are resolved regular files.",
      "status": [
        "ok",
        "invalid_input",
        "nonfinite_input",
        "overflow",
        "work_limit",
        "outside_domain",
        "singular",
        "not_positive_definite",
        "conditioning_budget_exceeded"
      ],
      "unknown_status": "The original emitter fallback 'unknown' is a transport refusal; raw bytes and failure diagnostics are retained.",
      "bounds": {
        "json_bytes": 1048576,
        "stderr_bytes": 1048576,
        "stdout_bytes": 1048576,
        "maximum_depth": 32,
        "maximum_total_items": 65536,
        "maximum_string_utf8_bytes": 8192,
        "maximum_numeric_token_characters": 256,
        "identifier_ascii_bytes": 256,
        "message_utf8_bytes": 4096
      }
    },
    "axes": {
      "synthetic_identity": "synthetic-fixed-bilinear-source-two-exposures-two-bands-one-common-two-state-optical-law",
      "grid": {
        "index": 0,
        "id": "synthetic-bilinear-ux16-v1",
        "source_origin": "original analytic synthetic control",
        "requested_mask": 4,
        "requested_group": "transmitted_photons",
        "omitted_groups": [
          "mean_flux",
          "energy"
        ]
      },
      "states": [
        {
          "index": 0,
          "id": "S0",
          "relative_mass": 1,
          "normalized_mass": "1/4"
        },
        {
          "index": 1,
          "id": "S1",
          "relative_mass": 3,
          "normalized_mass": "3/4"
        }
      ],
      "physical_bands": [
        {
          "index": 0,
          "id": "blue",
          "observed_wavelength_knots_metre": [
            2,
            3
          ]
        },
        {
          "index": 1,
          "id": "red",
          "observed_wavelength_knots_metre": [
            3,
            4
          ]
        }
      ],
      "prepared_band_order": [
        "S0/blue",
        "S0/red",
        "S1/blue",
        "S1/red"
      ],
      "exposures": [
        {
          "index": 0,
          "id": "E0",
          "observer_interval_second": [
            9,
            11
          ],
          "full_duration_second": 2,
          "covered_duration_second": 2,
          "expected_coverage": "full",
          "covered_fraction": "1"
        },
        {
          "index": 1,
          "id": "E1",
          "observer_interval_second": [
            12,
            15
          ],
          "full_duration_second": 3,
          "covered_duration_second": 2,
          "expected_coverage": "partial",
          "covered_fraction": "2/3"
        }
      ],
      "channels": [
        "E0/blue",
        "E0/red",
        "E1/blue",
        "E1/red"
      ],
      "primary_rows": [
        "S0/E0/blue",
        "S0/E0/red",
        "S0/E1/blue",
        "S0/E1/red",
        "S1/E0/blue",
        "S1/E0/red",
        "S1/E1/blue",
        "S1/E1/red"
      ],
      "row_mapping": "i=4*state+2*exposure+band; native band_index=2*state+band; detector state_index=i//4 and channel_index=i%4",
      "sigma_controls_electrons": [
        0,
        0.75
      ],
      "reference_record_order": [
        "discrete-all-detected",
        "discrete-selected-all-detected",
        "discrete-censored",
        "continuous-all-detected",
        "continuous-selected-all-detected",
        "continuous-censored"
      ],
      "per_sigma_record_indices": [
        0,
        1,
        2
      ],
      "detector_batch_slots": [
        "nominal",
        "lower-photon-sensitivity",
        "upper-photon-sensitivity"
      ],
      "detector_row_slots": [
        "joint-all-detected",
        "selected-control-joint-numerator",
        "joint-censored",
        "threshold-nondetection-probe0",
        "threshold-nondetection-probe1",
        "threshold-nondetection-probe2"
      ],
      "temporal_control_ids": [
        "no-time-overlap",
        "zero-duration",
        "invalid-redshift",
        "required-state-row-refusal",
        "temporal-work-refusal",
        "temporal-row-admission-refusal"
      ],
      "detector_control_ids": [
        "detector-QE-domain",
        "detector-work-refusal",
        "selected-nondetection",
        "zero-noise-threshold-boundary"
      ],
      "reducer_control_ids": [
        "required-QE-refusal",
        "required-work-refusal",
        "selected-censored-refusal"
      ],
      "record_measure": {
        "sigma0_record0": "joint discrete counting probability",
        "sigma0_record1": "selected-only discrete counting probability given same ALL-four-detected event",
        "sigma0_record2": "joint three detected counts times one nondetection probability",
        "sigma075_record0": "joint numerical ADU^-4 density",
        "sigma075_record1": "selected-only numerical ADU^-4 density given same ALL-four-detected event",
        "sigma075_record2": "joint numerical ADU^-3 density times one nondetection probability"
      },
      "selection_event": "ALL four channels detected, Y>=1.5ADU; zero-noise K>=1; one shared-state mixture denominator",
      "unit_policy": "photon expectation=count, lambda=expected electrons, durations=observer seconds, continuous density=caller ADU coordinate; reference log/difference measures follow record map. Native payload lacks these labels: immutable request/controller record binds them without modifying consumer bytes.",
      "alternative_records": "Never multiply record controls together as independent measurements; conditionally independent readouts share one optical state across both exposures.",
      "preservation": "State IDs/masses, rows and requested/omitted masks retained on failures; no dropping, reordering or renormalizing."
    },
    "native_types": {
      "Outcome": {
        "availability": "uint enum0=omitted,1=available,2=failed",
        "status": "status",
        "value": "nullable nonnegative number"
      },
      "TemporalRow": {
        "index": "uint[0,7]",
        "grid_index": "uint[0,0]",
        "band_index": "uint[0,4];4 only declared required-state-row-refusal",
        "source_epoch_second": "number exactly10",
        "observer_lower_second": "number",
        "observer_upper_second": "number",
        "admission_status": "status",
        "coverage": "enum unassessed,no_overlap,partial,full",
        "observer_duration_second": "nonnegative",
        "covered_observer_second": "nonnegative",
        "covered_fraction": "probability",
        "segment_work": "uint[0,512]",
        "photons": "Outcome",
        "mean_flux": "Outcome",
        "energy": "Outcome"
      },
      "TemporalBatch": {
        "status": "status",
        "segment_work": "uint[0,512]",
        "required_rows_admitted": "bool",
        "rows": "array(TemporalRow,0,8)"
      },
      "DetectorRow": {
        "status": "status",
        "detected": "bool",
        "threshold_adu": "number exactly1.5",
        "electron_count": "nullable uint[0,4294967295]",
        "measured_adu": "nullable number",
        "zero_probability": "bool",
        "log_value": "nullable number",
        "detection_probability": "nullable probability",
        "log_error": "nonnegative"
      },
      "DetectorBatch": {
        "status": "status",
        "photons": "nonnegative",
        "QE": "number; allow declared refusedQE1.01, never impose[0,1] on a retained failed-control source",
        "full_exposure_second": "number exactly2 or3 by channel/control",
        "sigma_electrons": "number exactly0 or0.75 by sigma/control",
        "poisson_terms": "uint[0,256]",
        "omitted_tail": "probability",
        "rows": "array(DetectorRow,0,6); primary successfulbatch6, single-detector controls0..1"
      },
      "DetectorAttempt": {
        "state_index": "uint[0,1]",
        "channel_index": "uint[0,3]",
        "lambda_electrons": "nullable nonnegative",
        "batches": "array(nullable DetectorBatch,3,3)"
      },
      "ConditionalState": {
        "state_id": "enum S0,S1 by slot",
        "mass": "number exactly0.25 or0.75 by slot",
        "status": "status",
        "log_record": "number",
        "log_event": "number",
        "record_error": "nonnegative",
        "event_error": "nonnegative"
      },
      "ReducerRecord": {
        "record_index": "uint[0,2]",
        "conditional": "array(ConditionalState,2,2)",
        "selected_only": "bool",
        "status": "status",
        "log_joint": "nullable number",
        "log_event": "nullable number",
        "log_value": "nullable number",
        "joint_log_error": "nonnegative",
        "event_log_error": "nonnegative",
        "value_log_error": "nonnegative",
        "product_per_exposure_log_record": "nullable number",
        "product_per_exposure_log_event": "nullable number",
        "record_relative_difference": "nullable number",
        "event_relative_difference": "nullable number"
      },
      "DetectorControl": {
        "sigma_electrons": "number exactly0 or0.75 by slot",
        "attempts": "array(DetectorAttempt,8,8)",
        "records": "array(ReducerRecord,3,3)"
      },
      "TemporalControl": {
        "id": "exacttemporalcontrolID by slot",
        "joint_withheld": "bool",
        "native": "TemporalBatch"
      },
      "SingleDetectorControl": {
        "id": "exactdetectorcontrolID by slot",
        "native": "DetectorBatch"
      },
      "ExternalReducerControl": {
        "id": "exactreducercontrolID by slot",
        "retained_state_masses": "array(number,2,2) exactly[0.25,0.75]",
        "replacement_batch": "nullable DetectorBatch",
        "reducer": "ReducerRecord"
      },
      "NativePayload": {
        "schema": "exact external-temporal-optical-sdk-control/v1",
        "model_id": "exact finite_bilinear_rest_spectral_time_zero_outside_full_observer_exposure_mean",
        "detector_model_id": "exact DETECTOR/Poisson-arrivals-fixed-QE-Gaussian-read/v1",
        "selection_id": "exact DETECTOR/threshold-joint-or-selected-censoring/v1",
        "sampled_empirical_sensitivity": "number; source-bound longdouble literal3e-12, exact emitted token preserved",
        "temporal_photon_endpoint_sensitivity": "number; source-bound longdouble literal3.2e-12, exact emitted token preserved",
        "prepared_status": "status",
        "grid_count": "uint[0,1]",
        "band_count": "uint[0,4]",
        "retained_temporal_payload_bytes": "nullable uint[0,1048576]",
        "prepare_seconds": "nonnegative",
        "temporal_evaluate_seconds": "nonnegative",
        "primary": "TemporalBatch",
        "detector_controls": "array(DetectorControl,0,2)",
        "primary_poisson_terms": "uint[0,12288]",
        "actual_temporal_controls": "array(TemporalControl,6,6)",
        "actual_detector_controls": "array(SingleDetectorControl,4,4)",
        "actual_external_reducer_controls": "array(ExternalReducerControl,0,3)",
        "all_poisson_terms": "uint[0,13824]",
        "all_temporal_segment_work": "uint[0,3584]",
        "detector_six_row_payload_bound_bytes": "nullable uint[0,1048576]",
        "total_native_seconds": "nonnegative"
      }
    },
    "native_semantic_rules": [
      "Primary required_rows_admitted exactly recomputes existing emitter rule: batch statusok, eight rows, each admissionok and photons availability1/statusok/finite nonnegative. This flag alone never establishes numerical acceptance.",
      "Outcome available1 requires statusok and finite nonnegative value; omitted0 and failed2 require null value. Requested mean_flux and energy remain omitted0/statusok/value=null on every row because mask4 requests photons only.",
      "Row identity and exposure bounds are exact from source/request; failed temporal rows may retain unassessed coverage and zero defaults. Do not apply successfulduration/coverage relations to withheld failed values.",
      "Native detector_controls has exactly2 groups and external_reducer_controls exactly3 iff primary is admitted; otherwise both are empty. Temporal6 and single-detector4 control families are emitted even when primary fails.",
      "Detector batch slot0 is mandatory; endpoint slots1/2 null only for exact nominal photons0. Positive nominal expectation requires both complete endpoint attempts, including their actual outside_domain/refusal records. Never infer skip from a missing/failed row.",
      "Nominal/endpoint six-row source: rows0/1 detected; row2 nondetected only channel3; rows3/4/5 threshold-only nondetection. sigma0 detected rows hold electron_count and nullmeasured_adu; sigma075 detected rows hold measured_adu and nullcount; nondetection holds bothnull.",
      "Detector failed status may coexist with zero rows or failed individual rows under batchstatusok. Structuralzero is statusok/zero_probabilitytrue/log_valuenull; it is not a numerical refusal. Failedstatus finite-log availability follows actual SDK; required reducers cannot consume a failedrow.",
      "Reducer statusok requires all seven aggregate/comparator scalar values finite; refusedstatus requires all sevennull. Conditional state logs/errors remain emitted even on refusal; status binds their validity and zero placeholders never become earned complete conditional data.",
      "Selected record1 divides its joint numerator by the same ALL-four-detected event exactly once. Selected censored record2 refusesinvalid_input and withholds aggregate/comparator groups. Replacement batch isnull only selected-censored-refusal; requiredQE/work controls retain the actual failedreplacement batch.",
      "All work sums recompute over actually retained nonnull batches/controls, without using historical measured work counts as expected constants. EachbatchPoissonterms<=256; primary48batches bound12288, all54batches bound13824; seven temporal calls each<=512 yield3584 totalcap.",
      "Exactly representable source values such as0.25/0.75/1.5 and integerduration retain exact numeric equality. QE source literal1.01 is binary64 and its full emitted decimal token is preserved; admit it only in declared failed-control contexts instead of asserting exactdecimal1.01. Longdouble sensitivity literal source identities are bound by exact unchanged compiled source and SDK headers, never rounded to a binary64 parser target.",
      "The unchanged C++ emitter can print a non-JSON inf/nan for an unexpected nonfinite diagnostic. Strict parsing must refuse such output and preserve complete bounded raw bytes/execution failure; do not sanitize it into a successful typed payload or alter the original consumer."
    ],
    "reference_types": {
      "FrequencyRow": {
        "polynomial": "decimal_string",
        "polynomial_error": "decimal_error",
        "frequency32": "decimal_string",
        "frequency64": "decimal_string",
        "frequency32_error": "decimal_error",
        "frequency64_error": "decimal_error"
      },
      "DiscreteDetectorReferenceRow": {
        "count_probability": "decimal_string",
        "count_probability_error": "decimal_error",
        "below": "decimal_string",
        "below_error": "decimal_error"
      },
      "ContinuousDetectorReferenceRow": {
        "density32_per_adu": "decimal_string",
        "density64_per_adu": "decimal_string",
        "density_error_per_adu": "decimal_error",
        "below32": "decimal_string",
        "below64": "decimal_string",
        "below_error": "decimal_error",
        "density_tail_per_adu": "decimal_error",
        "cdf_tail": "decimal_error"
      },
      "DetectorReferenceGroup": {
        "sigma_electrons": "number exactly0 or0.75 by slot",
        "rows": "array(DiscreteDetectorReferenceRow orContinuousDetectorReferenceRow accordingtosigma,8,8)"
      },
      "ReferenceRecord": {
        "sigma_index": "uint[0,1]",
        "record_index": "uint[0,2]",
        "conditional_record": "array(decimal_string,2,2)",
        "conditional_event": "array(decimal_string,2,2)",
        "log_joint": "decimal_string",
        "joint_log_error": "decimal_error",
        "log_event": "decimal_string",
        "event_log_error": "decimal_error",
        "log_value": "decimal_string",
        "value_log_error": "decimal_error",
        "product_per_exposure_log_record": "decimal_string",
        "product_per_exposure_log_event": "decimal_string",
        "product_per_exposure_record_log_error": "decimal_error",
        "product_per_exposure_event_log_error": "decimal_error",
        "record_relative_difference": "decimal_string",
        "event_relative_difference": "decimal_string",
        "record_signed_difference": "decimal_string",
        "record_signed_difference_error": "decimal_error",
        "event_signed_difference": "decimal_string",
        "event_signed_difference_error": "decimal_error",
        "record_relative_difference_error": "decimal_error",
        "event_relative_difference_error": "decimal_error",
        "positive_joint": "decimal_string",
        "positive_joint_error": "decimal_error",
        "positive_event": "decimal_string",
        "positive_event_error": "decimal_error"
      },
      "ReferencePayload": {
        "schema": "exact original-high-precision-temporal-optical-reference/v1",
        "mpmath_version": "exact1.3.0",
        "decimal_digits": "uint exactly60",
        "shared_ancestry": "array(ASCIIidentifier,2,2) exactly['declared equations and exact SI constants','mpmath arithmetic/GL node generator; no native kernels']",
        "frequency_seconds": "nonnegative",
        "elapsed_seconds": "nonnegative",
        "integrand_evaluations": "closed{time_frequency:uint[0,30720],characteristic_function:uint[0,73728]}; complete route exactly30720/73728 from original loops",
        "photons": "array(decimal_string,8,8)",
        "photon_errors": "array(decimal_error,8,8)",
        "frequency_rows": "array(FrequencyRow,8,8)",
        "lambda_electrons": "array(decimal_string,8,8)",
        "detector_rows": "array(DetectorReferenceGroup,2,2)",
        "records": "array(ReferenceRecord,6,6)"
      },
      "PositiveMixtureFailure": {
        "stage": "exact positive-mixture-error-gate",
        "value": "decimal_string",
        "error": "decimal_error",
        "conditional_values": "array(decimal_string,2,2)",
        "conditional_errors": "array(decimal_error,2,2)"
      },
      "ConditionalPositiveFailure": {
        "stage": "exact conditional-positive-reference-gate",
        "sigma_index": "uint[0,1]",
        "channel_index": "uint[0,7]",
        "lambda": "decimal_string",
        "density": "decimal_string",
        "density_error": "decimal_error",
        "nondetection": "decimal_string",
        "nondetection_error": "decimal_error",
        "coarse_refined_attempt": "DiscreteDetectorReferenceRow orContinuousDetectorReferenceRow accordingtosigma"
      },
      "TransportFailureDetails": {
        "operation": "boundedASCIIidentifier",
        "path": "nullable path",
        "reason": "boundedUTF8message"
      },
      "ReferenceError": {
        "kind": "enum admission,runtime,reference,io,resource",
        "stage": "boundedASCIIidentifier",
        "message": "boundedUTF8message",
        "message_truncated": "bool",
        "message_sha256": "sha256 of originaluntruncated UTF8message",
        "details": "nullable PositiveMixtureFailure orConditionalPositiveFailure orTransportFailureDetails"
      },
      "ReferencePartial": {
        "completed_groups": "array(enum frequency,photons,lambda,detector-sigma0,detector-sigma075,records,0,6); original group completion order",
        "photon_prefix": "array(decimal_string,0,8)",
        "photon_error_prefix": "array(decimal_error,0,8)",
        "frequency_prefix": "array(FrequencyRow,0,8)",
        "lambda_prefix": "array(decimal_string,0,8)",
        "detector_prefix": "array(DetectorReferenceGroup,0,2)",
        "active_sigma_index": "nullable uint[0,1]",
        "active_sigma_rows": "array(DiscreteDetectorReferenceRow orContinuousDetectorReferenceRow accordingtoactive_sigma_index,0,8)",
        "record_prefix": "array(ReferenceRecord,0,6)",
        "control_failure": "nullable PositiveMixtureFailure orConditionalPositiveFailure"
      },
      "ReferenceIdentities": {
        "request_sha256": "sha256",
        "original_contract_sha256": "sha256 exact originalv2",
        "admission_sha256": "sha256 exactpreserved admission",
        "original_reference_sha256": "sha256 exactf912 ancestry",
        "reference_script_sha256": "sha256 actualnewforksource",
        "current_source_sha256": "closed map of consumer.cpp,reference.py,controller.py,verify_sdk.py,contract.snapshot.json tosha256; actual admitted snapshot bytes"
      },
      "ReferenceWork": {
        "elapsed_seconds": "nonnegative",
        "time_frequency_evaluations": "uint[0,30720]",
        "characteristic_function_evaluations": "uint[0,73728]",
        "cpu_user_seconds": "nonnegative",
        "cpu_system_seconds": "nonnegative",
        "maximum_rss_kib": "uint[0,4294967295]"
      },
      "ReferenceQualification": {
        "scope": "exact empirical-synthetic-fixed-temporal-shared-optical-control/v1",
        "native_outputs_consumed": "bool exactlyfalse",
        "physical_qualification": "bool exactlyfalse",
        "observational_qualification": "bool exactlyfalse",
        "original_input_certificate": "bool exactlyfalse",
        "limits": "array(boundedUTF8message,1,16)"
      },
      "ReferenceEnvelope": {
        "schema_version": "uint exactly1",
        "interface_id": "exact temporal-shared-optical-reference-transport/v1",
        "status": "enum completed,refused",
        "error": "nullable ReferenceError",
        "identities": "nullable ReferenceIdentities untiladmissioncompleted",
        "axes": "exact immutable axes snapshot excluding documentaryunitnotes",
        "settings": "exact approvedrequest scientific/resource settings snapshot",
        "runtime_before": "nullable RuntimeFingerprint",
        "runtime_after": "nullable RuntimeFingerprint",
        "payload": "nullable ReferencePayload",
        "partial": "ReferencePartial",
        "work": "ReferenceWork",
        "qualification": "ReferenceQualification"
      }
    },
    "reference_semantic_rules": [
      "Complete payload/scientific data field names/order and arithmetic remain originalsource f912; the new transport envelope never supplies new scientific arguments or fitted values.",
      "Completed means scientific route emitted and source/runtime/transport identities unchanged; it is not controller comparison acceptance. Controller retains all original169 comparison semantics and original earned-error charging.",
      "Completed requires payloadnonnull, errornull, all full axes/array orders, runtimebefore/afternonnull equal, actual source/request identities equal. Refused preserves earned prefixes/details and actual work; complete payload may remain present if only after-seal identity fails but is unavailable for numerical acceptance.",
      "Partial photon/error/frequency prefixes have equal lengths. Full completed groups remain exact fullaxis; active_sigma_rows are an ordered prefix and control_failure retains failing row coarse/refined values. Records are earned in the original interleaved order: sigma0 rows then records0..2, sigma075 rows then records0..2. The records completed-group marker requires all six records, not the first three. Do not fabricate missing fine values, zeros or completed groups.",
      "All highprecision numbers stay60-digit scalarstrings; signed comparator differences may be negative, continuous density may exceed1, and errors are nonnegative. Complete positive densities/probabilities and mixtures must satisfy the original positive/error gate; failure details allow nonpositive failed values to remain diagnostic.",
      "Event probability is dimensionless and in(0,1]; continuous likelihood is ADU density, not constrained<=1. Positiveeventerror<positiveevent and positivejointerror<positivejoint are required. Selected value is same logjoint-logevent only record1.",
      "No independent event-dependence lower threshold is invented. Original numerator witness1e-3 and every event/comparator diagnostic retain separate roles and errors.",
      "Request SHA is actual admitted bytes echoed and independently matched by controller. Do not hardcode new requestSHA inside a reference that the request itself hashes: avoid cyclic source/request identity. Frozen interface/domain/budgets precede implementation; exact implemented script hashes and executable request follow source review. The approved closed axes/settings and shape/domain definitions are self-contained in public request bytes; ignored proposal hashes are documentary ancestry only.",
      "On killed/timeout process or unparseable/nonfinite exception payload, retain bounded raw output/partialfile and execution failure; absence of a typed final envelope cannot be relabeled as successful reference.",
      "No original acceptance values, historical fitted output, native photons/logs, comparison receipts or immutable attempts are consumed by reference arithmetic. Sealed old records are documentary lineage only."
    ],
    "runtime_types": {
      "FileIdentity": {
        "path": "path absolute resolved",
        "bytes": "uint[0,67108864]",
        "sha256": "sha256"
      },
      "PythonIdentity": {
        "version": "boundedUTF8message exactpinned3.14.8seal",
        "implementation": "exactCPython",
        "requested_executable": "path",
        "resolved_executable": "path",
        "executable_bytes": "uint[0,67108864]",
        "executable_sha256": "sha256 exact815b1275bf87e7595fe1cdf3f2efc6ae58e27b842262a7b530de8610ac2760ba",
        "prefix": "path",
        "base_prefix": "path",
        "pyvenv_cfg": "FileIdentity exactcfg11a0241468a7d4cb44419cd7c7a87cd6149c3e90e38b75be6e1a201279621e43",
        "dont_write_bytecode": "bool exactlytrue"
      },
      "MpmathIdentity": {
        "version": "exact1.3.0",
        "module_origin": "path",
        "module_root": "path",
        "backend": "exactpython",
        "decimal_digits": "uint exactly60",
        "binary_precision_bits": "uint actualmp context atdps60",
        "gmpy2_loaded": "bool exactlyfalse",
        "original_subset_map_sha256": "sha256 exact982dc09f2ddf1206cc9cd4467c7be1ce4238a612489b43bf5f5d196a8b2495a8"
      },
      "ImportedModuleIdentity": {
        "name": "boundedASCIIidentifier exactsys.moduleskey",
        "kind": "enum built-in,frozen,source,extension,namespace",
        "origin": "nullable path; source/extension require resolvedregular origin ininventory, built-in/frozen/namespace require null",
        "namespace_paths": "array(path,0,32); nonemptyonlynamespace, sortedunique"
      },
      "RuntimeFingerprint": {
        "schema_version": "uint exactly1",
        "interface_id": "exact temporal-shared-optical-reference-runtime/v1",
        "python": "PythonIdentity",
        "mpmath": "MpmathIdentity",
        "thread_environment": "closed OMP_NUM_THREADS,OPENBLAS_NUM_THREADS,MKL_NUM_THREADS,NUMEXPR_NUM_THREADS,BLIS_NUM_THREADS,VECLIB_MAXIMUM_THREADS,PYTHONDONTWRITEBYTECODE all exactstring'1'",
        "actual_thread_count": "uint exactly1 from/proc/self/status",
        "imported_modules": "array(ImportedModuleIdentity,1,4096) sorteduniquename; complete actualsys.modules set after metadata/helperimports",
        "mapped_libraries": "array(FileIdentity,0,256) sorteduniquepath; complete actual regular native .so mappings subset ofinventory",
        "inventory": "array(FileIdentity,1,4096) sorteduniquepath; totalbytes<=67108864",
        "inventory_sha256": "sha256 canonicalcompactUTF8 list sortedobjectkeys",
        "numeric_configuration": "closed{arithmetic_backend:exactmpmath-libmp-python,decimal_digits:uint60,binary_precision_bits:uintactual,active_native_numeric_backends:emptyarray,gmpy2_loaded:false,byteorder:enumlittle,big}"
      }
    },
    "runtime_semantic_rules": [
      "Seal original complete selected mpmath*.py/license/dist-info metadata subset map against its preserved canonical hash, not merely __init__.py/version. Add every relevant packagefile, importedstdlib source/extension, executable/cfg, currentreference source and mappedregular.so to expanded complete inventory.",
      "Record actual moduleorigin, BACKENDpython and actualprocessThreads1. Environment strings state policy only. Refuse a loaded gmpy2 or native numerical backend rather than silently changing arithmetic ancestry.",
      "Finish metadata/backend discovery before inventory. All mapped regular.so files from/proc/self/maps are independently hashed; no claim that native loaderresolution proves an actual processmapping.",
      "Fingerprint before/after source/science on success and failure. Runtime getter calls noGL nodegeneration, source integrals, Fourier evaluation or likelihood mixture; no package installation, SDK/library/runtime edits or bytecode writes.",
      "Expanded inventory canonical digest is SHA256(json.dumps(inventory,sort_keys=True,separators=(',',':'),allow_nan=False).encode('utf8')); paths are sorted unique resolved absolute strings. Every fileidentity hash is independently checked by controller.",
      "Separate complete imported_modules and mapped_libraries sets must also be exactly equal before/after and externalpreflight/terminal. Unionfileinventory equality alone cannot detect a new imported mpmath module whose file was already in the complete package inventory. Unsupported moduleorigin, disappeared/deleted mapping or unreadable library is an explicit runtime refusal, not an omitted identity.",
      "The controller records pinned native ELF loader --list output before native and in terminal finally as structural dependency resolution only. This does not establish actual native /proc mapping. Reference fingerprint derives its own actual mapped-library set directly from its current process."
    ],
    "unchanged_scientific_settings": {
      "decimal_digits": 60,
      "arithmetic_reservation": "1e-50*(1+abs(quantity))",
      "frequency_and_Fourier_GL_orders": [
        32,
        64
      ],
      "Fourier_panel_width_inverse_electron": "1/4",
      "Fourier_panels": 96,
      "Fourier_cutoff_inverse_electron": 24,
      "sigma_positive_electrons": "3/4",
      "state_masses": [
        "1/4",
        "3/4"
      ],
      "endpoint_relative_empirical_allocation": "3.2e-12",
      "native_photon_comparison": "1e-300+2e-12*abs(polynomial)",
      "frequency_comparison_and_refinement": "1e-300+2e-13*abs(polynomial)",
      "joint_event_selected_log_comparison": "2e-10*(1+abs(reference_log))",
      "reference_log_refinement_max_fraction": "0.05",
      "native_external_aggregate_log_diagnostic": "5e-12+1e-8*(1+abs(native_log))",
      "same_selection_identity": "1e-17 controllercheck",
      "numerator_dependence_witness": "1e-3, syntheticcriterion only",
      "event_dependence_minimum": null
    },
    "proposed_resource_limits": {
      "jobs": 1,
      "threads": 1,
      "compile_cpu_seconds": 180,
      "compile_wall_seconds": 180,
      "native_cpu_seconds": 180,
      "native_wall_seconds": 180,
      "reference_cpu_seconds": 180,
      "reference_wall_seconds": 180,
      "runtime_preflight_cpu_seconds": 180,
      "runtime_preflight_wall_seconds": 180,
      "address_limit_bytes": 1073741824,
      "stdout_limit_bytes_each": 1048576,
      "stderr_limit_bytes_each": 1048576,
      "reference_json_limit_bytes": 1048576,
      "source_tree_bytes": 2097152,
      "source_archive_bytes": 4194304,
      "source_file_bytes": 1048576,
      "attempt_store_bytes": 268435456
    },
    "gates": {
      "current": "proposal only; parent must approve closed contract and immutable request before any implementation",
      "compute": "zero allocated jobs; no imports, syntax checks, tests, builds, native/reference execution or scientific recomputation",
      "preserve": "All48 admitted files/original failedattempts/sealedreceipts remain byte-identical; public/SDK/runtime/package changes forbidden in this subtask",
      "interpretation": "Synthetic numerical control only; physical and observational source/calibration/detector/population/selection/covariance/NEXT17 gates remain open"
    }
  }
}
'''

if __name__ == "__main__":
    raise SystemExit(main())
