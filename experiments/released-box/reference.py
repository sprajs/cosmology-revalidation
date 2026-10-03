"""Bounded independent original-input box reference; empirical, never a certificate.

No native outputs, old fitted values, or old reference modules are read.  The
original binary32 FITS payloads and their exact binary64 transport are admitted
again.  LAPACK completion peers share Cholesky whitening; original-C wide
residuals and a separate eigenvalue determinant route calibrate their numerical
agreement.  Decimal 80/120 arithmetic concerns the reduced system and tails.
"""
from __future__ import annotations

import contextlib
import ctypes
import decimal
from decimal import Decimal as D, localcontext
import hashlib
import io
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import resource
import signal
import stat
import struct
import sys
import time

THREAD_NAMES = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
for _name in THREAD_NAMES:
    os.environ[_name] = "1"


def cap_resource(kind, limit):
    _,hard = resource.getrlimit(kind)
    effective = limit if hard == resource.RLIM_INFINITY else min(limit,hard)
    resource.setrlimit(kind,(effective,effective))


cap_resource(resource.RLIMIT_AS,3221225472)
cap_resource(resource.RLIMIT_CPU,900)
cap_resource(resource.RLIMIT_FSIZE,1048576)
signal.alarm(900)

import numpy as np
import scipy
import scipy.linalg as la
import scipy.linalg.blas
import scipy.linalg.lapack
try:
    import threadpoolctl
except ImportError:
    threadpoolctl = None

INTERFACE = "released-fixed44-box-reference/v1"
REQUEST_SHA256 = "fe6a68858122b0ef590d4f046cea508d3b76b33f3718d47cb60e51c86ae55228"
CONTRACT_SHA256 = "ee885213f48cdb04d3b6d0bada41fea3d65524feb9f33b0d3bda985645aa0741"
LINEAGE_SHA256 = "be07b74a1550dcc3501f332246213cea0a827e7ff1a7bb8ab3ac3f404401d7b6"
CONSTRAINED_SHA256 = "2b8fab46b097dec96162130f85b9e1e675cdeb52af0032a631c666b834e5d499"
STRUCTURE_SHA256 = "e41d7f45428056e8c05e2952dbac54020009422e9a9e4c583aa4158e7eb2c2b1"
ACTIVE = tuple(j for j in range(47) if j != 44)
N, P = 3492, 46
CANONICAL = {
    "C.f64": (97552512, "ef4c2703047e3f8f3b74a77889d9aa74d66df99a249021f2ea9b3757b0c46471"),
    "X.f64": (1312992, "7ac0bcf92b3658ff4af1002cfc75f49200fc1c45bf88dc1f91d19bb39a3ad288"),
    "y.f64": (27936, "e24516870e16f1fc4362f6c8837965163695ce21f425a7ffee7cd549cd43b61b"),
}
ASSETS = {
    "MCMC_utils.py": (2730, "e6ec3d83a9b126d7772ec6dd0d1b58ca757b820acd853fb89218c379cf873841"),
    "lstsq_results.txt": (721, "37d2d423d06b6a2100c47578eb9b1c566a575caf6ade0ed61f6b9586d2a35c95"),
    "allc_shoes_ceph_topantheonwt6.0_112221.fits": (48781440, "a42778672d25df7a559bd2949b1e412b99b2e35dd9a36895f6b38828be019172"),
    "alll_shoes_ceph_topantheonwt6.0_112221.fits": (659520, "9a2ce872dd20ed4fdf5005ce62a805d3eefbbb0b20b0c614013e7c0094009db8"),
    "ally_shoes_ceph_topantheonwt6.0_112221.fits": (17280, "10bb034ca3fff53f6625c9fe47ebd054b37ef45d1140af1c766831cfa6871433"),
    "run_mcmc.py": (2020, "00546b9de2d73cd6315c0333b5f7b07a8b948f9e73056a88e76abe0249464c85"),
}
SETTINGS = {
    "digits": [80, 120], "correction_steps": [1, 3], "diagnostic_digits": 36,
    "postcast_source_residual_digits": 2048, "postcast_correction_steps": 3,
    "observations": N, "original_columns": 47, "active_columns": P,
    "covariance_route": "original-C-LAPACK-Cholesky-with-longdouble-residual-refinement/v1",
    "completion_peers": ["LAPACK-gesdd-SVD", "LAPACK-pivoted-economic-QR"],
    "logdet_C_peer": "original-C-component-LAPACK-eigvalsh-evr/v1",
    "reduced_route": "stdlib-Decimal-LDLT-original-X-transpose-Cinv-X/v1",
    "tail_route": "correlation-valid-union-Chernoff-Decimal/v1",
    "median_route": "correlation-valid-union-bound-original46/v1",
    "reference_refinement_fraction": "0.05",
    "allocations": {"coefficient_absolute": "1e-8", "coefficient_relative": "1e-9",
        "variance_absolute": "2e-12", "variance_relative": "2e-12",
        "quadratic_absolute": "1e-7", "quadratic_relative": "1e-10",
        "log_absolute": "1e-7", "median_absolute": "5e-10"},
    "resources": {"jobs": 1, "threads": 1, "cpu_seconds": 900,
        "wall_seconds": 900, "address_limit_bytes": 3221225472,
        "output_limit_bytes": 1048576},
    "maximum_excluded_mass": "0.1",
    "maximum_covariance_forward_sensitivity": "1e-10",
    "minimum_longdouble_mantissa_bits": 64,
}
QUALIFICATION = {
    "scope": "empirical-original-input-comparison/v1",
    "original_input_certificate": False,
    "conditional_enclosure_scope": "conditional-Gaussian-union-mass-and-median-only/v1",
    "native_outputs_consumed": False,
    "observational_qualification": "blocked-original-contract-source-gaps",
    "limits": [
        "QR and SVD share source data and LAPACK Cholesky whitening; they are numerical peers, not observational independence.",
        "Residual sensitivity uses a LAPACK inverse-norm estimate, not a rigorous original-C inverse bound.",
        "Original-C determinant calibration compares independent Cholesky and eigenvalue algorithms in binary64; reduced Decimal precision alone does not certify C.",
        "Gaussian union bounds are valid for arbitrary correlations conditional on the computed completion; original-input completion errors are separately empirical.",
        "Normalization brackets use nearest Machin/log/base arithmetic and a heuristic Decimal cushion; their 80/120 calibration is empirical, not a certified normalization enclosure.",
        "LAPACK peer q_min values are approximate profile quadratics at reported binary64 peer coefficients, not exact-minimum oracles.",
        "Original source selection, coordinate44 physical identity, calibration, event and cross-covariance gaps remain unresolved.",
    ],
}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def regular(path, expected_size=None, expected_hash=None):
    path = Path(path)
    info = path.lstat()
    if not stat.S_ISREG(info.st_mode):
        raise ValueError("input must be a regular non-symlink file: " + str(path))
    if expected_size is not None and info.st_size != expected_size:
        raise ValueError("input size differs: " + str(path))
    value = digest(path)
    if expected_hash is not None and value != expected_hash:
        raise ValueError("input hash differs: " + str(path))
    return value


def closed_json(path, maximum_bytes):
    if Path(path).stat().st_size > maximum_bytes:
        raise ValueError("bounded JSON input exceeds allocation")
    def pairs(items):
        value = {}
        for key, item in items:
            if key in value:
                raise ValueError("duplicate JSON key")
            value[key] = item
        return value
    def constant(_):
        raise ValueError("nonfinite JSON token")
    return json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=pairs,
                      parse_constant=constant)


def mapped_libraries():
    files = set()
    for line in Path("/proc/self/maps").read_text().splitlines():
        parts = line.split(maxsplit=5)
        if len(parts) == 6 and parts[5].startswith("/") and ".so" in parts[5]:
            path = Path(parts[5])
            if not path.is_file():
                raise ValueError("mapped runtime library unavailable")
            files.add(path.resolve())
    return files


def actual_threadpools():
    if threadpoolctl is not None:
        pools = threadpoolctl.threadpool_info()
        for pool in pools:
            if not isinstance(pool.get("filepath"),str):
                raise ValueError("actual backend filepath is unavailable")
            pool["filepath"] = str(Path(pool["filepath"]).resolve())
        return pools
    # The pinned existing NumPy/SciPy environment lacks threadpoolctl.  Read
    # actual OpenBLAS getter/configuration symbols from already mapped library
    # files.  This does not perform a BLAS calculation or silently load a
    # replacement backend.  Unknown numerical backend families are refused.
    pools = []
    for path in sorted(mapped_libraries(),key=str):
        name = path.name.lower()
        if not name.startswith("lib"):
            continue
        if "openblas" not in name:
            if any(token in name for token in ("mkl","blis","blas","lapack")):
                raise ValueError("actual mapped numerical backend has no frozen introspection route")
            continue
        library = ctypes.CDLL(str(path))
        def symbol(names):
            for name in names:
                try:
                    return getattr(library,name)
                except AttributeError:
                    pass
            raise ValueError("OpenBLAS runtime getter symbol is unavailable")
        getter = symbol(("scipy_openblas_get_num_threads64_","scipy_openblas_get_num_threads",
                         "openblas_get_num_threads64_","openblas_get_num_threads"))
        getter.argtypes,getter.restype = [],ctypes.c_int
        config = symbol(("scipy_openblas_get_config64_","scipy_openblas_get_config",
                         "openblas_get_config64_","openblas_get_config"))
        config.argtypes,config.restype = [],ctypes.c_char_p
        raw = config()
        if raw is None:
            raise ValueError("OpenBLAS runtime configuration is unavailable")
        pools.append({"user_api": "blas", "internal_api": "openblas",
            "filepath": str(path), "num_threads": int(getter()),
            "configuration": raw.decode("ascii"), "introspection": "mapped-OpenBLAS-runtime-getters/v1"})
    return pools


def runtime_fingerprint():
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        np.show_config()
        scipy.show_config()
    pools = actual_threadpools()
    if not pools or any(type(pool.get("num_threads")) is not int or pool["num_threads"] != 1 for pool in pools):
        raise ValueError("actual numerical backend thread count is unavailable or exceeds one")
    files = set()
    for module in tuple(sys.modules.values()):
        name = getattr(module, "__file__", None)
        if name:
            path = Path(name)
            if path.suffix in (".pyc", ".pyo"):
                candidate = Path(importlib.util.source_from_cache(str(path)))
                if candidate.is_file():
                    path = candidate
            if path.is_file():
                files.add(path.resolve())
    for package in (np, scipy):
        root = Path(package.__file__).resolve().parent
        files.update(path.resolve() for path in root.rglob("*.so") if path.is_file())
        sibling = root.parent / (root.name + ".libs")
        if sibling.is_dir():
            files.update(path.resolve() for path in sibling.rglob("*") if path.is_file())
    files.update(mapped_libraries())
    inventory = [{"path": str(path), "bytes": path.stat().st_size, "sha256": digest(path)}
                 for path in sorted(files, key=str)]
    executable = Path(sys.executable)
    return {"python": {"version": sys.version, "implementation": platform.python_implementation(),
                "executable": {"path": str(executable), "resolved_path": str(executable.resolve()),
                    "bytes": executable.stat().st_size, "sha256": digest(executable)}},
        "packages": {"numpy": {"version": np.__version__, "module_root": str(Path(np.__file__).resolve().parent)},
            "scipy": {"version": scipy.__version__, "module_root": str(Path(scipy.__file__).resolve().parent)},
            "decimal": {"precision_policy": [80, 120], "module_path": str(Path(decimal.__file__).resolve())}},
        "thread_environment": {name: os.environ.get(name) for name in THREAD_NAMES},
        "inventory": inventory, "inventory_sha256": hashlib.sha256(canonical_json(inventory)).hexdigest(),
        "numeric_configuration": json.dumps({"show_config": buffer.getvalue(), "threadpools": pools,
            "decimal": {"version": decimal.__version__, "libmpdec_version": decimal.__libmpdec_version__,
                        "rounding": decimal.getcontext().rounding}}, sort_keys=True)}


def exact_float(value):
    return D.from_float(float(value))


def exact_wide(value):
    if not np.isfinite(value):
        raise ValueError("nonfinite wide result")
    fraction, exponent = np.frexp(np.longdouble(value))
    bits = int(np.finfo(np.longdouble).nmant) + 1
    mantissa = int(np.ldexp(fraction, bits))
    exponent = int(exponent) - bits
    with localcontext() as ctx:
        ctx.prec = max(160,bits+abs(exponent)+2)
        return D(mantissa) * (D(2) ** exponent)


def text_number(value):
    if not value.is_finite():
        raise ValueError("nonfinite reference decimal")
    return str(value)


def diagnostic_number(value, rounding=decimal.ROUND_CEILING):
    # Diagnostic residual/norm/error estimates are nonnegative.  Reporting
    # outward36 digits bounds record size without truncating fine arithmetic,
    # completion values or the empirical aggregate comparison errors.
    with localcontext() as ctx:
        ctx.prec,ctx.rounding = SETTINGS["diagnostic_digits"],rounding
        return text_number(+value)


def interval(lower, upper):
    if lower > upper:
        raise ValueError("reversed interval")
    return {"lower": text_number(lower), "upper": text_number(upper)}


def pi_decimal():
    def arctan_inverse(denominator):
        x = D(1) / D(denominator)
        term = x
        total = x
        tolerance = D(10) ** (-decimal.getcontext().prec - 5)
        for k in range(1, 1000):
            term *= -x*x
            addition = term / D(2*k + 1)
            total += addition
            if abs(addition) < tolerance:
                return total
        raise ValueError("bounded Machin pi series did not refine")
    return 16*arctan_inverse(5) - 4*arctan_inverse(239)


def gaussian_union_upper(means, variances, lower, upper, mean_errors=None, variance_errors=None):
    """Outward Chernoff union bound for the exact declared Decimal completion.

    All margin subtraction/division, variance inflation, square, exponent and
    sums are directed.  Decimal sqrt/exp are correctly rounded nearest even
    independently of context rounding: one outward next_plus covers each.
    This bounds tails conditional on declared mean/variance; optional errors
    remain empirical original-input calibration, never certified C bounds.
    """
    mean_errors = [D(0)]*len(means) if mean_errors is None else mean_errors
    variance_errors = [D(0)]*len(means) if variance_errors is None else variance_errors
    with localcontext() as ctx:
        bound,margins_lower,margins_upper = D(0),[],[]
        for mu,var,lo,hi,merror,verror in zip(means,variances,lower,upper,mean_errors,variance_errors):
            if var <= 0 or merror < 0 or verror < 0:
                raise ValueError("tail completion domain differs")
            ctx.rounding = decimal.ROUND_CEILING
            sigma_upper = (var+verror).sqrt().next_plus()
            for endpoint,is_lower,store in ((lo,True,margins_lower),(hi,False,margins_upper)):
                ctx.rounding = decimal.ROUND_FLOOR
                difference_lower = (mu-endpoint-merror) if is_lower else (endpoint-mu-merror)
                margin_lower = difference_lower/sigma_upper
                store.append(margin_lower)
                if margin_lower < 0:
                    term_upper = D(1)
                else:
                    square_lower = margin_lower*margin_lower
                    ctx.rounding = decimal.ROUND_CEILING
                    exponent_upper = -square_lower/D(2)
                    term_upper = exponent_upper.exp().next_plus()/D(2)
                ctx.rounding = decimal.ROUND_CEILING
                bound += term_upper
        return bound.next_plus(),margins_lower,margins_upper


def median_radius_upper(variance, excluded):
    # pi <22/7 supplies a rational upper constant, removing any dependence
    # of this theorem-bound on an unbounded nearest-rounded pi series value.
    with localcontext() as ctx:
        ctx.rounding = decimal.ROUND_CEILING
        sigma_upper = variance.sqrt().next_plus()
        constant_upper = (D(44)/D(7)).sqrt().next_plus()
        return sigma_upper*excluded*constant_upper


def ldlt(matrix):
    p = len(matrix)
    if p == 0 or any(len(row) != p for row in matrix):
        raise ValueError("reduced precision must be a nonempty square matrix")
    if any(not isinstance(value,D) or not value.is_finite() for row in matrix for value in row):
        raise ValueError("reduced precision must contain finite Decimal values")
    if any(matrix[i][j] != matrix[j][i] for i in range(p) for j in range(i)):
        raise ValueError("asymmetric reduced precision")
    lower = [[D(int(i == j)) for j in range(p)] for i in range(p)]
    diagonal = []
    for j in range(p):
        pivot = matrix[j][j] - sum((lower[j][k]**2*diagonal[k] for k in range(j)), D(0))
        if pivot <= 0:
            raise ValueError("reduced original-input precision is not SPD")
        diagonal.append(pivot)
        for i in range(j+1, p):
            lower[i][j] = (matrix[i][j] - sum((lower[i][k]*lower[j][k]*diagonal[k]
                                              for k in range(j)), D(0))) / pivot
    return lower, diagonal


def solve_ldlt(factor, vector):
    lower, diagonal = factor
    p = len(vector)
    if p != len(diagonal) or any(not isinstance(value,D) or not value.is_finite() for value in vector):
        raise ValueError("reduced solve vector domain differs")
    forward = []
    for i in range(p):
        forward.append(vector[i] - sum((lower[i][j]*forward[j] for j in range(i)), D(0)))
    output = [forward[i]/diagonal[i] for i in range(p)]
    for i in range(p-1, -1, -1):
        output[i] -= sum((lower[j][i]*output[j] for j in range(i+1, p)), D(0))
    return output


def fits_payload(path, shape):
    cards, ended, offset = {}, False, 0
    with Path(path).open("rb") as stream:
        for _ in range(32):
            block = stream.read(2880)
            if len(block) != 2880:
                raise ValueError("truncated FITS header")
            offset += 2880
            for j in range(0, 2880, 80):
                card = block[j:j+80].decode("ascii")
                key = card[:8].strip()
                if key == "END":
                    ended = True
                    break
                if card[8:10] == "= ":
                    if key in cards:
                        raise ValueError("duplicate FITS value card")
                    cards[key] = card[10:].split("/", 1)[0].strip()
            if ended:
                break
    if not ended or cards.get("SIMPLE") != "T" or cards.get("BITPIX") != "-32":
        raise ValueError("unsupported original FITS primary image")
    if int(cards.get("NAXIS", "-1")) != len(shape):
        raise ValueError("FITS axes differ")
    axes = tuple(int(cards.get("NAXIS"+str(k), "-1")) for k in range(1, len(shape)+1))
    if axes[::-1] != shape or cards.get("BSCALE", "1") not in ("1", "1.0") or cards.get("BZERO", "0") not in ("0", "0.0"):
        raise ValueError("FITS shape or scaling differs")
    if any(key in cards for key in ("BLANK", "GROUPS")):
        raise ValueError("unsupported FITS data convention")
    count = math.prod(shape)
    payload_end = offset + count*4
    if Path(path).stat().st_size != ((payload_end+2879)//2880)*2880:
        raise ValueError("FITS payload/trailing extent differs")
    with Path(path).open("rb") as stream:
        stream.seek(payload_end)
        if any(stream.read()):
            raise ValueError("nonzero FITS padding")
    return np.memmap(path, dtype=">f4", mode="r", offset=offset, shape=shape)


def admit(input_path, report, tracked):
    raw_hash = regular(input_path)
    args = closed_json(input_path, 32768)
    expected = {"schema_version", "request_sha256", "canonical_directory", "original_source_directory", "source_root"}
    if type(args) is not dict or set(args) != expected or type(args["schema_version"]) is not int or args["schema_version"] != 1:
        raise ValueError("closed reference input schema differs")
    if args["request_sha256"] != REQUEST_SHA256:
        raise ValueError("reference request pin differs")
    for key in ("canonical_directory", "original_source_directory", "source_root"):
        value = args[key]
        if type(value) is not str or not value or len(os.fsencode(value)) > 4096 or not Path(value).is_absolute():
            raise ValueError("reference directory must be bounded absolute path")
        if not Path(value).is_dir() or Path(value).is_symlink():
            raise ValueError("reference directory unavailable or symlink")
    identities = report["identities"]
    identities["reference_input_sha256"] = raw_hash
    identities["reference_script_sha256"] = regular(Path(__file__).resolve())
    adjacent = Path(input_path).resolve().parent
    request_path, contract_path = adjacent/"request.json", adjacent/"contract.snapshot.json"
    identities["request_sha256"] = regular(request_path, expected_hash=REQUEST_SHA256)
    identities["contract_sha256"] = regular(contract_path, expected_hash=CONTRACT_SHA256)
    request, contract = closed_json(request_path, 1048576), closed_json(contract_path, 1048576)
    lineage = Path(args["source_root"])/"experiments/released-ladder/lineage.json"
    constrained = Path(args["source_root"])/"experiments/released-ladder/constrained.json"
    identities["lineage_sha256"] = regular(lineage, expected_hash=LINEAGE_SHA256)
    identities["constrained_sha256"] = regular(constrained, expected_hash=CONSTRAINED_SHA256)
    tracked.update({Path(input_path): raw_hash, Path(__file__).resolve(): identities["reference_script_sha256"],
        request_path: REQUEST_SHA256, contract_path: CONTRACT_SHA256,
        lineage: LINEAGE_SHA256, constrained: CONSTRAINED_SHA256})
    if request["transport"]["active_original_indices"] != list(ACTIVE) or request["requested_marginal"]["original_index"] != 46 or request["requested_marginal"]["active_parameter_index"] != 45 or request["requested_marginal"]["cumulative_probability"] != .5:
        raise ValueError("ordered fixed44/original46 target differs")
    rows = contract["all_original_prior_rows"]
    if len(rows) != 47 or [row["original_index"] for row in rows] != list(range(47)):
        raise ValueError("original support order differs")
    lower = [rows[j]["lower_binary64_hex"] for j in ACTIVE]
    upper = [rows[j]["upper_binary64_hex"] for j in ACTIVE]
    if request["literal_support"] != [{"original_index": j, "lower_binary64_hex": lower[k], "upper_binary64_hex": upper[k]} for k,j in enumerate(ACTIVE)]:
        raise ValueError("literal support pin differs")
    for lo, hi in zip(lower, upper):
        if type(lo) is not str or type(hi) is not str:
            raise ValueError("support hex values must be strings")
        lf, hf = float.fromhex(lo), float.fromhex(hi)
        if not math.isfinite(lf) or not math.isfinite(hf) or lf >= hf or lf.hex() != lo or hf.hex() != hi:
            raise ValueError("support endpoints unsupported")
    report["target"] = {"active_original_indices": list(ACTIVE),
        "fixed_coordinate": {"original_index": 44, "value": 0.0},
        "marginal": {"original_index": 46, "active_index": 45, "cumulative_probability": .5},
        "lower_binary64_hex": lower, "upper_binary64_hex": upper,
        "parameter_measure": "original-active46-Lebesgue-fixed44-point-mass-zero/v1"}
    identities["original_source_sha256"] = {}
    sources = Path(args["original_source_directory"])
    for name, (size, sha) in ASSETS.items():
        path = sources/name
        identities["original_source_sha256"][name] = regular(path, size, sha)
        tracked[path] = sha
    folder = Path(args["canonical_directory"])
    for name, (size, sha) in CANONICAL.items():
        path = folder/name
        identities["canonical_"+name[0]+"_sha256"] = regular(path, size, sha)
        tracked[path] = sha
    c = np.memmap(folder/"C.f64", dtype="<f8", mode="r", shape=(N,N))
    x47 = np.memmap(folder/"X.f64", dtype="<f8", mode="r", shape=(N,47))
    y = np.memmap(folder/"y.f64", dtype="<f8", mode="r", shape=(N,))
    tiny = np.finfo(np.float64).tiny
    for matrix in (c,x47,y):
        if not np.isfinite(matrix).all() or np.any((matrix != 0) & (np.abs(matrix) < tiny)):
            raise ValueError("canonical inputs contain nonfinite or subnormal values")
    original_c = fits_payload(sources/"allc_shoes_ceph_topantheonwt6.0_112221.fits", (N,N))
    original_l = fits_payload(sources/"alll_shoes_ceph_topantheonwt6.0_112221.fits", (47,N))
    original_y = fits_payload(sources/"ally_shoes_ceph_topantheonwt6.0_112221.fits", (N,))
    for start in range(0,N,32):
        stop = min(N,start+32)
        if not np.array_equal(c[start:stop], original_c[start:stop].astype(np.float64)) or not np.array_equal(x47[start:stop], original_l[:,start:stop].T.astype(np.float64)) or not np.array_equal(y[start:stop], original_y[start:stop].astype(np.float64)):
            raise ValueError("canonical bytes do not exactly decode original binary32 source")
        if not np.array_equal(c[start:stop], c[:,start:stop].T):
            raise ValueError("asymmetric original covariance")
    if np.any(np.diag(c) <= 0):
        raise ValueError("nonpositive original covariance diagonal")
    report["gates"]["source_identity"] = True
    return c, np.asarray(x47[:,ACTIVE]), np.asarray(y), lower, upper


def source_components(c, report):
    unseen = set(range(N))
    components = []
    general_sizes = []
    while unseen:
        seed = min(unseen)
        unseen.remove(seed)
        rows, queue = [seed], [seed]
        for row in queue:
            for neighbour in np.flatnonzero(c[row] != 0).tolist():
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    rows.append(neighbour)
                    queue.append(neighbour)
        rows.sort()
        components.append(rows)
        if len(rows) > 1:
            first = float(c[rows[0],rows[1]])
            constant = all(np.all(np.delete(c[row,rows],k) == first) for k,row in enumerate(rows))
            if not constant:
                general_sizes.append(len(rows))
    sizes = [len(rows) for rows in components]
    if sizes != [2593,55,339,143,354]+[1]*8:
        raise ValueError("exact-zero original source structure differs from frozen inventory")
    report["source_structure"] = {"inventory_sha256": STRUCTURE_SHA256,
        "component_count": len(components), "maximum_component_size": max(sizes),
        "maximum_general_component_size": max(general_sizes), "component_sizes": sizes,
        "dense_factor_multiply_add_pairs_per_precision": sum(m*(m-1)*(m+1)//6 for m in general_sizes)}
    return components


def source_peers(c, x, y, components, report):
    factor = la.cholesky(c, lower=True, check_finite=False)
    report["work"]["covariance_factorizations"] += 1
    a = la.solve_triangular(factor,x,lower=True,check_finite=False)
    z = la.solve_triangular(factor,y,lower=True,check_finite=False)
    u,s,vt = la.svd(a,full_matrices=False,lapack_driver="gesdd",check_finite=False)
    rank = int(np.count_nonzero(s > np.finfo(float).eps*max(a.shape)*s[0]))
    if rank != P:
        raise ValueError("original active46 SVD rank is not full")
    svd_mean = (vt.T/s)@(u.T@z)
    svd_variance = np.sum((vt.T/s)**2,axis=1)
    q,r,pivot = la.qr(a,mode="economic",pivoting=True,check_finite=False)
    qr_mean = np.empty(P)
    qr_mean[pivot] = la.solve_triangular(r,q.T@z,check_finite=False)
    inverse_r = la.solve_triangular(r,np.eye(P),check_finite=False)
    qr_variance = np.empty(P)
    qr_variance[pivot] = np.sum(inverse_r**2,axis=1)
    condition = exact_float(s[0]/s[-1])
    def peer(mean,variance,logdet):
        residual = y-x@mean
        white = la.solve_triangular(factor,residual,lower=True,check_finite=False)
        return {"mean": [text_number(exact_float(v)) for v in mean],
            "variance": [text_number(exact_float(v)) for v in variance],
            "q_min": text_number(exact_float(white@white)), "logdet_H": text_number(exact_float(logdet)),
            "rank": rank, "condition2": text_number(condition)}
    report["diagnostics"]["original_input_LAPACK_SVD"] = peer(svd_mean,svd_variance,2*np.log(s).sum())
    report["diagnostics"]["original_input_LAPACK_pivoted_QR"] = peer(qr_mean,qr_variance,2*np.log(np.abs(np.diag(r))).sum())
    with localcontext() as ctx:
        ctx.prec = 140
        chol_logdet = sum((exact_float(v).ln()*2 for v in np.diag(factor)),D(0))
        eig_logdet = D(0)
        component_calibration = []
        for component_index,rows in enumerate(components):
            original_block = np.array(c[np.ix_(rows,rows)],order="F")
            norm = float(np.max(np.sum(np.abs(original_block),axis=1)))
            block_factor = np.array(factor[np.ix_(rows,rows)],order="F")
            rcond,info = la.get_lapack_funcs("pocon",(block_factor,))(block_factor,norm,uplo="L")
            if info != 0 or not np.isfinite(rcond) or rcond <= 0:
                raise ValueError("original-C component inverse norm estimator refused")
            component_calibration.append({"component_index": component_index, "rows": rows,
                "norm": exact_float(norm), "rcond": exact_float(rcond),
                "invnorm": exact_float(1/(norm*rcond))})
            eigenvalues = la.eigvalsh(original_block,
                                      driver="evr",overwrite_a=True,check_finite=False)
            if not np.isfinite(eigenvalues).all() or np.any(eigenvalues <= 0):
                raise ValueError("independent original-C eigenvalue route refused SPD")
            eig_logdet += sum((exact_float(v).ln() for v in eigenvalues),D(0))
            report["work"]["covariance_eigenvalue_components"] += 1
        report["diagnostics"]["logdet_C_peer"] = {"cholesky": text_number(chol_logdet),
            "eigenvalue": text_number(eig_logdet), "absolute_difference": text_number(abs(chol_logdet-eig_logdet)),
            "component_count": len(components)}
        anorm = max(item["norm"] for item in component_calibration)
        invnorm = max(item["invnorm"] for item in component_calibration)
        report["diagnostics"]["covariance_condition_estimate"] = {"norm_infinity": diagnostic_number(anorm),
            "reciprocal_condition_estimate": diagnostic_number(1/(anorm*invnorm),decimal.ROUND_FLOOR),
            "inverse_infinity_norm_estimate": diagnostic_number(invnorm),
            "estimation_scope": "empirical-LAPACK-pocon-not-certified-bound/v1",
            "components": [{"component_index": item["component_index"], "rows": len(item["rows"]),
                "norm_infinity": diagnostic_number(item["norm"]),
                "reciprocal_condition_estimate": diagnostic_number(item["rcond"],decimal.ROUND_FLOOR),
                "inverse_infinity_norm_estimate": diagnostic_number(item["invnorm"])} for item in component_calibration]}
    return factor, chol_logdet, component_calibration


def residual_record(cwide, rhswide, inverse, components, steps, report, completion_digits=None):
    residual = np.array(rhswide,copy=True)
    rhs_count = rhswide.shape[1]
    component_forward,component_records = [],[]
    with localcontext() as ctx:
        ctx.prec = 140
        arithmetic = D(2) ** (-int(np.finfo(np.longdouble).nmant))
        for item in components:
            rows,cnorm,invnorm = item["rows"],item["norm"],item["invnorm"]
            # Exact-zero source components retain every row.  An isolated row
            # with zero X/y has zero inverse basis and zero accumulation error;
            # its large inverse norm cannot contaminate a different block.
            residual[rows] -= cwide[np.ix_(rows,rows)]@inverse[rows]
            norms = [exact_wide(value) for value in np.max(np.abs(residual[rows]),axis=0)]
            wnorm = [exact_wide(value) for value in np.max(np.abs(inverse[rows]),axis=0)]
            bnorm = [exact_wide(value) for value in np.max(np.abs(rhswide[rows]),axis=0)]
            gamma = D(len(rows)+1)*arithmetic/(1-D(len(rows)+1)*arithmetic)
            rounding = [gamma*(cnorm*w+bnorm[j]) for j,w in enumerate(wnorm)]
            forward = [(norms[j]+rounding[j])*invnorm for j in range(rhs_count)]
            backward = [norms[j]/(cnorm*wnorm[j]+bnorm[j]) if cnorm*wnorm[j]+bnorm[j] else D(0) for j in range(rhs_count)]
            # Propagate the outward-reported forward allowances, maintaining
            # identical max identities in top/component diagnostics.
            component_forward.append([D(diagnostic_number(value)) for value in forward])
            component_records.append({"component_index": item["component_index"], "rows": len(rows),
                "residual_infinity_by_rhs": list(map(diagnostic_number,norms)),
                "backward_error_by_rhs": list(map(diagnostic_number,backward)),
                "estimated_forward_error_by_rhs": list(map(diagnostic_number,forward)),
                "inverse_infinity_norm_estimate": diagnostic_number(invnorm)})
        norms = [max(D(row["residual_infinity_by_rhs"][j]) for row in component_records) for j in range(rhs_count)]
        forward = [max(row[j] for row in component_forward) for j in range(rhs_count)]
        backward = [max(D(row["backward_error_by_rhs"][j]) for row in component_records) for j in range(rhs_count)]
        record = {"correction_steps": steps, "residual_infinity_by_rhs": list(map(text_number,norms)),
            "backward_error_by_rhs": list(map(text_number,backward)),
            "estimated_forward_error_by_rhs": list(map(text_number,forward)),
            "inverse_infinity_norm_estimate": diagnostic_number(max(item["invnorm"] for item in components)),
            "component_estimates": component_records}
    if completion_digits is None:
        report["work"]["wide_residual_evaluations"] += 1
        report["diagnostics"]["wide_residual_history"].append(record)
    else:
        report["work"]["postcast_wide_residual_evaluations"] += 1
        record["completion_digits"] = completion_digits
        report["diagnostics"]["postcast_residual_history"].append(record)
    return residual, component_forward


def source_postcast(factor, cwide, components, y, row_entries, beta, digits, unrounded_q, report):
    """Independently solve the actual source residual rounded to binary64.

    Scalar Decimal products form the exact binary32-source/binary64-beta
    residual before its nearest binary64 cast.  A separate allowance covers
    the native wide residual accumulation followed by its binary64 cast.
    No native residual/output is consumed.
    """
    with localcontext() as ctx:
        ctx.prec = SETTINGS["postcast_source_residual_digits"]
        exact_residual = [exact_float(y[i])-sum((v*beta[j] for j,v in row_entries[i]),D(0)) for i in range(N)]
        absolute_terms = [abs(exact_float(y[i]))+sum((abs(v*beta[j]) for j,v in row_entries[i]),D(0)) for i in range(N)]
        cast_residual = np.array([float(value) for value in exact_residual],dtype=np.float64)
    if not np.isfinite(cast_residual).all():
        raise ValueError("postcast binary64 adjusted residual is nonfinite")
    rhswide = np.asarray(cast_residual[:,None],dtype=np.longdouble)
    inverse = np.asarray(la.cho_solve((factor,True),cast_residual[:,None],check_finite=False),dtype=np.longdouble)
    report["work"]["postcast_covariance_solves"] += 1
    residual,_ = residual_record(cwide,rhswide,inverse,components,0,report,digits)
    for steps in range(1,SETTINGS["postcast_correction_steps"]+1):
        inverse += np.asarray(la.cho_solve((factor,True),np.asarray(residual,dtype=np.float64),check_finite=False),dtype=np.longdouble)
        report["work"]["postcast_correction_steps_completed"] += 1
        residual,forward = residual_record(cwide,rhswide,inverse,components,steps,report,digits)
    with localcontext() as ctx:
        ctx.prec = digits
        source_residual = [exact_float(value) for value in cast_residual]
        applied = [exact_wide(inverse[i,0]) for i in range(N)]
        quadratic = sum((source_residual[i]*applied[i] for i in range(N)),D(0))
        if quadratic < 0:
            raise ValueError("negative original-input rounded-residual quadratic")
        inverse_error,native_rounding_allowance = D(0),D(0)
        ctx.rounding = decimal.ROUND_CEILING
        arithmetic = D(2)**(-int(np.finfo(np.longdouble).nmant))
        gamma = D(P+1)*arithmetic/(1-D(P+1)*arithmetic)
        binary64_rounding = D(2)**(-52)
        binary64_subnormal = D(2)**(-1074)
        for index,item in enumerate(components):
            rows = item["rows"]
            inverse_error += sum((abs(source_residual[i]) for i in rows),D(0))*forward[index][0]
            # Include native longdouble dot/subtraction error, its nearest
            # binary64 residual cast (normal/subnormal), and this route's
            # exact-source residual-to-binary64 cast difference.  This works
            # through cancellation/binade crossings, without assuming a
            # fixed ULP at the unrounded residual.
            delta = max((gamma*absolute_terms[i]
                         +binary64_rounding*(abs(exact_residual[i])+gamma*absolute_terms[i])
                         +binary64_subnormal+(max(exact_residual[i],source_residual[i])-min(exact_residual[i],source_residual[i]))
                         if absolute_terms[i] != 0 else D(0) for i in rows),default=D(0))
            applied_norm = max(abs(applied[i]) for i in rows)+forward[index][0]
            native_rounding_allowance += 2*len(rows)*delta*applied_norm+len(rows)*item["invnorm"]*delta**2
        rounding_allowance = max(quadratic,unrounded_q)-min(quadratic,unrounded_q)+native_rounding_allowance
        return quadratic,inverse_error+rounding_allowance,rounding_allowance


def reduced_completion(c, factor_c, cwide, components, x, y, inverse, forward, lo_hex, hi_hex, logdet_c, digits, steps, report):
    with localcontext() as ctx:
        ctx.prec = digits
        ctx.Emin,ctx.Emax = -999999,999999
        xcols = [[(i,exact_float(x[i,j])) for i in np.flatnonzero(x[:,j])] for j in range(P)]
        wide_rows = [[exact_wide(inverse[i,j]) for j in range(P+1)] for i in range(N)]
        hraw = [[sum((v*wide_rows[i][k] for i,v in xcols[j]),D(0)) for k in range(P)] for j in range(P)]
        b = [sum((v*wide_rows[i][P] for i,v in xcols[j]),D(0)) for j in range(P)]
        asymmetry = max(abs(hraw[j][k]-hraw[k][j]) for j in range(P) for k in range(P))
        h = [[(hraw[j][k]+hraw[k][j])/2 for k in range(P)] for j in range(P)]
        factor = ldlt(h)
        report["work"]["reduced_decimal_factorizations"] += 1
        mean = solve_ldlt(factor,b)
        covariance_columns = [solve_ldlt(factor,[D(int(i==j)) for i in range(P)]) for j in range(P)]
        variance = [covariance_columns[j][j] for j in range(P)]
        if any(v <= 0 for v in variance):
            raise ValueError("nonpositive original reduced variance")
        hinvnorm = max(sum(abs(covariance_columns[j][i]) for j in range(P)) for i in range(P))
        row_component = [None]*N
        for component_index,item in enumerate(components):
            for row in item["rows"]:
                row_component[row] = component_index
        xnorms = [[sum((abs(v) for i,v in column if row_component[i] == component_index),D(0))
                   for column in xcols] for component_index in range(len(components))]
        delta_h = max(sum(sum((xnorms[index][j]*forward[index][k] for index in range(len(components))),D(0))
                          for k in range(P)) for j in range(P)) + P*asymmetry/2
        delta_b = max(sum((xnorms[index][j]*forward[index][P] for index in range(len(components))),D(0)) for j in range(P))
        sensitivity = hinvnorm*delta_h
        if sensitivity >= 1:
            raise ValueError("original-C residual sensitivity does not identify reduced completion")
        mean_error = hinvnorm*(delta_b+delta_h*max(map(abs,mean)))/(1-sensitivity)
        variance_error = hinvnorm**2*delta_h/(1-sensitivity)
        logdet_h_error = -P*(1-sensitivity).ln()
        row_entries = [[] for _ in range(N)]
        for j,column in enumerate(xcols):
            for i,v in column:
                row_entries[i].append((j,v))
        def source_q(beta):
            residual = [exact_float(y[i])-sum((v*beta[j] for j,v in row_entries[i]),D(0)) for i in range(N)]
            applied = [wide_rows[i][P]-sum((wide_rows[i][j]*beta[j] for j in range(P)),D(0)) for i in range(N)]
            value = sum((residual[i]*applied[i] for i in range(N)),D(0))
            estimate = D(0)
            for index,item in enumerate(components):
                inverse_error = forward[index][P]+sum((abs(beta[j])*forward[index][j] for j in range(P)),D(0))
                estimate += sum((abs(residual[i]) for i in item["rows"]),D(0))*inverse_error
            if value < 0:
                raise ValueError("negative original-input quadratic")
            return value,estimate
        qmin,qerror = source_q(mean)
        cast_beta = [exact_float(float(v)) for v in mean]
        qpost_unrounded,qpost_unrounded_error = source_q(cast_beta)
        qpost,qposterror,qpost_rounding_allowance = source_postcast(factor_c,cwide,components,y,row_entries,
            cast_beta,digits,qpost_unrounded,report)
        qposterror += qpost_unrounded_error
        hnorm = max(sum(map(abs,row)) for row in h)
        qerror += P*(hnorm+delta_h)*mean_error**2
        logdet_h = sum((v.ln() for v in factor[1]),D(0))
        pi = pi_decimal()
        log2pi = (2*pi).ln()
        lower = [exact_float(float.fromhex(v)) for v in lo_hex]
        upper = [exact_float(float.fromhex(v)) for v in hi_hex]
        excluded,margins_lower,margins_upper = gaussian_union_upper(mean,variance,lower,upper)
        if excluded <= 0 or excluded >= D(SETTINGS["maximum_excluded_mass"]):
            raise ValueError("positive correlated-valid excluded-mass bound exceeds frozen analytic route")
        with localcontext() as outward:
            outward.rounding = decimal.ROUND_FLOOR
            probability_lower = 1-excluded
            log_probability_lower = probability_lower.ln().next_minus()
        log_prior = sum(((upper[j]-lower[j]).ln() for j in range(P)),D(0))
        base = -qmin/2 + P*log2pi/2-logdet_h/2
        # If E<=0.1, z=E sqrt(2pi) obeys Phi(z)-1/2 >= z phi(z)
        # = E exp(-pi E^2) >= E/2.  Thus this median bracket is valid
        # for every correlation pattern, conditional on this completion.
        median_radius = median_radius_upper(variance[45],excluded)
        cushion = D(10)**(-digits+10)*max(D(1),abs(base),abs(log_prior),abs(logdet_c),max(map(abs,mean)))
        with localcontext() as outward:
            outward.rounding = decimal.ROUND_FLOOR
            median_lower = max(lower[45],mean[45]-median_radius-cushion)
            outward.rounding = decimal.ROUND_CEILING
            median_upper = min(upper[45],mean[45]+median_radius+cushion)
        report["gates"]["finite_box_tail"] = True
        output = {"digits": digits, "correction_steps": steps,
            "reduced_forward_sensitivity": text_number(sensitivity),
            "mean": list(map(text_number,mean)), "variance": list(map(text_number,variance)),
            "q_min": text_number(qmin), "q_postcast": text_number(qpost),
            "q_postcast_unrounded_residual": text_number(qpost_unrounded),
            "q_postcast_residual_rounding_allowance": text_number(qpost_rounding_allowance),
            "logdet_H": text_number(logdet_h), "logdet_C": text_number(+logdet_c),
            "endpoint_lower_margins": list(map(text_number,margins_lower)),
            "endpoint_upper_margins": list(map(text_number,margins_upper)),
            "excluded_mass_upper": text_number(excluded),
            "box_probability": interval(probability_lower,D(1)),
            "log_box_probability": interval(log_probability_lower-cushion,D(0)),
            "log_prior_volume": interval(log_prior-cushion,log_prior+cushion),
            "log_Z_relative_box": interval(base+log_probability_lower-cushion,base+cushion),
            "log_prior_normalized_relative_evidence": interval(base+log_probability_lower-log_prior-2*cushion,base-log_prior+2*cushion),
            "log_observation_normalized_evidence": interval(base+log_probability_lower-log_prior-(logdet_c+N*log2pi)/2-3*cushion,base-log_prior-(logdet_c+N*log2pi)/2+3*cushion),
            "median_original46": interval(median_lower,median_upper)}
        calibration = {"mean": mean_error+cushion, "variance": variance_error+cushion,
            "q_min": qerror+cushion, "q_postcast": qposterror+cushion,
            "logdet_H": logdet_h_error+cushion, "cushion": cushion,
            "sensitivity": sensitivity, "median_radius": median_radius+cushion}
        return output,calibration


def mid(value):
    return (D(value["lower"])+D(value["upper"]))/2


def width(value):
    return D(value["upper"])-D(value["lower"])


def endpoint_difference(first, second):
    return max(abs(D(first[key])-D(second[key])) for key in ("lower","upper"))


def qualify(report, calibration):
    coarse,fine = report["coarse"],report["fine"]
    peers = [report["diagnostics"][name] for name in ("original_input_LAPACK_SVD", "original_input_LAPACK_pivoted_QR")]
    with localcontext() as ctx:
        ctx.prec = 140
        mean_error = [max(abs(D(coarse["mean"][j])-D(fine["mean"][j])),
            *(abs(D(peer["mean"][j])-D(fine["mean"][j])) for peer in peers))+calibration["mean"] for j in range(P)]
        variance_error = [max(abs(D(coarse["variance"][j])-D(fine["variance"][j])),
            *(abs(D(peer["variance"][j])-D(fine["variance"][j])) for peer in peers))+calibration["variance"] for j in range(P)]
        errors = {"estimation_scope": "empirical-original-input-calibration-not-certificate/v1",
            "mean_absolute_errors": list(map(text_number,mean_error)),
            "variance_absolute_errors": list(map(text_number,variance_error))}
        for key in ("q_min", "q_postcast", "logdet_H"):
            measured = abs(D(coarse[key])-D(fine[key]))
            if key != "q_postcast":
                measured = max(measured,*(abs(D(peer[key])-D(fine[key])) for peer in peers))
            errors[key+"_absolute_error"] = text_number(measured+calibration[key])
        logc_error = D(report["diagnostics"]["logdet_C_peer"]["absolute_difference"])+abs(D(coarse["logdet_C"])-D(fine["logdet_C"]))+calibration["cushion"]
        errors["logdet_C_absolute_error"] = text_number(logc_error)
        prior_error = endpoint_difference(coarse["log_prior_volume"],fine["log_prior_volume"])+width(fine["log_prior_volume"])/2
        errors["log_prior_volume_absolute_error"] = text_number(prior_error)
        completion_log_error = (D(errors["q_min_absolute_error"])+D(errors["logdet_H_absolute_error"]))/2
        # Inflate every tail margin by the aggregate empirical completion
        # errors.  This calibrates the dependence of finite-box integrals and
        # the median on the original-input completion; it does not turn those
        # empirical errors into certified original-C enclosures.
        robust_excluded,_,_ = gaussian_union_upper(list(map(D,fine["mean"])),list(map(D,fine["variance"])),
            [exact_float(float.fromhex(v)) for v in report["target"]["lower_binary64_hex"]],
            [exact_float(float.fromhex(v)) for v in report["target"]["upper_binary64_hex"]],mean_error,variance_error)
        if robust_excluded <= 0 or robust_excluded >= D(SETTINGS["maximum_excluded_mass"]):
            raise ValueError("empirical completion errors exceed frozen finite-box tail route")
        tail_log_error = abs((1-robust_excluded).ln()-D(fine["log_box_probability"]["lower"]))
        for key,extra in (("log_Z_relative_box",D(0)),
                          ("log_prior_normalized_relative_evidence",prior_error),
                          ("log_observation_normalized_evidence",prior_error+logc_error/2)):
            errors[key+"_absolute_error"] = text_number(endpoint_difference(coarse[key],fine[key])+width(fine[key])/2+completion_log_error+extra+tail_log_error)
        # The original46 mean and variance errors also widen the conditional
        # box-median radius; this remains an empirical original-input estimate.
        with localcontext() as outward:
            outward.rounding = decimal.ROUND_CEILING
            robust_radius = median_radius_upper(D(fine["variance"][45])+variance_error[45],robust_excluded)
        median_error = endpoint_difference(coarse["median_original46"],fine["median_original46"])+max(width(fine["median_original46"])/2,robust_radius)+mean_error[45]
        errors["median_original46_absolute_error"] = text_number(median_error)
        report["refinement"] = errors
        fraction = D("0.05")
        allocation = SETTINGS["allocations"]
        means_ok = all(mean_error[j] <= fraction*(D(allocation["coefficient_absolute"])+D(allocation["coefficient_relative"])*abs(D(fine["mean"][j]))) for j in range(P))
        vars_ok = all(variance_error[j] <= fraction*(D(allocation["variance_absolute"])+D(allocation["variance_relative"])*abs(D(fine["variance"][j]))) for j in range(P))
        q_ok = all(D(errors[key+"_absolute_error"]) <= fraction*(D(allocation["quadratic_absolute"])+D(allocation["quadratic_relative"])*abs(D(fine[key]))) for key in ("q_min","q_postcast"))
        log_budget = fraction*D(allocation["log_absolute"])
        logh_ok = D(errors["logdet_H_absolute_error"]) <= log_budget
        logc_ok = logc_error <= log_budget
        sensitivity_ok = calibration["sensitivity"] <= D(SETTINGS["maximum_covariance_forward_sensitivity"])
        report["gates"]["original_input_empirical_completion"] = means_ok and vars_ok and q_ok and logh_ok and sensitivity_ok
        report["gates"]["original_input_empirical_covariance_logdet"] = logc_ok
        report["gates"]["original_input_empirical_normalization"] = all(D(errors[key+"_absolute_error"]) <= log_budget for key in ("log_Z_relative_box","log_prior_normalized_relative_evidence","log_observation_normalized_evidence"))
        report["gates"]["original46_median_reference"] = median_error <= D(allocation["median_absolute"])
        report["gates"]["reference_refinement"] = all(report["gates"][key] for key in ("original_input_empirical_completion","original_input_empirical_covariance_logdet","original_input_empirical_normalization","original46_median_reference"))
    for key,value in report["gates"].items():
        if key == "runtime_identity":
            continue
        if not value:
            report["blockers"].append("frozen gate failed: "+key)


def empty_report():
    return {"schema_version": 1, "interface_id": INTERFACE, "status": "refused", "accepted": False,
        "error": None, "blockers": [], "identities": {}, "target": None, "settings": SETTINGS,
        "runtime_before": None, "runtime_after": None, "source_structure": None,
        "diagnostics": {"ancestry": "Original released C/X/y; LAPACK QR/SVD share Cholesky whitening. Native factor/QR, source-wide residuals, Decimal reduced arithmetic and source-C eigenvalue determinant provide distinct numerical routes with shared data and mathematical ancestry.",
            "wide_residual_history": [], "postcast_residual_history": []}, "coarse": None, "fine": None, "refinement": None,
        "gates": {key: False for key in ("source_identity","runtime_identity","original_input_empirical_completion",
            "original_input_empirical_covariance_logdet","finite_box_tail","reference_refinement",
            "original_input_empirical_normalization","original46_median_reference")},
        "work": {"elapsed_seconds": 0.0, "cpu_user_seconds": 0.0, "cpu_system_seconds": 0.0,
            "maximum_rss_kib": 0, "covariance_factorizations": 0, "covariance_eigenvalue_components": 0,
            "wide_residual_evaluations": 0, "correction_steps_completed": 0,
            "postcast_covariance_solves": 0, "postcast_wide_residual_evaluations": 0,
            "postcast_correction_steps_completed": 0,
            "reduced_decimal_factorizations": 0, "original_input_rows": N, "native_outputs_consumed": False},
        "qualification": QUALIFICATION}


def main(argv):
    if argv == ["--runtime-fingerprint"]:
        print(json.dumps(runtime_fingerprint(),sort_keys=True,allow_nan=False))
        return 0
    report,tracked,started = empty_report(),{},time.monotonic()
    try:
        if len(argv) != 1:
            raise ValueError("expected one closed reference input JSON or --runtime-fingerprint")
        if np.__version__ != "2.5.3" or scipy.__version__ != "1.18.1":
            raise ValueError("frozen reference NumPy/SciPy versions differ")
        if np.finfo(np.longdouble).nmant+1 < 64 or np.finfo(np.longdouble).maxexp < 16384:
            raise ValueError("wide original-input residual arithmetic is unsupported")
        # Fingerprint imports/backend after admission helpers are loaded and
        # before any scientific route.  Newly loaded modules/libraries fail the
        # final equality gate rather than silently changing runtime ancestry.
        report["runtime_before"] = runtime_fingerprint()
        c,x,y,lower,upper = admit(Path(argv[0]),report,tracked)
        components = source_components(c,report)
        factor,logc,component_calibration = source_peers(c,x,y,components,report)
        rhs = np.column_stack((x,y))
        inverse = np.asarray(la.cho_solve((factor,True),rhs,check_finite=False),dtype=np.longdouble)
        cwide,rhswide = np.asarray(c,dtype=np.longdouble),np.asarray(rhs,dtype=np.longdouble)
        residual,_ = residual_record(cwide,rhswide,inverse,component_calibration,0,report)
        for steps in range(1,4):
            correction = la.cho_solve((factor,True),np.asarray(residual,dtype=np.float64),check_finite=False)
            inverse += np.asarray(correction,dtype=np.longdouble)
            report["work"]["correction_steps_completed"] = steps
            residual,forward = residual_record(cwide,rhswide,inverse,component_calibration,steps,report)
            if steps in (1,3):
                label,digits = ("coarse",80) if steps == 1 else ("fine",120)
                report[label],calibration = reduced_completion(c,factor,cwide,component_calibration,x,y,inverse,forward,lower,upper,logc,digits,steps,report)
        qualify(report,calibration)
    except Exception as error:
        report["error"] = type(error).__name__+": "+str(error)
        report["blockers"].append(report["error"])
    finally:
        for path,sha in tracked.items():
            try:
                regular(path,expected_hash=sha)
            except Exception as error:
                report["gates"]["source_identity"] = False
                report["blockers"].append("post-run source identity: "+str(error))
        try:
            report["runtime_after"] = runtime_fingerprint()
            report["gates"]["runtime_identity"] = report["runtime_before"] == report["runtime_after"]
        except Exception as error:
            report["blockers"].append("post-run runtime identity: "+str(error))
        if not report["gates"]["runtime_identity"]:
            report["blockers"].append("runtime identities differ or incomplete")
        usage = resource.getrusage(resource.RUSAGE_SELF)
        report["work"].update({"elapsed_seconds": time.monotonic()-started,
            "cpu_user_seconds": usage.ru_utime, "cpu_system_seconds": usage.ru_stime,
            "maximum_rss_kib": usage.ru_maxrss})
    report["accepted"] = all(report["gates"].values()) and not report["blockers"] and report["error"] is None
    report["status"] = "accepted" if report["accepted"] else "refused"
    encoded = json.dumps(report,sort_keys=True,allow_nan=False)
    if len(encoded.encode("utf-8")) > SETTINGS["resources"]["output_limit_bytes"]:
        raise ValueError("reference output exceeds frozen capture allocation")
    print(encoded)
    return 0 if report["accepted"] else 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
