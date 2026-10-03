"""Official clik primary reference transport; no spectra or likelihood substitute.

Caller admits exact products/runtime and owns bounded process/failure capture.
Released clik constants are preserved. Calibration prior is applied ONCE to the
three-component combination, with a named source, coordinate measure and an
explicit normalization convention. Module import loads no scientific runtime.
"""
from __future__ import annotations

import hashlib
import json
import math
from numbers import Integral, Real
import os
from pathlib import Path
import stat
import struct

SPECTRA = ("TT", "EE", "BB", "TE", "TB", "EB")
PRODUCTS = {
    "commander_TT": "low_l/commander/commander_dx12_v3_2_29.clik",
    "simall_EE": "low_l/simall/simall_100x143_offlike5_EE_Aplanck_B.clik",
    "plik_TTTEEE": "hi_l/plik/plik_rd12_HM_v22b_TTTEEE.clik",
    "plik_lite_TTTEEE": "hi_l/plik_lite/plik_lite_v22_TTTEEE.clik",
}
PRIOR_SOURCE = {
    "id": "Planck2018-V-section3.3.4-equation39-yP",
    "url": "https://www.aanda.org/articles/aa/full_html/2020/09/aa36386-19/aa36386-19.html",
    "parameter": "A_planck", "source_coordinate": "y_P", "measure": "dA_planck",
    "mean": 1.0, "sigma": 0.0025,
    "A_planck_to_y_P_mapping": "unverified; experiment declaration only",
}
SIMALL_METADATA = {
    "free_calib": ("str", "A_planck"), "lmin": ("int", "2"),
    "nell": ("int", "28"), "nstepsEE": ("int", "3000"),
    "lmax": ("int", "29"), "stepEE": ("float", "0.0001"),
    "lkl_type": ("str", "simall"), "pipeid": ("str", "simall_EE_BB_TE"),
    "unit": ("int", "1"),
}


class LikelihoodRefusal(ValueError):
    """Known earned prefix is available even when a combined score is withheld."""
    def __init__(self, message, record):
        super().__init__(message)
        self.record = record


def finite_real(value, label):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{label}: expected finite real value")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{label}: nonfinite value")
    return result


def fixed_calibration(nuisance):
    if type(nuisance) is not dict or set(nuisance) != {"A_planck"}:
        raise ValueError("first route requires only fixed A_planck=1")
    value = finite_real(nuisance["A_planck"], "A_planck")
    if value != 1.0:
        raise ValueError("first route requires fixed A_planck=1")
    return value


def read_simall_support(product):
    """Read the selected released table's small CLDF metadata, never its arrays."""
    path = Path(product) / "clik/lkl_0/_mdb"
    identity = None
    try:
        if any(p.is_symlink() for p in (path, *path.parents)):
            raise ValueError("SimAll metadata symlink")
        with os.fdopen(os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK), "rb") as stream:
            before = os.fstat(stream.fileno())
            if not stat.S_ISREG(before.st_mode) or before.st_size > 4096:
                raise ValueError("SimAll metadata requires a regular file <=4096 bytes")
            raw = stream.read(4097)
            after = os.fstat(stream.fileno())
        identity = {"path": str(path), "bytes": len(raw),
                    "sha256": hashlib.sha256(raw).hexdigest()}
        if (len(raw) > 4096 or len(raw) != before.st_size
                or (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns,
                    before.st_ctime_ns) !=
                   (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns,
                    after.st_ctime_ns)):
            raise ValueError("SimAll metadata changed during bounded read")
        entries = {}
        for line in raw.decode("ascii").splitlines():
            fields = line.split()
            if len(fields) != 3 or fields[0] in entries:
                raise ValueError("SimAll metadata malformed or duplicate key")
            entries[fields[0]] = tuple(fields[1:])
        if entries != SIMALL_METADATA:
            raise ValueError("SimAll metadata differs from selected released EE table")
        # CLDF reads double with %lg; clik_simall stores it in a C float.
        step = struct.unpack("=f", struct.pack("=f", float(entries["stepEE"][1])))[0]
        return {"metadata_identity": identity, "lmin": 2, "nell": 28,
                "nstepsEE": 3000, "stepEE_native_float32": step, "unit": 1,
                "free_calib": "A_planck", "conditioning": "fixed A_planck=1",
                "rule": "0 <= Cl_uK2*ell*(ell+1)/2/pi/native_step < nstepsEE"}
    except Exception as exc:
        raise LikelihoodRefusal("SimAll support metadata refused", {
            "phase": "simall-support-metadata", "status": "failed",
            "metadata_identity": identity, "error": str(exc)[:2048]}) from exc


def guard_simall(spectra, nuisance, support):
    """Guard the selected C table's index BEFORE the component's compute call.

    At A_planck=1, unit=1 and the admitted unbinned/no-window release, the
    lklbs calibration/selection does not change Cl. Operation order matches
    simall_lkl's double arithmetic; the denominator is its stored C float.
    """
    record = {"phase": "simall-support", "status": "failed",
              "metadata_identity": support["metadata_identity"], "ell": None,
              "Cl_observation": None, "index_argument": None}
    try:
        fixed_calibration(nuisance)
        if "EE" not in spectra or len(spectra["EE"]) < 30:
            raise ValueError("SimAll requires EE ell=0..29")
        ratios = []
        for ell in range(2, 30):
            value = spectra["EE"][ell]
            record.update(ell=ell, Cl_observation=None, index_argument=None)
            if isinstance(value, Real) and not isinstance(value, bool):
                record["Cl_observation"] = scalar_observation(value)
            cl = finite_real(value, f"EE[{ell}]")
            dl = cl * ell * (ell + 1) / 2.0 / math.pi
            ratio = dl / support["stepEE_native_float32"]
            record["index_argument"] = scalar_observation(ratio)
            if not math.isfinite(ratio) or not 0.0 <= ratio < support["nstepsEE"]:
                raise ValueError(f"SimAll EE[{ell}] index outside [0,3000)")
            ratios.append(ratio)
        return {"status": "passed", "metadata_identity": support["metadata_identity"],
                "ells_checked": 28, "minimum_index_argument": min(ratios),
                "maximum_index_argument": max(ratios)}
    except Exception as exc:
        record["error"] = str(exc)[:2048]
        raise LikelihoodRefusal("SimAll support refused before native compute", record) from exc


def validate_contract(lmax, extra_names):
    if len(lmax) != 6 or any(isinstance(n, bool) or not isinstance(n, Integral)
                           or not -1 <= n <= 10000 or int(n) != n for n in lmax):
        raise ValueError("lmax: require six bounded integer maxima")
    names = tuple(extra_names)
    if (len(names) > 256
            or any(not isinstance(n, str) or not n or len(n) > 128 for n in names)
            or len(set(names)) != len(names)):
        raise ValueError("nuisance: invalid runtime name list")
    return names


def build_clik_vector(spectra, nuisance, lmax, extra_names):
    """Cl microkelvin², ell0..lmax TT EE BB TE TB EB, then ordered nuisance."""
    names = validate_contract(lmax, extra_names)
    vector = []
    for name, maximum in zip(SPECTRA, lmax):
        if maximum == -1:
            continue
        if name not in spectra or len(spectra[name]) < maximum + 1:
            raise ValueError(f"{name}: missing ell=0..{maximum} Cl block")
        for ell in range(maximum + 1):
            vector.append(finite_real(spectra[name][ell], f"{name}[{ell}]"))
    for name in names:
        if name not in nuisance:
            raise ValueError(f"nuisance: missing {name}")
        vector.append(finite_real(nuisance[name], name))
    return vector


def calibration_prior(A_planck, convention):
    """ONE declared Gaussian in dA_planck; the source-coordinate alias is open.

    The Gaussian law is on the real line; positive calibration evaluation is
    the likelihood domain. No truncated-normal renormalization is invented.
    Relative penalty omits the displayed constant. It is not a normalized
    density even at A_planck=1, where its value is exactly zero.
    """
    value = finite_real(A_planck, "A_planck")
    if value <= 0 or convention not in ("relative_penalty", "normalized_density"):
        raise ValueError("require positive A_planck and explicit prior convention")
    z = (value - PRIOR_SOURCE["mean"]) / PRIOR_SOURCE["sigma"]
    penalty = -0.5 * z * z
    constant = -math.log(PRIOR_SOURCE["sigma"] * math.sqrt(2.0 * math.pi))
    term = penalty + (constant if convention == "normalized_density" else 0.0)
    finite_real(term, "calibration prior")
    return {"kind": "shared_A_planck_gaussian", "source": dict(PRIOR_SOURCE),
            "value": value, "convention": convention, "relative_logpenalty": penalty,
            "normalization_log_constant": constant, "logprior_term": term,
            "applications": 1}


def combine_terms(component_loglikes, prior=None):
    """Raw likelihood is always distinct; no unnamed zero prior/posterior."""
    if not component_loglikes:
        raise ValueError("no component likelihoods")
    terms = {str(k): finite_real(v, str(k)) for k, v in component_loglikes.items()}
    if any(v <= -1e30 for v in terms.values()):
        raise ValueError("clik invalid sentinel")
    total = math.fsum(terms.values())
    finite_real(total, "combined likelihood")
    finite_real(-2.0 * total, "minus2_loglike")
    target = None
    if prior is not None:
        if (type(prior) is not dict or set(prior) !=
                {"kind", "source", "value", "convention", "relative_logpenalty",
                 "normalization_log_constant", "logprior_term", "applications"}):
            raise ValueError("require named calibration prior record")
        if prior != calibration_prior(prior["value"], prior["convention"]):
            raise ValueError("calibration prior identity/value mismatch")
        target = finite_real(total + prior["logprior_term"], "combined target")
    return {"components": terms, "loglike": total, "calibration_prior": prior,
            "logprior_term": None if prior is None else prior["logprior_term"],
            "logtarget": target, "logposterior": None,
            "posterior_normalization": None,
            "minus2_loglike": -2.0 * total}


def vector_identity(vector):
    """Exact canonical IEEE754 binary64 little-endian input scalar sequence."""
    raw = b"".join(struct.pack("<d", v) for v in vector)
    return {"encoding": "IEEE754-binary64-little-endian", "scalars": len(vector),
            "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


def scalar_observation(value):
    """Nonfinite diagnostics use actual bits, never invalid JSON literals."""
    number = float(value)
    return {"raw_value": number if math.isfinite(number) else None,
            "binary64_le_hex": struct.pack("<d", number).hex(),
            "classification": "finite" if math.isfinite(number) else
                              "nan" if math.isnan(number) else "infinity"}


class PlanckPrimary:
    """One retained owner per batch; constructors perform official selfchecks."""
    def __init__(self, plc_root, highl):
        if highl not in ("plik_TTTEEE", "plik_lite_TTTEEE"):
            raise ValueError("select full Plik or explicitly Plik-lite")
        root = Path(plc_root).resolve(strict=True)
        if not root.is_dir():
            raise ValueError("plc_root: require admitted plc_3.0 directory")
        self.root, self.highl = root, highl
        self.order = ("commander_TT", "simall_EE", highl)
        self.objects, self.contracts = {}, {}
        prefix = {"phase": "initialization", "status": "failed", "components": [],
                  "combined": None, "runtime_version": None}
        try:
            import clik
            prefix["runtime_version"] = clik.version()
            for key in self.order:
                row = {"id": key, "status": "started", "path": str(root / PRODUCTS[key]),
                       "contract": None}
                prefix["components"].append(row)
                if not (root / PRODUCTS[key]).exists():
                    raise FileNotFoundError(root / PRODUCTS[key])
                support = read_simall_support(root / PRODUCTS[key]) if key == "simall_EE" else None
                if support is not None:
                    row["support_metadata"] = support
                obj = clik.clik(str(root / PRODUCTS[key]))
                self.objects[key] = obj
                maxima = tuple(obj.get_lmax())
                names = tuple(obj.get_extra_parameter_names())
                validate_contract(maxima, names)
                # Released wrappers may return numpy integral scalars. Admission
                # precedes lossless conversion; floats and booleans stay refused.
                maxima = tuple(int(n) for n in maxima)
                if key == "simall_EE" and (maxima != (-1, 29, -1, -1, -1, -1)
                                          or names != ("A_planck",)):
                    raise ValueError("SimAll runtime axes differ from selected EE support")
                contract = {"lmax": maxima, "extra_names": names, "path": row["path"],
                            "input_units": "Cl_microkelvin_squared"}
                if support is not None:
                    contract["support_metadata"] = support
                self.contracts[key] = contract
                row.update(status="completed", contract=contract)
            self.runtime_version = prefix["runtime_version"]
        except Exception as exc:
            for row in prefix["components"]:
                if row["status"] == "started":
                    row["status"] = "refused"
                    detail = getattr(exc, "record", None)
                    if isinstance(detail, dict) and detail.get("phase") == "simall-support-metadata":
                        row["support_metadata_error"] = detail
                        row["native_error"] = None
                    else:
                        row["native_error"] = detail
            prefix["error"] = str(exc)[:2048]
            raise LikelihoodRefusal("Planck initialization refused", prefix) from exc

    def evaluate(self, spectra, nuisance, prior):
        rows = [{"id": k, "status": "not_started", "contract": self.contracts[k],
                 "input_vector": None, "raw_observation": None, "error": None}
                for k in self.order]
        record = {"phase": "evaluation", "status": "failed", "component_order": self.order,
                  "component_attempts": rows, "components": {}, "combined": None,
                  "runtime_version": self.runtime_version, "calibration_prior": prior}
        try:
            fixed_calibration(nuisance)
            if prior != calibration_prior(1.0, "relative_penalty"):
                raise ValueError("first route requires declared relative calibration prior at A_planck=1")
            for row in rows:
                key = row["id"]; row["status"] = "started"
                contract = self.contracts[key]
                if key == "simall_EE":
                    row["support_check"] = guard_simall(
                        spectra, nuisance, contract["support_metadata"])
                vector = build_clik_vector(spectra, nuisance, contract["lmax"], contract["extra_names"])
                row["input_vector"] = vector_identity(vector)
                # Pinned clik __call__ converts the supplied list to contiguous double.
                returned = self.objects[key](vector)
                if len(returned) != 1:
                    raise ValueError(f"{key}: expected one clik log likelihood")
                raw = returned[0]
                row["raw_observation"] = scalar_observation(raw)
                value = finite_real(raw, key)
                if value <= -1e30:
                    raise ValueError(f"{key}: clik invalid sentinel")
                record["components"][key] = value; row["status"] = "completed"
            result = combine_terms(record["components"], prior)
            record.update(status="completed", combined=result)
            result["component_attempts"] = rows
            result["runtime_version"] = self.runtime_version
            return result
        except Exception as exc:
            if any(r["status"] == "started" for r in rows):
                failed = next(r for r in rows if r["status"] == "started")
                failed.update(status="refused", error=str(exc)[:2048],
                              native_error=getattr(exc, "record", None))
                if isinstance(exc, LikelihoodRefusal) and exc.record.get("phase") == "simall-support":
                    failed["support_check"] = exc.record
                    failed["native_error"] = None
            record["error"] = str(exc)[:2048]
            raise LikelihoodRefusal("Planck evaluation refused", record) from exc


def validate_config(config):
    if (type(config) is not dict or set(config) !=
            {"backend", "highl", "plc_root", "nuisance", "calibration_prior"}):
        raise ValueError("require exact clik configuration fields")
    if config["backend"] != "clik" or config["highl"] != "plik_lite_TTTEEE":
        raise ValueError("first route is official clik Plik-lite TTTEEE")
    value = fixed_calibration(config["nuisance"])
    if config["calibration_prior"] != {"convention": "relative_penalty"}:
        raise ValueError("first route requires explicit relative calibration prior")
    return calibration_prior(value, "relative_penalty")


def prepare(config):
    """Call once outside the case loop, after caller's data/runtime admission."""
    validate_config(config)
    return PlanckPrimary(config["plc_root"], config["highl"])


def write_record(output_dir, name, record):
    # Serialize before exclusive creation; no invalid/nonfinite success JSON.
    raw = json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n"
    with (Path(output_dir) / name).open("x", encoding="utf-8") as stream:
        stream.write(raw)


def evaluate(cmb, config, output_dir, *, owner):
    """Reuse retained owner; preserve failed component/input/runtime prefix."""
    try:
        prior = validate_config(config)
        if (owner.highl != config["highl"]
                or owner.root != Path(config["plc_root"]).resolve(strict=True)):
            raise ValueError("retained owner/config mismatch")
        if cmb.get("lensed") is not True:
            raise ValueError("require lensed primary CMB spectra")
        ell = cmb["ell"]
        if (not ell or len(ell) > 9999 or any(type(n) is not int for n in ell)
                or list(ell) != list(range(2, len(ell) + 2))):
            raise ValueError("require contiguous integer ell=2 onward")
        spectra = {}
        for name in ("TT", "EE", "TE"):
            values = cmb["spectra"][name.lower()]["Cl_uK2"]
            if len(values) != len(ell):
                raise ValueError(f"{name}: Cl/ell length mismatch")
            spectra[name] = [0.0, 0.0] + [finite_real(v, name) for v in values]
        result = owner.evaluate(spectra, config["nuisance"], prior)
        result.update(backend="official_clik", highl=config["highl"],
                      contracts=owner.contracts, nuisance=dict(config["nuisance"]),
                      conditioning="explicit fixed nuisance; no profiling or inference",
                      lensing_likelihood_included=False)
    except Exception as exc:
        record = getattr(exc, "record", {"phase": "input", "status": "failed",
                    "component_attempts": [], "components": {}, "combined": None,
                    "error": str(exc)[:2048]})
        try:
            write_record(output_dir, "likelihood-failure.json", record)
        except Exception as output_error:
            # Never replace an old file; retain scientific failure AND write failure.
            raise LikelihoodRefusal("likelihood failure capture refused", {
                "earned": record, "output_error": str(output_error)[:2048]}) from exc
        raise LikelihoodRefusal("likelihood refused; prefix retained", record) from exc
    try:
        write_record(output_dir, "likelihood.json", result)
    except Exception as exc:
        raise LikelihoodRefusal("likelihood success capture refused", {
            "phase": "output", "status": "failed", "earned": result,
            "error": str(exc)[:2048]}) from exc
    return result
