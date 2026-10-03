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
from pathlib import Path
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
    """ONE sourced Gaussian in dA_planck, not a fitted-summary likelihood.

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
            "logtarget": target, "logposterior": target,
            "posterior_normalization": "not integrated",
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
                obj = clik.clik(str(root / PRODUCTS[key]))
                self.objects[key] = obj
                maxima = tuple(obj.get_lmax())
                names = tuple(obj.get_extra_parameter_names())
                validate_contract(maxima, names)
                # Released wrappers may return numpy integral scalars. Admission
                # precedes lossless conversion; floats and booleans stay refused.
                maxima = tuple(int(n) for n in maxima)
                contract = {"lmax": maxima, "extra_names": names, "path": row["path"],
                            "input_units": "Cl_microkelvin_squared"}
                self.contracts[key] = contract
                row.update(status="completed", contract=contract)
            self.runtime_version = prefix["runtime_version"]
        except Exception as exc:
            for row in prefix["components"]:
                if row["status"] == "started":
                    row["status"] = "refused"
                    row["native_error"] = getattr(exc, "record", None)
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
            for row in rows:
                key = row["id"]; row["status"] = "started"
                contract = self.contracts[key]
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
            record["error"] = str(exc)[:2048]
            raise LikelihoodRefusal("Planck evaluation refused", record) from exc


def validate_config(config):
    if (type(config) is not dict or set(config) !=
            {"backend", "highl", "plc_root", "nuisance", "calibration_prior"}):
        raise ValueError("require exact clik configuration fields")
    if config["backend"] != "clik" or config["highl"] != "plik_lite_TTTEEE":
        raise ValueError("first route is official clik Plik-lite TTTEEE")
    if type(config["nuisance"]) is not dict or "A_planck" not in config["nuisance"]:
        raise ValueError("require explicit A_planck/nuisance dictionary")
    if (type(config["calibration_prior"]) is not dict
            or set(config["calibration_prior"]) != {"convention"}):
        raise ValueError("require explicit calibration prior convention")
    return calibration_prior(config["nuisance"]["A_planck"],
                             config["calibration_prior"]["convention"])


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
