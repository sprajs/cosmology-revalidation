"""Score four retained CLASS products with one official Planck primary owner.

Parent admits exact data/runtime/configuration and an outer 900 s watchdog.
This source calls no CLASS solver and never combines BAO/SN or fits parameters.
"""
import argparse
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import re
import resource
import signal
import stat
import sys
import time
import types

CASES = [("anchor", "0.9660499"), ("precision", "0.9660499"),
         ("ns-minus", "0.9610499"), ("ns-plus", "0.9710499")]
PRODUCTS = [("commander_TT", "low_l/commander/commander_dx12_v3_2_29.clik"),
            ("simall_EE", "low_l/simall/simall_100x143_offlike5_EE_Aplanck_B.clik"),
            ("plik_lite_TTTEEE", "hi_l/plik_lite/plik_lite_v22_TTTEEE.clik")]
LIMITS = {"case_wall_seconds": 180, "total_wall_seconds": 900,
          "address_bytes": 2147483648, "attempt_bytes": 2147483648,
          "file_bytes": 134217728, "max_files": 128}
SELFCHECK_CRITERION = {"max_abs_printed_difference": 1e-6}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def pairs(items):
    result = {}
    for key, value in items:
        require(key not in result, "duplicate JSON key")
        result[key] = value
    return result


def read_config(path):
    with Path(path).open("rb") as stream:
        raw = stream.read(2097153)
    require(len(raw) <= 2097152, "score configuration byte bound")
    value = json.loads(raw, object_pairs_hook=pairs,
                       parse_constant=lambda x: (_ for _ in ()).throw(ValueError("nonfinite JSON")))
    return value, raw


def checked_file(pin, deadline, *, content=False):
    require(type(pin) is dict and set(pin) == {"path", "bytes", "sha256"}, "exact file pin")
    require(type(pin["bytes"]) is int and 0 <= pin["bytes"] < 2 ** 63, "file pin byte domain")
    require(type(pin["sha256"]) is str and re.fullmatch(r"[0-9a-f]{64}", pin["sha256"]), "file pin digest")
    path = Path(pin["path"])
    require(path.is_absolute() and ".." not in path.parts and len(str(path).encode()) <= 4096,
            "absolute file pin path")
    require(all(p.is_dir() and not p.is_symlink() for p in path.parents), "file pin ancestors")
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    try:
        before = os.fstat(fd)
        require(stat.S_ISREG(before.st_mode) and before.st_size == pin["bytes"], "file pin type/size")
        digest, count, chunks = hashlib.sha256(), 0, []
        if content:
            require(pin["bytes"] <= LIMITS["file_bytes"], "consumed JSON/source byte bound")
        while count <= pin["bytes"]:
            require(time.monotonic() < deadline, "whole scorer input-verification deadline")
            chunk = os.read(fd, min(65536, pin["bytes"] + 1 - count))
            if not chunk:
                break
            digest.update(chunk); count += len(chunk)
            if content:
                chunks.append(chunk)
        after, linked = os.fstat(fd), path.lstat()
        facts = lambda x: (x.st_dev, x.st_ino, x.st_size, x.st_mtime_ns, x.st_ctime_ns)
        observed = {"path": str(path), "bytes": count, "sha256": digest.hexdigest()}
        if observed != pin or facts(before) != facts(after) or facts(after) != facts(linked):
            exc = ValueError("file consumed identity/stat differs")
            exc.observed_identity, exc.expected_identity = observed, pin
            raise exc
        return (b"".join(chunks) if content else None), observed
    finally:
        os.close(fd)


def load_module(name, pin, deadline):
    require(name not in sys.modules, "module already present before admission: " + name)
    require(pin["bytes"] <= 262144, "adapter source byte bound")
    raw, observed = checked_file(pin, deadline, content=True)
    module = types.ModuleType(name)
    module.__file__, module.__package__ = pin["path"], ""
    sys.modules[name] = module
    # Execute the exact admitted RAM bytes, not a second mutable path read.
    exec(compile(raw, pin["path"], "exec"), module.__dict__)
    return module, observed


def pinned_json(pin, deadline):
    raw, _ = checked_file(pin, deadline, content=True)
    return json.loads(raw, object_pairs_hook=pairs,
                      parse_constant=lambda x: (_ for _ in ()).throw(ValueError("nonfinite JSON")))


def product_tree(plc_root, files, deadline):
    """Refuse undeclared CLDF keys/subtrees before their optional native lookup."""
    roots = [Path(plc_root) / relative for _, relative in PRODUCTS]
    require(type(files) is list and 0 < len(files) <= 8192, "released file inventory bound")
    expected_files, expected_dirs = set(), set(roots)
    for pin in files:
        path = Path(pin["path"])
        owners = [root for root in roots if root in path.parents]
        require(len(owners) == 1 and path not in expected_files, "released file ownership/duplicate")
        expected_files.add(path)
        parent = path.parent
        while parent != owners[0]:
            expected_dirs.add(parent)
            parent = parent.parent
    observed_files, observed_dirs, pending, entries = set(), set(), list(roots), 0
    while pending:
        require(time.monotonic() < deadline, "released directory admission deadline")
        folder = pending.pop()
        require(stat.S_ISDIR(folder.lstat().st_mode), "released product directory type")
        observed_dirs.add(folder)
        with os.scandir(folder) as stream:
            for entry in stream:
                entries += 1
                require(entries <= 16384, "released directory entry bound")
                path, mode = Path(entry.path), entry.stat(follow_symlinks=False).st_mode
                require(len(str(path).encode()) <= 4096, "released member path bound")
                if stat.S_ISDIR(mode):
                    require(path in expected_dirs, "undeclared released directory")
                    pending.append(path)
                else:
                    require(stat.S_ISREG(mode) and path in expected_files,
                            "undeclared or nonregular released file")
                    observed_files.add(path)
    require(observed_files == expected_files and observed_dirs == expected_dirs,
            "released directory inventory differs")
    return {"files": sorted(map(str, observed_files)), "directories": sorted(map(str, observed_dirs))}


def error_record(exc):
    result = {"kind": type(exc).__name__, "message": str(exc)[:2048]}
    for name in ("record", "observed_identity", "expected_identity"):
        if hasattr(exc, name):
            result[name] = getattr(exc, name)
    return result


def arm(deadline):
    seconds = min(LIMITS["case_wall_seconds"], deadline - time.monotonic())
    require(seconds > 0, "expired scorer phase deadline")
    # SIG_DFL is deliberately a kernel termination, also while inside a C call.
    # No returned prefix is invented if this kills initialization/evaluation.
    signal.signal(signal.SIGALRM, signal.SIG_DFL)
    signal.setitimer(signal.ITIMER_REAL, seconds)


def disarm():
    signal.setitimer(signal.ITIMER_REAL, 0)


def initialization(owner_factory, attempt, transport, c_api, deadline):
    log = attempt / "initialization.stdout.log"
    saved = os.dup(1)
    old_size = resource.getrlimit(resource.RLIMIT_FSIZE)
    owner = None
    try:
        require(old_size[1] == resource.RLIM_INFINITY or old_size[1] >= 1048576, "initialization log hard limit")
        sys.stdout.flush(); c_api.flush()
        with log.open("xb", buffering=0) as stream:
            os.dup2(stream.fileno(), 1)
            resource.setrlimit(resource.RLIMIT_FSIZE, (1048576, old_size[1]))
            arm(deadline)
            try:
                owner = owner_factory()
            finally:
                c_api.flush()
                disarm()
        return owner
    except BaseException as exc:
        # A returned owner remains earned if only subsequent capture/flush fails.
        exc.initialized_owner = owner
        raise
    finally:
        os.dup2(saved, 1); os.close(saved)
        resource.setrlimit(resource.RLIMIT_FSIZE, old_size)
        if log.exists():
            log.chmod(0o444)


def selfchecks(raw, plc_root, criterion):
    require(type(criterion) is dict and criterion == SELFCHECK_CRITERION,
            "unchanged initialization selfcheck criterion")
    limit = criterion["max_abs_printed_difference"]
    require(type(limit) in (int, float) and math.isfinite(limit) and limit >= 0,
            "selfcheck difference criterion domain")
    pattern = r"Checking likelihood '([^']+)' on test data\. got (\S+) expected (\S+) \(diff ([^)]+)\)"
    found = re.findall(pattern, raw.decode("utf-8", "strict"))
    expected = [str(Path(plc_root) / path) for _, path in PRODUCTS]
    require([x[0] for x in found] == expected, "three ordered released initialization selfchecks required")
    rows = []
    for (name, got, target, difference), (identifier, _) in zip(found, PRODUCTS):
        values = [float(got), float(target), float(difference)]
        require(all(math.isfinite(v) for v in values), "finite printed selfcheck values")
        rows.append({"id": identifier, "path": name, "got_printed": got, "expected_printed": target,
                     "difference_printed": difference, "criterion": criterion,
                     "passed": abs(values[2]) <= limit, "format": "official C %g diagnostic",
                     "unrounded_selfcheck_values_available": False})
    return rows


def score(config_path, attempt_path, deadline_utc):
    started = time.monotonic()
    config, raw_config = read_config(config_path)
    require(type(config) is dict and set(config) == {"schema", "scorer", "transport", "adapter", "c_api",
            "library", "python", "prediction_record", "cases", "likelihood", "admission",
            "admission_files", "selfcheck_criteria", "limits"}, "closed scoring configuration")
    require(config["schema"] == "lcdm-retained-primary-score-config/v1" and config["limits"] == LIMITS,
            "scorer schema/unchanged limits")
    require(config["selfcheck_criteria"] == SELFCHECK_CRITERION,
            "unchanged initialization selfcheck criterion")
    require([(x["id"], x["ns"]) for x in config["cases"]] == CASES and
            all(set(x) == {"id", "ns", "products"} for x in config["cases"]), "exact four retained cases")
    required = [config[k] for k in ("scorer", "transport", "adapter", "c_api", "library", "python",
                                   "admission", "admission_files", "selfcheck_criteria")]
    require(all(x is not None for x in required), "actual runtime/source/admission pins or criterion not supplied")
    utc = datetime.datetime.fromisoformat(deadline_utc)
    require(utc.tzinfo is not None, "root deadline UTC offset")
    remaining = (utc - datetime.datetime.now(datetime.timezone.utc)).total_seconds()
    require(remaining > 3, "expired root score grant")
    deadline = min(started + LIMITS["total_wall_seconds"], started + remaining)
    require(Path(config["scorer"]["path"]) == Path(__file__).absolute(), "scorer source path")
    checked_file(config["scorer"], deadline)
    require(Path(config["python"]["path"]) == Path("/proc/self/exe").resolve(strict=True),
            "selected Python executable identity")
    checked_file(config["python"], deadline)
    require(set(config["transport"]) == {"run", "theory"}, "frozen prediction transport sources")
    _, _ = load_module("theory", config["transport"]["theory"], deadline)
    transport, _ = load_module("retained_class_transport", config["transport"]["run"], deadline)
    prediction = pinned_json(config["prediction_record"], deadline)
    require(prediction["schema"] == "lcdm-external-class-attempt/v1" and prediction["status"] == "completed"
            and prediction["gates"]["execution"] == "passed", "completed retained CLASS attempt")
    require(len(prediction["cases"]) == 4, "completed four-case prediction cardinality")
    for case, original in zip(config["cases"], prediction["cases"]):
        require((original["id"], original["ns"]) == (case["id"], case["ns"])
                and original["status"] == "completed" and original["gates"]["theory"] == "passed"
                and original["products"] == case["products"], "actual case/product record binding")
    require(config["likelihood"]["calibration_prior"] == {"convention": "relative_penalty"},
            "one explicit relative calibration prior")
    require(config["likelihood"]["nuisance"] == {"A_planck": 1.0}, "explicit fixed calibration")
    attempt = transport.path_without_symlinks(attempt_path)
    attempt.mkdir(mode=0o700)
    record = {"schema": "lcdm-retained-primary-score-attempt/v1", "status": "failed",
              "config": {"path": str(Path(config_path).absolute()), "bytes": len(raw_config),
                         "sha256": hashlib.sha256(raw_config).hexdigest()},
              "prediction_record": config["prediction_record"],
              "cases": [{**x, "status": "not_started", "likelihood": None} for x in config["cases"]],
              "failures": [],
              "initialization": None, "selfchecks": None, "cleanup": None,
              "admission_receipt": config["admission"], "inputs_before": [], "inputs_after": [],
              "gates": {"execution": "unassessed", "numerical": "unassessed", "inference": "blocked",
                        "interpretation": "blocked"},
              "scope": "three official primary terms plus one calibration prior, fixed other coordinates",
              "BAO_SN_combination": False, "posterior_normalization": None,
              "runtime_ABI_independent_qualification": None, "outer_watchdog_required": True}
    c_api, adapter, owner, product_files = None, None, None, None
    pins = []
    try:
        transport.write_new(attempt / "config.snapshot.json", raw_config)
        require(type(config["admission_files"]) is list and 0 < len(config["admission_files"]) <= 8192,
                "bounded parent-owned data/runtime inventory")
        pins = [record["config"], config["scorer"], config["python"], config["transport"]["run"],
                config["transport"]["theory"], config["adapter"], config["c_api"], config["admission"],
                config["library"], config["prediction_record"]] + config["admission_files"] + [
                x["products"] for x in config["cases"]]
        unique = {}
        for pin in pins:
            require(pin["path"] not in unique or pin == unique[pin["path"]], "conflicting repeated file pin")
            unique[pin["path"]] = pin
        pins = list(unique.values())
        for pin in pins:
            record["inputs_before"].append(checked_file(pin, deadline)[1])
        # Receipt hash is source/data admission lineage, not inferred qualification.
        admission = pinned_json(config["admission"], deadline)
        plc_root = Path(config["likelihood"]["plc_root"])
        require(admission["library"] == config["library"] and admission["python"] == config["python"]
                and admission["plc_root"] == str(plc_root)
                and admission["file_inventory"] in config["admission_files"],
                "runtime/data admission binding differs")
        product_manifest = pinned_json(admission["file_inventory"], deadline)
        require(product_manifest["kind"] == "exact-three-product-file-inventory"
                and product_manifest["plc_root"] == str(plc_root)
                and product_manifest["order"] == [identifier for identifier, _ in PRODUCTS],
                "three-product inventory identity/order")
        product_files = product_manifest["files"]
        declared = [pin for pin in config["admission_files"] if any(
            plc_root / relative in Path(pin["path"]).parents for _, relative in PRODUCTS)]
        require(sorted(declared, key=lambda pin: pin["path"]) ==
                sorted(product_files, key=lambda pin: pin["path"]), "released file pins differ from admission")
        record["released_products_before"] = product_tree(plc_root, product_files, deadline)
        resource.setrlimit(resource.RLIMIT_AS, (LIMITS["address_bytes"],) * 2)
        for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
                     "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
            os.environ[name] = "1"
        c_api, _ = load_module("clik", config["c_api"], deadline)
        record["C_transport"] = c_api.configure(config["library"])
        adapter, _ = load_module("likelihood", config["adapter"], deadline)
        transport.write_new(attempt / "initialization-start.json", transport.json_bytes({
            "status": "started_no_return", "component_order": [x[0] for x in PRODUCTS],
            "kernel_wall_seconds": min(180, deadline - time.monotonic()), "combined": None}))
        try:
            owner = initialization(lambda: adapter.prepare(config["likelihood"]), attempt,
                                   transport, c_api, deadline)
            record["initialization"] = {"status": "returned", "contracts": owner.contracts,
                                        "runtime_version": owner.runtime_version}
        except BaseException as exc:
            record["initialization"] = error_record(exc)
            returned_owner = getattr(exc, "initialized_owner", None)
            if returned_owner is not None:
                record["initialization"]["earned_owner"] = {
                    "contracts": returned_owner.contracts, "runtime_version": returned_owner.runtime_version}
            transport.write_new(attempt / "initialization-failure.json", transport.json_bytes(record["initialization"]))
            raise
        raw, _ = transport.file_bytes(attempt / "initialization.stdout.log", 1048576)
        record["selfchecks"] = selfchecks(raw, config["likelihood"]["plc_root"], config["selfcheck_criteria"])
        record["gates"]["numerical"] = ("printed-selfchecks-only" if all(x["passed"] for x in record["selfchecks"])
                                         else "failed")
        require(all(x["passed"] for x in record["selfchecks"]), "released printed selfcheck criterion refused")
        transport.write_new(attempt / "initialization.json", transport.json_bytes({
            "initialization": record["initialization"], "selfchecks": record["selfchecks"]}))
        for case, row in zip(config["cases"], record["cases"]):
            product = pinned_json(case["products"], deadline)
            folder = attempt / case["id"]
            folder.mkdir(mode=0o700)
            row["status"] = "started_no_return"
            transport.write_new(folder / "evaluation-start.json", transport.json_bytes({
                **row, "kernel_wall_seconds": min(180, deadline - time.monotonic())}))
            arm(deadline)
            try:
                result = adapter.evaluate(product["cmb"], config["likelihood"], folder, owner=owner)
                row.update(status="completed", likelihood={key: result[key] for key in
                           ("components", "loglike", "calibration_prior", "logprior_term", "logtarget")})
            except BaseException as exc:
                row.update(status="refused", failure=error_record(exc))
                raise
            finally:
                disarm()
            transport.inventory(attempt, LIMITS)
        anchor = record["cases"][0]["likelihood"]
        record["conditional_ns_scan"] = [{"id": row["id"], "ns": row["ns"],
            "delta_loglike_from_anchor": row["likelihood"]["loglike"] - anchor["loglike"],
            "delta_logtarget_from_anchor": row["likelihood"]["logtarget"] - anchor["logtarget"],
            "component_deltas": {key: value - anchor["components"][key]
                                 for key, value in row["likelihood"]["components"].items()}}
            for row in record["cases"] if row["id"] != "precision"]
        record["precision_likelihood_difference"] = {
            "component_deltas": {key: value - anchor["components"][key]
                for key, value in record["cases"][1]["likelihood"]["components"].items()},
            "delta_loglike": record["cases"][1]["likelihood"]["loglike"] - anchor["loglike"],
            "certified_numerical_error_bound": None}
        record["status"] = "completed"
        record["gates"].update(execution="passed", numerical="empirical-precision-difference-only")
    except BaseException as exc:
        record["failures"].append(error_record(exc))
        record["gates"]["execution"] = "failed"
    finally:
        disarm()
        if c_api is not None:
            try:
                arm(deadline)
                record["cleanup"] = c_api.close_all()
                require(all(x["status"] in ("returned", "already_closed") for x in record["cleanup"]),
                        "native cleanup refused")
            except BaseException as exc:
                record["status"] = "failed"; record["gates"]["execution"] = "failed"
                record["failures"].append(error_record(exc))
            finally:
                disarm()
        for pin in pins:
            try:
                record["inputs_after"].append(checked_file(pin, deadline)[1])
            except BaseException as exc:
                record["status"] = "failed"; record["gates"]["execution"] = "failed"
                record["failures"].append(error_record(exc))
                break
        try:
            if product_files is not None:
                record["released_products_after"] = product_tree(plc_root, product_files, deadline)
                require(record.get("released_products_before") == record["released_products_after"],
                        "terminal released directory drift")
            require(record["inputs_before"] == record["inputs_after"], "terminal score source/data drift")
            record["output_inventory"] = transport.seal_inventory(attempt, LIMITS)
            require(not record["output_inventory"]["errors"], "terminal score output sealing")
            require(time.monotonic() < deadline, "whole score final deadline")
        except BaseException as exc:
            record["status"] = "failed"; record["gates"]["execution"] = "failed"
            record["failures"].append(error_record(exc))
        record["wall_seconds"] = time.monotonic() - started
        record["wall_scope"] = "through source/data recheck and output seal, before terminal JSON/write"
        observations = {"before": record["inputs_before"], "after": record["inputs_after"]}
        observed_raw = transport.json_bytes(observations)
        require(len(observed_raw) <= LIMITS["file_bytes"], "input-observation sidecar bound")
        record["input_observations"] = transport.write_new(attempt / "input-observations.json", observed_raw)
        record["inputs_before"] = {"observed_files": len(observations["before"])}
        record["inputs_after"] = {"observed_files": len(observations["after"])}
        transport.inventory(attempt, LIMITS)
        if time.monotonic() >= deadline:
            record["status"] = "failed"; record["gates"]["execution"] = "failed"
            record["failures"].append({"kind": "TimeoutError", "message": "whole score deadline after identity sidecar"})
        record["wall_seconds"] = time.monotonic() - started
        record["wall_scope"] = "through source/data recheck, raw output seal and identity sidecar, before terminal JSON/write"
        raw = transport.json_bytes(record)
        require(len(raw) <= 2097152, "score terminal receipt byte reserve")
        transport.write_new(attempt / "record.json", raw)
    print(json.dumps({"status": record["status"], "record": str(attempt / "record.json"),
                      "wall_seconds_through_record_write": time.monotonic() - started}, allow_nan=False))
    return 0 if record["status"] == "completed" else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--attempt", required=True)
    parser.add_argument("--deadline-utc", required=True)
    args = parser.parse_args()
    return score(args.config, args.attempt, args.deadline_utc)


if __name__ == "__main__":
    raise SystemExit(main())
