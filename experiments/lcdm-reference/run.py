"""Bounded external CLASS predictions and a small frozen ns scan.

Run only after the source/configuration review and parent's explicit job grant.
No build, installation or cosmological solver implementation occurs here.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import resource
import selectors
import signal
import stat
import subprocess
import sys
import time

import theory

REVISION = "0ceb7a9a4c1e444ef5d5d56a8328a0640be91b18"
POINT_SHA = "549160924345dad88be62b182c53e13f419a47dfabc032fe5303ad57a2ff2d7d"
PRECISION_SHA = "b6614adedad02f6e22cd42cc9274d048778e5c546c7a8fcfa282b30e24997199"
CASE_DESIGNS = [{"id": "anchor", "ns": "0.9660499", "precision": False},
                {"id": "precision", "ns": "0.9660499", "precision": True},
                {"id": "ns-minus", "ns": "0.9610499", "precision": False},
                {"id": "ns-plus", "ns": "0.9710499", "precision": False}]


def require(condition, message):
    if not condition:
        raise ValueError(message)


class FileIdentityError(ValueError):
    def __init__(self, message, observed, expected=None):
        super().__init__(message)
        self.observed_identity = observed
        self.expected_identity = expected


def failure(exc, limit=2048):
    result = {"kind": type(exc).__name__, "message": str(exc)[:limit]}
    if isinstance(exc, FileIdentityError):
        result["observed_identity"] = exc.observed_identity
        result["expected_identity"] = exc.expected_identity
    return result


def pairs(items):
    value = {}
    for key, item in items:
        require(key not in value, "duplicate JSON key")
        value[key] = item
    return value


def json_bytes(value):
    return (json.dumps(value, sort_keys=True, allow_nan=False, indent=2) + "\n").encode()


def path_without_symlinks(path):
    path = Path(path).absolute()
    require(".." not in path.parts and len(str(path).encode()) <= 4096, "invalid path")
    for parent in reversed(path.parents):
        require(parent.is_dir() and not parent.is_symlink(), "symlink/non-directory ancestor")
    return path


def file_bytes(path, limit, expected=None):
    path = path_without_symlinks(path)
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    try:
        before = os.fstat(fd)
        require(stat.S_ISREG(before.st_mode) and before.st_size <= limit, "regular file byte bound")
        chunks, count = [], 0
        while True:
            chunk = os.read(fd, min(65536, limit + 1 - count))
            if not chunk:
                break
            chunks.append(chunk)
            count += len(chunk)
            require(count <= limit, "consumed file byte bound")
        after = os.fstat(fd)
        linked = path.lstat()
        facts = lambda x: (x.st_dev, x.st_ino, x.st_size, x.st_mtime_ns, x.st_ctime_ns)
        raw = b"".join(chunks)
        identity = {"path": str(path), "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest(),
                    "device": after.st_dev, "inode": after.st_ino,
                    "mode": oct(stat.S_IMODE(after.st_mode)), "mtime_ns": after.st_mtime_ns}
        if facts(before) != facts(after) or facts(after) != facts(linked) or count != after.st_size:
            raise FileIdentityError("file changed during read", identity, expected)
        if expected is not None:
            require(type(expected["bytes"]) is int and expected["bytes"] >= 0,
                    "declared file byte count type")
            if not all(identity.get(k) == expected[k] for k in ("path", "bytes", "sha256")):
                raise FileIdentityError("declared file identity differs", identity, expected)
        return raw, identity
    finally:
        os.close(fd)


def write_new(path, raw):
    path = path_without_symlinks(path)
    with path.open("xb") as stream:
        stream.write(raw)
    path.chmod(0o444)
    return {"path": str(path), "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


def load_json(path, limit=2097152, expected=None):
    raw, identity = file_bytes(path, limit, expected)
    value = json.loads(raw, object_pairs_hook=pairs,
                       parse_constant=lambda x: (_ for _ in ()).throw(ValueError("nonfinite JSON")))
    return value, identity, raw


def validate_reference(value):
    require(set(value) == {"schema", "question", "origin", "class", "parameters", "cases", "limits",
                           "first_likelihood_target", "limits_of_claim"},
            "closed scientific reference configuration")
    require(value["schema"] == "lcdm-external-class-reference-design/v1", "reference schema")
    require(value["parameters"] == theory.class_parameters(), "explicit physical parameters differ")
    count = len(value["cases"])
    require(1 <= count <= 4 and json_bytes(value["cases"]) == json_bytes(CASE_DESIGNS[:count]),
            "fixed ordered typed cases")
    require(set(value["class"]) == {"repository", "revision", "tag", "supplied_point", "precision_variant"} and
            value["class"]["revision"] == REVISION and
            value["class"]["supplied_point"]["sha256"] == POINT_SHA and
            value["class"]["precision_variant"]["sha256"] == PRECISION_SHA,
            "supplied CLASS point or precision identity")
    limits = value["limits"]
    require(set(limits) == {"case_wall_seconds", "total_wall_seconds", "address_bytes", "case_cpu_seconds",
                            "attempt_bytes", "combined_logs_bytes", "file_bytes", "max_files", "max_table_rows"},
            "closed resource limits")
    ceilings = {"case_wall_seconds": 180, "total_wall_seconds": 900, "address_bytes": 2147483648,
                "case_cpu_seconds": 180, "attempt_bytes": 2147483648, "combined_logs_bytes": 16777216,
                "file_bytes": 134217728, "max_files": 128, "max_table_rows": 100000}
    require(all(type(limits[k]) is int and 0 < limits[k] <= ceilings[k] for k in ceilings),
            "resource limit domain or ceiling")
    return value


def engine_identities(engine, deadline=None):
    require(set(engine) == {"revision", "source_root", "source_manifest", "binary", "build_receipt"},
            "closed CLASS engine identity")
    require(engine["revision"] == REVISION, "pinned CLASS source revision")
    root = path_without_symlinks(engine["source_root"])
    require(root.is_dir() and not root.is_symlink(), "CLASS source directory")
    identities = []
    for key, limit in (("binary", 67108864), ("build_receipt", 2097152), ("source_manifest", 2097152)):
        _, identity = file_bytes(engine[key]["path"], limit, engine[key])
        identities.append({"role": key, **identity})
    manifest, _, _ = load_json(engine["source_manifest"]["path"], expected=engine["source_manifest"])
    # Parent supplies exact regular tracked source bytes; this is no Git/ABI claim.
    require(manifest["source_revision"] == REVISION and manifest["source_root"] == str(root),
            "source manifest origin/root")
    members = manifest["source_files"]
    require(type(members) is list and 0 < len(members) <= 512, "source manifest cardinality")
    seen, total = set(), 0
    for member in members:
        if deadline is not None:
            require(time.monotonic() < deadline, "source verification deadline")
        relative = Path(member["path"])
        require(not relative.is_absolute() and ".." not in relative.parts and
                member["path"] not in seen, "source manifest relative path")
        seen.add(member["path"])
        pin = {"path": str(root / relative), "bytes": member["bytes"], "sha256": member["sha256"]}
        _, identity = file_bytes(pin["path"], 16777216, pin)
        identities.append({"role": "source/" + member["path"], **identity})
        total += identity["bytes"]
        require(total <= 67108864, "source total byte bound")
    return identities


def inventory(attempt, limits):
    rows, total = [], 0
    for directory, subdirs, files in os.walk(attempt, followlinks=False):
        for name in subdirs:
            require(not (Path(directory) / name).is_symlink(), "output directory symlink")
        for name in files:
            path = Path(directory) / name
            facts = path.lstat()
            require(stat.S_ISREG(facts.st_mode) and facts.st_size <= limits["file_bytes"],
                    "unexpected output type/size")
            require(len(rows) < limits["max_files"] - 1, "output file-count bound with terminal reserve")
            total += facts.st_size
            require(total <= limits["attempt_bytes"] - 2097152, "attempt terminal reserve")
            rows.append({"path": str(path), "bytes": facts.st_size})
    return {"bytes": total, "files": rows, "transient_peak_proven": False}


def controller_sources():
    folder = Path(__file__).absolute().parent
    require(Path(theory.__file__).absolute() == folder / "theory.py", "theory module origin")
    return [file_bytes(folder / name, 262144)[1] for name in ("run.py", "theory.py")]


def seal_inventory(attempt, limits):
    observed = inventory(attempt, limits)
    rows, errors = [], []
    for item in observed["files"]:
        try:
            path = Path(item["path"])
            path.chmod(0o444)
            rows.append(file_bytes(path, limits["file_bytes"])[1])
        except BaseException as exc:
            errors.append({"path": item["path"], **failure(exc, 512)})
    return {"bytes": observed["bytes"], "files": rows, "errors": errors,
            "record_excluded_to_avoid_self_hash": True, "transient_peak_proven": False}


def run_class(argv, source_root, case_dir, attempt, limits, deadline, remaining_logs):
    """One bounded child; retain raw streams, actual wait4 and cleanup failures."""
    def child_limits():
        resource.setrlimit(resource.RLIMIT_AS, (limits["address_bytes"],) * 2)
        resource.setrlimit(resource.RLIMIT_CPU, (limits["case_cpu_seconds"],) * 2)
        resource.setrlimit(resource.RLIMIT_FSIZE, (limits["file_bytes"],) * 2)
        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    env = {"PATH": "/usr/bin:/bin", **{name: "1" for name in
           ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")}}
    row = {"argv": argv, "cwd": source_root, "environment": env, "status": "not_started",
           "pid": None, "returncode": None, "wait4": None, "cleanup": [], "error": None,
           "observed_log_bytes": 0, "retained_log_bytes": 0}
    process, waited = None, False
    streams, selector = {}, None
    start = time.monotonic()
    try:
        for name in ("stdout", "stderr"):
            streams[name] = (case_dir / (name + ".log")).open("xb")
        process = subprocess.Popen(argv, cwd=source_root, env=env, stdin=subprocess.DEVNULL,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                   start_new_session=True, preexec_fn=child_limits)
        row.update(status="started", pid=process.pid)
        selector = selectors.DefaultSelector()
        for name in streams:
            pipe = getattr(process, name)
            os.set_blocking(pipe.fileno(), False)
            selector.register(pipe, selectors.EVENT_READ, name)
        end = min(deadline - 3, start + limits["case_wall_seconds"])
        while not waited or selector.get_map():
            require(time.monotonic() < end, "CLASS child deadline")
            for key, _ in selector.select(0.05):
                raw = os.read(key.fileobj.fileno(), 65536)
                if not raw:
                    selector.unregister(key.fileobj)
                    continue
                row["observed_log_bytes"] += len(raw)
                available = remaining_logs - row["retained_log_bytes"]
                retained = raw[:max(0, available)]
                require(streams[key.data].write(retained) == len(retained), "log write return")
                row["retained_log_bytes"] += len(retained)
                require(len(retained) == len(raw), "combined log cap")
            if not waited:
                pid, status, usage = os.wait4(process.pid, os.WNOHANG)
                if pid:
                    waited = True
                    process.returncode = os.waitstatus_to_exitcode(status)
                    row.update(returncode=process.returncode,
                               wait4={"status": status, "user_seconds": usage.ru_utime,
                                      "system_seconds": usage.ru_stime, "max_rss_kib": usage.ru_maxrss})
            inventory(attempt, limits)
        require(row["returncode"] == 0, "CLASS returned failure")
        row["status"] = "completed"
    except BaseException as exc:
        row["status"], row["error"] = "failed", failure(exc, 1024)
    finally:
        if process is not None:
            # An unreaped leader reserves its PID. After wait4, the numeric group
            # identity may be reused: do not signal or probe it again.
            group_absent = None
            if waited:
                row["cleanup"].append({"operation": "group-cleanup", "status": "not_attempted_leader_reaped"})
            else:
                for signum in (signal.SIGTERM, signal.SIGKILL):
                    try:
                        os.killpg(process.pid, signum)
                        row["cleanup"].append({"signal": signum, "status": "sent"})
                    except ProcessLookupError:
                        group_absent = True
                        row["cleanup"].append({"signal": signum, "status": "absent"})
                        break
                    except OSError as exc:
                        row["status"] = "failed"
                        row["cleanup"].append({"signal": signum, "status": "unavailable", "message": str(exc)[:512]})
                if group_absent is None:
                    try:
                        os.killpg(process.pid, 0)
                        group_absent = False
                    except ProcessLookupError:
                        group_absent = True
                    except OSError as exc:
                        row["cleanup"].append({"operation": "group-probe", "status": "unavailable",
                                               "message": str(exc)[:512]})
            cleanup_end = min(deadline, time.monotonic() + 2)
            while not waited and time.monotonic() < cleanup_end:
                try:
                    pid, status, usage = os.wait4(process.pid, os.WNOHANG)
                    if pid:
                        waited = True
                        process.returncode = os.waitstatus_to_exitcode(status)
                        row.update(returncode=process.returncode,
                                   wait4={"status": status, "user_seconds": usage.ru_utime,
                                          "system_seconds": usage.ru_stime, "max_rss_kib": usage.ru_maxrss})
                except OSError as exc:
                    row["cleanup"].append({"operation": "wait4", "status": "unavailable", "message": str(exc)[:512]})
                    break
                if not waited:
                    time.sleep(0.01)
            row["group_observed_absent"] = group_absent
            row["descendant_cleanup_proven"] = False
            if not waited:
                row["status"] = "failed"
                row["cleanup"].append({"operation": "reap", "status": "unavailable"})
            for name in ("stdout", "stderr"):
                try:
                    getattr(process, name).close()
                except BaseException as exc:
                    row["status"] = "failed"
                    row["cleanup"].append({"operation": name + "-pipe-close", "status": "unavailable",
                                           "message": str(exc)[:512]})
        if selector is not None:
            try:
                selector.close()
            except BaseException as exc:
                row["status"] = "failed"
                row["cleanup"].append({"operation": "selector-close", "status": "unavailable",
                                       "message": str(exc)[:512]})
        for name, stream in streams.items():
            try:
                stream.close()
                (case_dir / (name + ".log")).chmod(0o444)
                _, row[name + "_identity"] = file_bytes(case_dir / (name + ".log"),
                                                       limits["combined_logs_bytes"])
            except BaseException as exc:
                row["status"] = "failed"
                row["cleanup"].append({"operation": name + "-log-seal", "status": "unavailable",
                                       "message": str(exc)[:512]})
    row["wall_seconds"] = time.monotonic() - start
    return row


def products(case_dir, parameters, limits, binary_identity):
    identities = {}
    def table(name):
        raw, pin = file_bytes(case_dir / name, limits["file_bytes"])
        identities[name] = pin
        return theory.parse_table(raw, limits["max_table_rows"])
    cmb = theory.cmb_spectra(table("class_cl_lensed.dat"), float(parameters["T_cmb"]), 2508)
    background = table("class_background.dat")
    thermo = table("class_thermodynamics.dat")
    stdout, identities["stdout"] = file_bytes(case_dir / "stdout.log", 16777216)
    drag = theory.predicted_drag(stdout)
    linear, nonlinear = [], []
    for index, z in enumerate(theory.REDSHIFTS, 1):
        linear.append(theory.matter_spectrum(table(f"class_z{index}_pk.dat"),
                                             float(parameters["H0"]) / 100, z))
        nonlinear.append(theory.matter_spectrum(table(f"class_z{index}_pk_nl.dat"),
                                                float(parameters["H0"]) / 100, z))
    state = {"parameters": parameters, "background": background, "drag": drag,
             "source_identities": {**identities, "binary": binary_identity}}
    result = {"cmb": cmb, "linear_matter": linear, "nonlinear_matter": nonlinear,
              "accuracy_scope": {"table_format": "%.12e: 12 decimal places in mantissa",
                                 "solver_certified_error_bound": None,
                                 "background_interpolation_error_bound": None},
              "nonlinear_method": "halofit", "drag": drag,
              "background": theory.background_at_many(background, theory.REDSHIFTS),
              "thermodynamics": {"columns": thermo["columns"], "rows": len(thermo["rows"]),
                                 "raw_identity": identities["class_thermodynamics.dat"]},
              "source_identities": state["source_identities"]}
    return result, state


def run(config_path, attempt_path, deadline_utc):
    start = time.monotonic()
    config, config_pin, config_raw = load_json(config_path)
    require({"reference", "engine", "likelihood"} <= set(config) <=
            {"reference", "engine", "likelihood", "bao"}, "closed runtime configuration")
    require(config["likelihood"] is None and config.get("bao") is None,
            "first prediction-only controller requires likelihood:null and bao:null")
    reference, reference_pin, reference_raw = load_json(config["reference"]["path"], expected=config["reference"])
    validate_reference(reference)
    limits = reference["limits"]
    utc = datetime.datetime.fromisoformat(deadline_utc)
    require(utc.tzinfo is not None, "root deadline needs UTC offset")
    seconds = (utc - datetime.datetime.now(datetime.timezone.utc)).total_seconds()
    require(seconds > 3, "expired root deadline")
    deadline = min(start + limits["total_wall_seconds"], start + seconds)
    attempt = path_without_symlinks(attempt_path)
    attempt.mkdir(mode=0o700)  # Existing attempts are never reused.
    record = {"schema": "lcdm-external-class-attempt/v1", "status": "failed", "config": config_pin,
              "reference": reference_pin, "engine_before": None, "engine_after": None,
              "controller_before": None, "controller_after": None,
              "cases": [{**x, "status": "not_started", "gates": {"theory": "unassessed",
                         "likelihood": "not_requested" if config["likelihood"] is None else "unassessed"}}
                        for x in reference["cases"]],
              "precision_comparison": None, "failures": [],
              "gates": {"execution": "unassessed", "numerical": "unassessed",
                        "inference": "blocked", "interpretation": "blocked"}}
    results, log_bytes, active_row = {}, 0, None
    try:
        record["controller_before"] = controller_sources()
        write_new(attempt / "config.snapshot.json", config_raw)
        write_new(attempt / "reference.snapshot.json", reference_raw)
        for key in ("supplied_point", "precision_variant"):
            declared = reference["class"][key]
            relative = Path(declared["path"])
            require(not relative.is_absolute() and ".." not in relative.parts, "source example path")
            pin = {"path": str(Path(config["engine"]["source_root"]) / relative),
                   "bytes": declared["bytes"], "sha256": declared["sha256"]}
            raw, _ = file_bytes(pin["path"], 65536, pin)
            write_new(attempt / ("supplied-point.ini" if key == "supplied_point" else "supplied-precision.pre"), raw)
        record["engine_before"] = engine_identities(config["engine"], deadline)
        binary = config["engine"]["binary"]
        for row in record["cases"]:
            active_row = row
            require(time.monotonic() < deadline - 3, "whole prediction deadline")
            row["status"] = "preparing"
            row["gates"]["theory"] = "failed"
            engine_identities(config["engine"], deadline)
            case_dir = attempt / row["id"]
            case_dir.mkdir(mode=0o700)
            parameters = theory.class_parameters(row["ns"])
            ini = case_dir / "input.ini"
            # input_set_root appends '_' even with overwrite_root=yes.
            row["input"] = write_new(ini, theory.render_ini(parameters, case_dir / "class"))
            argv = [binary["path"], str(ini)]
            if row["precision"]:
                argv.append(str(attempt / "supplied-precision.pre"))
            process = run_class(argv, config["engine"]["source_root"], case_dir, attempt, limits,
                                deadline, limits["combined_logs_bytes"] - log_bytes)
            row["process"] = process
            row["status"] = process["status"]
            log_bytes += process["retained_log_bytes"]
            require(process["status"] == "completed", "CLASS case failed; dependent cases withheld")
            row["status"] = "processing_products"
            result, state = products(case_dir, parameters, limits, binary)
            state["source_identities"]["input_ini"] = row["input"]
            result["source_identities"]["input_ini"] = row["input"]
            row["products"] = write_new(case_dir / "products.json", json_bytes(result))
            row["gates"]["theory"] = "passed"
            row["likelihood"] = None
            row["background_state_available"] = True
            row["bao"] = None
            results[row["id"]] = result
            inventory(attempt, limits)
            require(time.monotonic() < deadline, "whole reference deadline")
            row["status"] = "completed"
            active_row = None
        if "anchor" in results and "precision" in results:
            record["precision_comparison"] = theory.product_difference(results["anchor"], results["precision"])
        record["status"] = "completed"
        record["gates"]["execution"] = "passed"
        record["gates"]["numerical"] = ("empirical-differences-recorded-no-certified-bound"
                                          if record["precision_comparison"] is not None else "unassessed")
    except BaseException as exc:
        if active_row is not None:
            active_row["status"] = "failed"
        record["failures"].append(failure(exc))
        record["gates"]["execution"] = "failed"
    finally:
        try:
            record["engine_after"] = engine_identities(config["engine"])
            require(record["engine_before"] == record["engine_after"], "terminal CLASS/source drift")
            file_bytes(config_pin["path"], 2097152, config_pin)
            file_bytes(reference_pin["path"], 2097152, reference_pin)
        except BaseException as exc:
            record["status"] = "failed"
            record["gates"]["execution"] = "failed"
            record["failures"].append(failure(exc))
        try:
            record["controller_after"] = controller_sources()
            require(record["controller_before"] == record["controller_after"], "terminal controller source drift")
        except BaseException as exc:
            record["status"] = "failed"
            record["gates"]["execution"] = "failed"
            record["failures"].append(failure(exc))
        try:
            record["output_inventory"] = seal_inventory(attempt, limits)
            require(not record["output_inventory"]["errors"], "raw output sealing refused")
        except BaseException as exc:
            record["status"] = "failed"
            record["gates"]["execution"] = "failed"
            record["failures"].append(failure(exc))
        if time.monotonic() >= deadline:
            record["status"] = "failed"
            record["gates"]["execution"] = "failed"
            record["failures"].append({"kind": "TimeoutError", "message": "whole controller deadline at final sealing"})
        record["wall_seconds"] = time.monotonic() - start
        record["wall_seconds_scope"] = "start through final source checks and output sealing; excludes terminal record serialization/write"
        record["outer_watchdog_required"] = True
        raw = json_bytes(record)
        require(len(raw) <= 2097152, "terminal record reserved byte bound")
        write_new(attempt / "record.json", raw)
    print(json.dumps({"status": record["status"], "record": str(attempt / "record.json"),
                      "wall_seconds_through_record_write": time.monotonic() - start}, allow_nan=False))
    return 0 if record["status"] == "completed" else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--attempt", required=True)
    parser.add_argument("--deadline-utc", required=True, help="parent job grant's absolute timezone-aware deadline")
    args = parser.parse_args()
    previous = {}
    def interrupted(signum, frame):
        raise InterruptedError(f"controller received signal {signum}")
    for signum in (signal.SIGTERM, signal.SIGINT):
        previous[signum] = signal.getsignal(signum)
        signal.signal(signum, interrupted)
    try:
        return run(args.config, args.attempt, args.deadline_utc)
    finally:
        for signum, handler in previous.items():
            signal.signal(signum, handler)


if __name__ == "__main__":
    raise SystemExit(main())
