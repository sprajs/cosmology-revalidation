#!/usr/bin/env python3
"""Execute one reviewed Irreducible request; retain every admitted attempt locally."""
import argparse
from datetime import datetime, timezone
from importlib.metadata import version
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import uuid

from packet import ROOT, load, read_packet, sha256, verify_inputs, within


def utc():
    return datetime.now(timezone.utc).isoformat()


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def source_identity():
    def git(*args):
        return subprocess.check_output(["git", "-C", str(ROOT), *args], text=True).strip()
    return {"revision": git("rev-parse", "HEAD"),
            "status": git("status", "--porcelain", "--untracked-files=all"),
            "runner_sha256": sha256(__file__),
            "packet_code_sha256": sha256(Path(__file__).with_name("packet.py"))}


def execute(folder, binary, name=None):
    folder = within(ROOT / "experiments", folder)
    packet, request = read_packet(folder)
    if request is None:
        raise ValueError("Experiment is blocked: " + "; ".join(packet["blockers"]))
    verify_inputs(packet)
    resolved = shutil.which(str(binary))
    if resolved is None:
        raise ValueError("Irreducible executable not found; build it separately and pass --irred")
    executable = Path(resolved).resolve()
    executable_hash = sha256(executable)
    discovery = subprocess.run([str(executable), "describe", "--json"],
                               capture_output=True, timeout=30, check=True)
    results = ROOT / "results"
    results.mkdir(exist_ok=True)
    label = name or (datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid.uuid4().hex[:8])
    store = within(within(results, packet['id']), label)
    store.mkdir(parents=True, exist_ok=False)
    record = {"schema_version": 1, "experiment": packet["id"], "started_utc": utc(),
              "execution": "incomplete", "returncode": None,
              "reproducible": source_identity(),
              "python": sys.version, "irred_executable_sha256": executable_hash,
              "jsonschema_version": version("jsonschema"),
              "lock_sha256": sha256(ROOT / "uv.lock") if (ROOT / "uv.lock").is_file() else None,
              "schema_sha256": sha256(ROOT / "schemas/experiment.schema.json"),
              "inputs": packet["inputs"],
              "qualification": {"numerical": "not_assessed", "inference": "not_assessed",
                                "interpretation": "not_assessed"}}
    code = 1
    try:
        write(store / "run.json", record)
        (store / "describe.json").write_bytes(discovery.stdout)
        describe = load(store / "describe.json")
        build = describe.get("build", {})
        record["irred_build_id"] = build.get("build_id")
        record["irred_source_revision"] = build.get("git_head")
        record["irred_source_status"] = build.get("git_status")
        if (describe.get("product") != "Irreducible"
                or build.get("git_head") != packet["execution"]["engine_revision"]
                or build.get("git_status") != ""):
            raise ValueError("Irreducible must be built from the packet's clean, pinned revision")
        if not any(c.get("id") == packet["execution"]["operation"]
                   and c.get("implementation") == "implemented"
                   for c in describe.get("capabilities", [])):
            raise ValueError("Requested operation is absent from this build")
        shutil.copyfile(folder / "experiment.json", store / "experiment.json")
        shutil.copyfile(request, store / "request.json")
        record["packet_sha256"] = sha256(store / "experiment.json")
        record["request_sha256"] = sha256(store / "request.json")
        argv = [str(executable), "run", str(store / "request.json"), str(store / "engine"),
                "--assurance", packet["execution"]["assurance"]]
        record["command"] = argv
        write(store / "run.json", record)
        with (store / "stdout.json").open("wb") as stdout, (store / "stderr.txt").open("wb") as stderr:
            process = subprocess.run(argv, cwd=ROOT, stdout=stdout, stderr=stderr,
                                     timeout=packet["execution"]["timeout_seconds"])
        record["returncode"] = process.returncode
        code = process.returncode if process.returncode >= 0 else 1
        record["execution"] = "completed" if code == 0 else "failed"
        if (store / "stdout.json").stat().st_size:
            reply = load(store / "stdout.json")
            receipt = reply.get("receipt", {})
            record["qualification"] = {k: receipt.get(k, "not_assessed")
                                       for k in ("numerical", "inference", "interpretation")}
            record["engine_accepted"] = receipt.get("accepted")
            record["engine_execution"] = receipt.get("execution", "not_assessed")
            record["engine_accepted_scope"] = receipt.get("accepted_scope")
            if (receipt.get("executable_digest") != executable_hash
                    or receipt.get("build_id") != build["build_id"]
                    or receipt.get("source_revision") != build["git_head"]
                    or receipt.get("source_status") != ""
                    or receipt.get("input_digest") != record["request_sha256"]):
                raise ValueError("Engine receipt identity differs from the admitted execution")
            if code == 0 and receipt.get("accepted") is not True:
                raise ValueError("Exit 0 without an accepted engine receipt")
        elif code == 0:
            raise ValueError("Exit 0 without a receipt")
        verify_inputs(packet)
        if sha256(executable) != executable_hash:
            raise ValueError("Executable changed during execution")
    except subprocess.TimeoutExpired as error:
        record.update(execution="timed_out", error=str(error))
        code = 124
    except KeyboardInterrupt:
        record.update(execution="interrupted", error="Interrupted by user")
        code = 130
    except Exception as error:
        record.update(execution="failed", error=str(error))
        code = 1
    finally:
        record["finished_utc"] = utc()
        record["output_sha256"] = {p.name: sha256(p) for p in store.iterdir()
                                    if p.is_file() and p.name != "run.json"}
        write(store / "run.json", record)
    print(json.dumps({"run": str(store.relative_to(ROOT)), "execution": record["execution"],
                      "qualification": record["qualification"], "error": record.get("error")}))
    return code


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("experiment", help="Directory name below experiments/")
    parser.add_argument("--irred", default=os.environ.get("IRRED_BINARY", "irred"))
    parser.add_argument("--name", help="New relative run name; never overwrite an existing run")
    args = parser.parse_args()
    try:
        return execute(args.experiment, args.irred, args.name)
    except Exception as error:
        print(str(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
