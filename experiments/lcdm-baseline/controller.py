#!/usr/bin/env python3
"""One bounded LCDM native experiment; no configurable commands or physics."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import resource
import subprocess
import sys
import uuid

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from packet import load, parse, read_packet, sha256, verify_inputs, within
import reference

FOLDER = Path(__file__).resolve().parent
OP = {"DM_over_rs": 0, "DH_over_rs": 1, "DV_over_rs": 2}
PRODUCER = {"distance_absolute_mpc":1e-12,"distance_relative":2e-14,"ratio_absolute":1e-14,"ratio_relative":2e-13,"sound_absolute_mpc":1e-12,"sound_relative":2e-14}
BUDGETS = {"distance_absolute_mpc": 1e-9, "distance_relative": 2e-11,
           "ratio_absolute": 1e-11, "ratio_relative": 5e-11,
           "density_absolute": 1e-8, "reference_fraction": .05}


def keys(value, expected):
    if not isinstance(value, dict) or set(value) != set(expected):
        raise ValueError("Unexpected typed fields")


def finite(value):
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError("Expected finite typed number")
    return value


def validate_request(q):
    keys(q, ("schema_version", "operation", "model", "redshifts", "rows", "budgets", "engine_identity", "resources"))
    if q["schema_version"] != 2 or q["operation"] != "lcdm-baseline.native":
        raise ValueError("Wrong native experiment")
    keys(q["model"], ("h0_km_s_mpc", "omega_m", "omega_r", "omega_b", "omega_gamma", "z_drag", "drag_origin"))
    m = q["model"]
    for k, v in m.items():
        if k != "drag_origin":
            finite(v)
    if (m["h0_km_s_mpc"] <= 0 or m["omega_r"] <= 0 or m["omega_gamma"] <= 0
            or not 0 <= m["omega_b"] <= m["omega_m"] or m["omega_gamma"] > m["omega_r"]
            or m["omega_m"]+m["omega_r"] > 1 or m["z_drag"] < 0):
        raise ValueError("Unsupported native physical domain")
    origin = m["drag_origin"]
    if not isinstance(origin, str) or not origin or len(origin.encode()) > 4096 or any(ord(c)<32 or ord(c)>126 for c in origin):
        raise ValueError("Invalid supplied drag origin")
    if not isinstance(q["redshifts"], list) or not 1 <= len(q["redshifts"]) <= 16:
        raise ValueError("Background row quota")
    for z in q["redshifts"]:
        if finite(z) < 0 or z > 5:
            raise ValueError("Background redshift domain")
    if not isinstance(q["rows"], list) or len(q["rows"]) != 13:
        raise ValueError("Exactly 13 ordered source rows required")
    for i, row in enumerate(q["rows"]):
        keys(row, ("id", "z", "observable"))
        if row["observable"] not in OP or row["id"] != f'desi-dr2-{i:02d}-{row["observable"]}' or not 0 <= finite(row["z"]) <= 5:
            raise ValueError("Wrong typed row axis")
    if q["budgets"] != BUDGETS:
        raise ValueError("Frozen comparison allocation changed")
    keys(q["resources"], ("maximum_input_bytes", "maximum_output_bytes", "timeout_seconds", "reference_panels"))
    if q["resources"] != {"maximum_input_bytes": 65536, "maximum_output_bytes": 65536, "timeout_seconds": 120, "reference_panels": [128, 256]}:
        raise ValueError("Frozen bounded resource policy changed")
    keys(q["engine_identity"], ("revision", "build_id", "manifest_sha256", "archive_sha256", "cli_sha256"))
    for k, v in q["engine_identity"].items():
        size = 40 if k == "revision" else 64
        if not isinstance(v, str) or len(v) != size or any(c not in "0123456789abcdef" for c in v):
            raise ValueError("Malformed pinned engine identity")


def validate_admission(packet,q):
    execution=packet["execution"]
    if (execution["engine_revision"]!=q["engine_identity"]["revision"]
            or execution["operation"]!=q["operation"]
            or execution.get("interface")!="native_lcdm_baseline"
            or execution["assurance"]!="numerical_contract"
            or execution["timeout_seconds"]!=q["resources"]["timeout_seconds"]):
        raise ValueError("Packet/request engine/interface/assurance/timeout identity differs")


def scientific_inputs(packet, q):
    verify_inputs(packet)
    expected_paths=["data/bao/desi_gaussian_bao_ALL_GCcomb_mean.txt","data/bao/desi_gaussian_bao_ALL_GCcomb_cov.txt"]
    if [x["path"] for x in packet["inputs"]]!=expected_paths or any(x["role"] != "released_fitted_summary" for x in packet["inputs"]):
        raise ValueError("Expected exact unique ordered mean/covariance released compression inputs")
    manifest = load(ROOT / "sources/legacy-inputs.json")
    historical = {x["path"]: x for x in manifest["files"] if x.get("group") == "bao"}
    for x in packet["inputs"]:
        if x["path"] not in historical or any(x[k] != historical[x["path"]][k] for k in ("bytes", "sha256")):
            raise ValueError("Source manifest identity differs")
    mean = within(ROOT, "data/bao/desi_gaussian_bao_ALL_GCcomb_mean.txt")
    covariance = within(ROOT, "data/bao/desi_gaussian_bao_ALL_GCcomb_cov.txt")
    source_specs={x["path"]:x for x in packet["inputs"]}
    def admitted_bytes(path):
        item=source_specs[str(path.relative_to(ROOT))]
        if path.stat().st_size!=item["bytes"] or item["bytes"]>65536:
            raise ValueError("Consumed input byte count differs")
        payload=path.read_bytes()
        if len(payload)!=item["bytes"] or hashlib.sha256(payload).hexdigest()!=item["sha256"]:
            raise ValueError("Consumed input bytes/hash differ from admitted source")
        return payload.decode("utf-8")
    mean_text,covariance_text=admitted_bytes(mean),admitted_bytes(covariance)
    observed, axes = [], []
    for line in mean_text.splitlines():
        if line.startswith("#"):
            continue
        z, y, kind = line.split()
        axes.append((float(z), kind))
        observed.append(finite(float(y)))
    if axes != [(x["z"], x["observable"]) for x in q["rows"]]:
        raise ValueError("Exact mean axes/order differs from reviewed full covariance")
    matrix = [[finite(float(v)) for v in line.split()] for line in covariance_text.splitlines() if line.strip() and not line.startswith("#")]
    if len(matrix) != 13 or any(len(row) != 13 for row in matrix) or any(matrix[i][j] != matrix[j][i] for i in range(13) for j in range(13)):
        raise ValueError("Exact full covariance shape/symmetry differs")
    return observed, matrix


def git(source, *args):
    return subprocess.check_output(["git", "-C", str(source), *args], text=True, timeout=30).strip()


def verify_headers(sdk,build):
    expected_headers={Path(name).name for name in build["sources"] if name.startswith("cpp/include/irred/")}
    actual_headers={p.name for p in (sdk/"include/irred").iterdir() if p.is_file()}
    if expected_headers!=actual_headers:
        raise ValueError("Installed SDK header inventory differs")
    for header in (sdk / "include/irred").iterdir():
        if not header.is_file() or sha256(header) != build["sources"].get("cpp/include/irred/"+header.name):
            raise ValueError("Installed SDK header differs")

def fingerprint(source, sdk, q):
    expected = q["engine_identity"]
    if git(source, "rev-parse", "HEAD") != expected["revision"] or git(source, "status", "--porcelain", "--untracked-files=all"):
        raise ValueError("Engine source must be clean at pinned revision")
    manifest = sdk / "build-manifest.json"
    if sha256(manifest) != expected["manifest_sha256"]:
        raise ValueError("SDK manifest hash differs")
    build = load(manifest)
    if build["git_head"] != expected["revision"] or build["git_status"] != "" or build["build_id"] != expected["build_id"]:
        raise ValueError("SDK clean build identity differs")
    content = dict(build)
    for k in ("build_id", "git_head", "git_status"):
        content.pop(k)
    if hashlib.sha256(json.dumps(content, sort_keys=True, separators=(",", ":")).encode()).hexdigest() != build["build_id"]:
        raise ValueError("Build manifest identity does not reconstruct")
    identity = {"reference_python_executable":sha256(Path(sys.executable))}
    for name, digest in build["sources"].items():
        path = within(source, name)
        if sha256(path) != digest:
            raise ValueError("Engine source hash differs: " + name)
        identity["source/"+name] = digest
    verify_headers(sdk,build)
    for relative in ("lib/libirred_core.a", "bin/irred"):
        path = sdk / relative
        digest = expected["archive_sha256" if relative.startswith("lib") else "cli_sha256"]
        if sha256(path) != digest:
            raise ValueError("SDK archive/executable hash differs")
    for path, digest in ((Path("/usr/bin/c++"), build["compiler_executable_digest"]),
                         (Path(build["standard_library"]), build["standard_library_digest"])):
        if sha256(path) != digest:
            raise ValueError("Release tool/archive executable differs")
        identity[str(path)] = digest
    for path in sorted(sdk.rglob("*")):
        if path.is_file():
            identity["sdk/"+str(path.relative_to(sdk))] = sha256(path)
    for path in sorted(FOLDER.iterdir()):
        if path.is_file():
            identity["packet/"+path.name] = sha256(path)
    for path in (ROOT / "scripts/packet.py", ROOT / "schemas/experiment.schema.json", ROOT / "sources/legacy-inputs.json", ROOT / "uv.lock"):
        identity[str(path.relative_to(ROOT))] = sha256(path)
    return identity, build


def bounded(command, store, label, timeout, limit=65536, stdin=None, expected_status=0):
    def limits():
        resource.setrlimit(resource.RLIMIT_FSIZE, (limit, limit))
        resource.setrlimit(resource.RLIMIT_AS, (1073741824, 1073741824))
    with (store / (label+".out")).open("xb") as out, (store / (label+".err")).open("xb") as err:
        process = subprocess.run(command, input=stdin, stdout=out, stderr=err, cwd=ROOT,
                                 timeout=timeout, preexec_fn=limits)
    if process.returncode != expected_status:
        raise ValueError(f"{label} failed with status {process.returncode}")
    return (store / (label+".out")).read_text()


def typed_transport(q, observed, covariance):
    m = q["model"]
    lines = [" ".join(str(m[k]) for k in ("h0_km_s_mpc", "omega_m", "omega_r", "omega_b", "omega_gamma", "z_drag"))+" "+json.dumps(m["drag_origin"]),
             f'13 {len(q["redshifts"])}', " ".join(map(str, q["redshifts"]))]
    lines.extend(f'{json.dumps(row["id"])} {row["z"]} {OP[row["observable"]]} {y}' for row, y in zip(q["rows"], observed))
    lines.extend(" ".join(map(str, row)) for row in covariance)
    payload = ("\n".join(lines)+"\n").encode()
    if len(payload) > q["resources"]["maximum_input_bytes"]:
        raise ValueError("Native payload quota exceeded")
    return payload


def check_outputs(output, q):
    keys(output, ("schema_version", "producer_policy", "density", "background", "callbacks"))
    if output["producer_policy"]!=PRODUCER:
        raise ValueError("Native compiled producer policy differs")
    if output["schema_version"] != 1 or type(output["callbacks"]) is not int or not 0 <= output["callbacks"] <= 4000000:
        raise ValueError("Native schema/callback status differs")
    if len(output["density"]) != 2 or len(output["background"]) != 2:
        raise ValueError("Native model count differs")
    for row in output["density"]:
        keys(row, ("status", "predictions", "quadratic", "log_determinant", "normalization", "log_density", "projection_estimate", "callbacks"))
        if row["status"] != "ok" or len(row["predictions"]) != 13 or not 0 <= finite(row["projection_estimate"]) <= 1e-8:
            raise ValueError("Native density/prediction status failed")
        for value in row["predictions"]+[row[k] for k in ("quadratic", "log_determinant", "normalization", "log_density")]:
            finite(value)
    for model in output["background"]:
        if len(model) != len(q["redshifts"]):
            raise ValueError("Native background length differs")
        for row, z in zip(model, q["redshifts"]):
            keys(row, ("status", "z", "E", "DM", "DL"))
            if row["status"] != "ok" or row["z"] != z:
                raise ValueError("Native background status/axis failed")
            for name in ("E", "DM", "DL"):
                finite(row[name])


def validate_reference(value,q):
    from decimal import Decimal as D
    keys(value,("background","predictions","quadratic","log_determinant","normalization","log_density"))
    if not isinstance(value["background"],list) or len(value["background"])!=len(q["redshifts"]) or not isinstance(value["predictions"],list) or len(value["predictions"])!=13:
        raise ValueError("Independent reference exact vector lengths differ")
    def scalar(x):
        if not isinstance(x,D) or not x.is_finite():
            raise ValueError("Independent reference requires finite Decimal fields")
    for row,z in zip(value["background"],q["redshifts"]):
        keys(row,("z","E","DM","DL"))
        for x in row.values():scalar(x)
        if row["z"]!=D.from_float(float(z)):
            raise ValueError("Independent reference redshift axis differs")
    for x in value["predictions"]+[value[k] for k in ("quadratic","log_determinant","normalization","log_density")]:scalar(x)


def comparisons(output, coarse, fine, q):
    validate_reference(coarse,q)
    validate_reference(fine,q)
    checks = []
    def compare(name, actual, expected, budget, refinement=None):
        from decimal import Decimal as D
        actual = D.from_float(float(actual))
        budget = D(str(budget))
        error = abs(actual-expected)
        if refinement is not None and refinement > budget*D(".05"):
            raise ValueError("Independent reference refinement allocation failed: " + name)
        if error > budget:
            raise ValueError("Scientific comparison failed: " + name)
        checks.append({"name": name, "error": str(error), "budget": str(budget),
                       "reference_refinement": str(refinement) if refinement is not None else None})
    for i, (c, f) in enumerate(zip(coarse["background"], fine["background"])):
        for name in ("E", "DM", "DL"):
            budget = (2e-11*abs(float(f[name]))) if name == "E" else 1e-9+2e-11*abs(float(f[name]))
            compare(f"background[{i}].{name}", output["background"][0][i][name], f[name], budget, abs(c[name]-f[name]))
    for i, (c, f) in enumerate(zip(coarse["predictions"], fine["predictions"])):
        compare(f"BAO[{i}]", output["density"][0]["predictions"][i], f, 1e-11+5e-11*abs(float(f)), abs(c-f))
    for name in ("quadratic", "log_determinant", "normalization", "log_density"):
        compare(name, output["density"][0][name], fine[name], 1e-8, abs(coarse[name]-fine[name]))
    from decimal import Decimal as D
    for i, (a, b) in enumerate(zip(output["density"][0]["predictions"], output["density"][1]["predictions"])):
        compare(f"H0_ratio[{i}]", b, D.from_float(a), 1e-11+5e-11*abs(a))
    for name in ("quadratic", "log_determinant", "normalization", "log_density"):
        compare("H0_"+name, output["density"][1][name], D.from_float(output["density"][0][name]), 1e-8)
    for i, (a, b) in enumerate(zip(*output["background"])):
        for name in ("E", "DM", "DL"):
            expected = D.from_float(a[name])/(1 if name == "E" else 2)
            compare(f"H0_background[{i}].{name}", b[name], expected, 1e-9+2e-11*abs(float(expected)))
    if len(checks)!=6*len(q["redshifts"])+34:
        raise ValueError("Incomplete named comparison coverage")
    return checks


def execute(source, sdk, name=None):
    label = name or datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid.uuid4().hex[:8]
    store = within(ROOT / "results/lcdm-baseline", label)
    store.mkdir(parents=True, exist_ok=False)
    record = {"schema_version": 1, "started_utc": datetime.now(timezone.utc).isoformat(),
              "execution": "failed", "qualification": {"numerical": "not_assessed", "inference": "not_assessed", "interpretation": "unqualified"},
              "python": sys.version,
              "reproducible_revision": git(ROOT, "rev-parse", "HEAD"),
              "reproducible_status": git(ROOT, "status", "--porcelain", "--untracked-files=all"),
              "blocked_full_planck": "0.06 eV massive-neutrino evolution, predicted drag/thermal history and CMB perturbations unavailable"}
    def write():
        (store / "run.json").write_text(json.dumps(record, indent=2, allow_nan=False)+"\n")
    write()
    try:
        for name in ("experiment.json","request.json","candidate.json"):
            path=FOLDER/name
            if path.exists() and (not path.is_file() or path.stat().st_size>65536):
                raise ValueError("Public packet input quota exceeded")
        packet, path, admitted = read_packet(FOLDER)
        if packet["execution"].get("interface") != "native_lcdm_baseline" or path is None:
            raise ValueError("Wrong native packet admission")
        request_bytes=path.read_bytes()
        if hashlib.sha256(request_bytes).hexdigest()!=admitted["request"] or sha256(FOLDER/"experiment.json")!=admitted["packet"]:
            raise ValueError("Packet/request changed after admission")
        q = parse(request_bytes.decode("utf-8"))
        validate_request(q)
        validate_admission(packet,q)
        observed, covariance = scientific_inputs(packet, q)
        before, build = fingerprint(source, sdk, q)
        record.update(admitted=admitted, identities_before=before, inputs=packet["inputs"], engine_build_id=build["build_id"], model=q["model"])
        for filename in ("experiment.json", "request.json", "candidate.json"):
            p = FOLDER / filename
            if p.is_file():
                (store / filename).write_bytes(p.read_bytes())
        if sha256(store/"experiment.json")!=admitted["packet"] or sha256(store/"request.json")!=admitted["request"]:
            raise ValueError("Copied packet/request changed after admission")
        if packet["origin"]["kind"]=="prospector_candidate" and sha256(store/"candidate.json")!=packet["origin"]["sha256"]:
            raise ValueError("Copied candidate changed after admission")
        description = parse(bounded([str(sdk / "bin/irred"), "describe", "--json"], store, "describe", 30, 1048576))
        if description.get("product") != "Irreducible" or description["build"]["build_id"] != build["build_id"] or description["build"]["git_head"] != q["engine_identity"]["revision"] or description["build"]["git_status"] != "":
            raise ValueError("Actual CLI discovery build differs")
        executable = store / "consumer"
        command = ["/usr/bin/c++", "-std=c++20", "-O2", "-Wall", "-Wextra", "-Wpedantic", "-fno-fast-math", "-ffp-contract=off", str(FOLDER / "consumer.cpp"), "-I", str(sdk / "include"), str(sdk / "lib/libirred_core.a"), "-o", str(executable)]
        record["compiler_command"] = command
        bounded(command, store, "compile", 120, 2097152)
        record["consumer_executable_sha256"] = sha256(executable)
        payload = typed_transport(q, observed, covariance)
        (store / "native-input.txt").write_bytes(payload)
        bounded([str(executable),"--underbudget-producer-control"],store,"underbudget-producer",120,stdin=payload,expected_status=2)
        if "prediction 8" not in (store/"underbudget-producer.err").read_text():
            raise ValueError("Under-budget control failed for an unrelated cause")
        record["underbudget_producer_control"]={"status":"conditioning_budget_exceeded","original_ratio_relative":2e-14,"final_ratio_relative":2e-13}
        output = parse(bounded([str(executable)], store, "native", 120, stdin=payload))
        record["actual_compiled_producer_policy"]=output["producer_policy"]
        check_outputs(output, q)
        reference_input=store/"reference-input.json"
        reference_input.write_text(json.dumps({"request":q,"observed":observed,"covariance":covariance},allow_nan=False)+"\n")
        refs=reference.deserialize(parse(bounded([sys.executable,str(FOLDER/"reference.py"),str(reference_input)],store,"reference",120)))
        coarse,fine=refs["coarse"],refs["fine"]
        checks = comparisons(output, coarse, fine, q)
        (store / "checks.json").write_text(json.dumps(checks, indent=2)+"\n")
        after, _ = fingerprint(source, sdk, q)
        scientific_inputs(packet, q)
        if before != after or sha256(executable) != record["consumer_executable_sha256"]:
            raise ValueError("Source/SDK/request/origin/reference/consumer changed during attempt")
        record["identities_after"] = after
        record.update(execution="completed", checks=len(checks), qualification={"numerical": "named_comparisons_passed", "inference": "not_assessed", "interpretation": "conditional_massless_variant_only"})
    except Exception as error:
        record["error"] = str(error)
    finally:
        if "before" in locals():
            try:
                final_identity,_=fingerprint(source,sdk,q)
                scientific_inputs(packet,q)
                record["identities_after"]=final_identity
                if final_identity!=before:
                    raise ValueError("Inputs/source/SDK changed during failed or completed attempt")
            except Exception as error:
                record["execution"]="failed"
                record["error"]=record.get("error", "")+"; after-check: "+str(error)
                record["qualification"]["numerical"]="not_accepted"
        record["finished_utc"] = datetime.now(timezone.utc).isoformat()
        record["output_sha256"] = {p.name: sha256(p) for p in store.iterdir() if p.is_file() and p.name != "run.json"}
        write()
        for p in store.iterdir():
            if p.is_file():
                p.chmod(0o555 if p.name == "consumer" else 0o444)
    print(json.dumps({"run": str(store.relative_to(ROOT)), "execution": record["execution"], "error": record.get("error")}))
    return 0 if record["execution"] == "completed" else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine-source", type=Path, required=True)
    parser.add_argument("--sdk", type=Path, required=True)
    parser.add_argument("--name")
    args = parser.parse_args()
    sys.exit(execute(args.engine_source.resolve(), args.sdk.resolve(), args.name))
