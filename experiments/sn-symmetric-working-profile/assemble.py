"""Exact source assembly and retained CLASS mean mapping; no factor or score."""
import argparse
from array import array
from fractions import Fraction
import hashlib
import importlib.machinery
import math
import os
from pathlib import Path
import struct
import sys
import time
import types
import transport as t

def rat(x):
    return Fraction.from_float(x)

def fraction_record(x):
    return {"numerator": str(x.numerator), "denominator": str(x.denominator)}

def rn_even(x):
    # One rational -> binary64 cast, then verify closest-neighbor/tie parity.
    v = float(x)
    t.need(math.isfinite(v), "assembly overflow")
    distance = abs(rat(v) - x)
    parity = struct.unpack(">Q", struct.pack(">d", v))[0] & 1
    for adjacent in (math.nextafter(v, -math.inf), math.nextafter(v, math.inf)):
        if math.isfinite(adjacent):
            other = abs(rat(adjacent) - x)
            t.need(distance <= other and (distance != other or parity == 0), "rational cast not RN-even")
    return v

def pair_projection(left,right):
    if left==right:
        return left,right,None  # preserve each original bit, including signed zero
    exact=(rat(left)+rat(right))/2
    value=rn_even(exact)
    return value,value,exact

def projection_updates(b,n):
    t.need(len(b)==n*n,"projection matrix shape")
    for i in range(n):
        t.need(b[i*n+i]>0,"positive original diagonal")
        for j in range(i):
            left,right=b[i*n+j],b[j*n+i]
            if left==right:continue
            chosen,mirrored,exact=pair_projection(left,right)
            t.need(chosen.hex()==mirrored.hex(),"identical projected pair bits")
            yield i,j,left,right,chosen,exact

def assembly(request, out, receipt):
    raw, original_authority = t.read(request["ancestry_authority"])
    ancestry = t.document(raw)
    pins = ancestry["files"]
    t.need(len(pins) == 34, "exact original34 input/code roles")
    roles = {pin["role"]: pin for pin in pins}
    t.need(len(roles) == 34, "unique ancestry roles")
    helper_buffers={}
    for pin in pins:
        keep=pin["role"] in ("code/theory.py","code/bao.py","code/sn.py")
        raw, actual = t.read(pin, keep)
        if keep:helper_buffers[pin["role"]]=raw
        receipt["sources_before"].append({"role": pin["role"], **actual})
    selection_raw, selection_pin = t.read(request["selection"])
    selection = t.document(selection_raw)
    original_rows = selection["original_rows"]
    indices = selection["selected_original_indices"]
    t.need(len(original_rows) == 1701 and len(indices) == t.N and indices == sorted(set(indices)), "original/selected order")
    t.need(indices == [r["original_index"] for r in original_rows if r["zHD"] > 0.01 or r["is_calibrator"]], "exact source selection")
    t.need([r["original_index"] for r in original_rows] == list(range(1701)), "original occurrence axis")
    t.need(sum(original_rows[i]["is_calibrator"] for i in indices) == 77 and len(selection["excluded_original_indices"]) == 44,
           "calibrator/complement counts")
    slots = {original: i for i, original in enumerate(indices)}
    covariance_raw, covariance_source = t.read(roles["sn/sn-covariance.cov"])
    position = covariance_raw.find(b"\n") + 1
    t.need(position > 0 and covariance_raw[:position].strip() == b"1701", "original C header")
    b = array("d")
    for i in range(1701):
        for j in range(1701):
            end = covariance_raw.find(b"\n", position)
            if end < 0:
                end = len(covariance_raw)
            raw = covariance_raw[position:end]
            t.need(raw and len(raw)+(end < len(covariance_raw)) <= 128, "C cell bound/truncation")
            value = t.number(raw.decode("ascii").strip())
            if i in slots and j in slots:
                b.append(value)
            position = end + (end < len(covariance_raw))
    t.need(position == len(covariance_raw) and len(b) == t.N*t.N, "full original and selected C shape")
    del covariance_raw
    receipt["original_selected_B"] = t.f64_payload(out / "B.f64be", b)
    shat = array("d", b)
    row_r, row_a, row_lower, row_upper = ([Fraction(0)] * t.N for _ in range(4))
    ledger = []
    for i,j,left,right,chosen,exact in projection_updates(b,t.N):
        shat[i*t.N+j] = shat[j*t.N+i] = chosen
        r, a = rat(chosen)-exact, (rat(left)-rat(right))/2
        e_lower, e_upper = rat(left)-rat(chosen), rat(right)-rat(chosen)
        for row in (i, j):
            row_r[row] += abs(r)
            row_a[row] += abs(a)
            row_lower[row] += abs(e_lower)
            row_upper[row] += abs(e_upper)
        ledger.append({"selected_i": i, "selected_j": j, "original_i": indices[i], "original_j": indices[j],
                           "B_i_j": left.hex(), "B_j_i": right.hex(), "S_hat_both_cells": chosen.hex(),
                           "S0_exact_average": fraction_record(exact), "R_i_j": fraction_record(r), "A_i_j": fraction_record(a),
                           "delta_from_B_i_j": fraction_record(rat(chosen)-rat(left)),
                           "delta_from_B_j_i": fraction_record(rat(chosen)-rat(right))})
    t.need(len(ledger) == 361, "exact361 unequal selected pairs; no source replacement")
    t.need(all(shat[i*t.N+j] == shat[j*t.N+i] for i in range(t.N) for j in range(i)), "derived exact symmetry")
    receipt["S_hat"] = t.f64_payload(out / "S_hat.f64be", shat)
    receipt["assembly_ledger"] = t.emit(out / "assembly-ledger.json", t.encoded({"target_id": t.TARGET,
        "derivation_id": t.DERIVATION, "original_B": receipt["original_selected_B"], "derived_S_hat": receipt["S_hat"],
        "unequal_pairs": ledger, "R_norm_inf_exact": fraction_record(max(row_r)), "A_norm_inf_exact": fraction_record(max(row_a)),
        "lower_completion_minus_S_hat_norm_inf_exact": fraction_record(max(row_lower)),
        "upper_completion_minus_S_hat_norm_inf_exact": fraction_record(max(row_upper)),
        "source_uncertainty": None, "triangle_variants_executed": False}), 1024**2)
    # Preserve exact historical mean mapping via already pinned source, without
    # calling old parse_data, assembling physics, or running CLASS again.
    t.need(not any(name in sys.modules for name in ("sn","bao","theory")),"no preloaded unadmitted mean helper")
    receipt["helper_code_provenance"]=[]
    # Execute only the consumed verified buffers; disk/cache may never supply
    # these three modules. Explicit dependency order resolves their imports.
    for name in ("theory","bao","sn"):
        role="code/"+name+".py";source=roles[role];buffer=helper_buffers[role]
        provenance={"name":name,"filename":source["path"],"bytes":len(buffer),"executed_source_sha256":hashlib.sha256(buffer).hexdigest(),"attempted":True,"completed":False,"method":"consumed-pinned-buffer/builtin-compile-exec/no-bytecode-cache"}
        receipt["helper_code_provenance"].append(provenance)
        t.need(provenance["executed_source_sha256"]==source["sha256"],"consumed helper code identity")
        module=types.ModuleType(name);module.__file__=source["path"];module.__package__=""
        module.__spec__=importlib.machinery.ModuleSpec(name,loader=None,origin=source["path"])
        sys.modules[name]=module
        exec(compile(buffer,source["path"],"exec",dont_inherit=True),module.__dict__)
        provenance["completed"]=True
    sn,theory=sys.modules["sn"],sys.modules["theory"]
    rows = [{**original_rows[i], "id": sn.INPUTS[0]["sha256"] + ":row:" + str(i)} for i in indices]
    data = {"schema": "pantheon-shoes-source-selected-data/v1", "selection": sn.SELECTION, "original_row_count": 1701,
            "rows": rows, "excluded_original_indices": selection["excluded_original_indices"],
            "table_sha256": sn.INPUTS[0]["sha256"], "source_identities": [covariance_source, selection_pin]}
    old_record_raw, old_record_pin = t.read(roles["class_record"])
    reference_raw, reference_pin = t.read(roles["class_reference"])
    old_record, reference = t.document(old_record_raw), t.document(reference_raw)
    t.need(old_record["reference"]["sha256"] == reference_pin["sha256"], "CLASS reference identity")
    receipt["case_vectors"] = []
    background_cache = {}
    all_mu = []
    for index, case_id in enumerate(t.CASES):
        case = old_record["cases"][index]
        t.need(case["id"] == case_id and case["status"] == "completed", "CLASS case order")
        parameters = theory.class_parameters(case["ns"])
        t.need(parameters == {**reference["parameters"], "n_s": case["ns"]}, "same physical point")
        ini, ini_pin = t.read(roles[case_id+"/input"])
        t.need(ini == theory.render_ini(parameters, Path(ini_pin["path"]).parent / "class"), "CLASS INI identity")
        source_ids = {"input_ini": ini_pin, "binary": roles["class_binary"], "retained_record": old_record_pin,
                      "scientific_reference": reference_pin}
        argv = [roles["class_binary"]["path"], ini_pin["path"]]
        t.need(type(case["precision"]) is bool and case["precision"] == reference["cases"][index]["precision"]
               and case["ns"] == reference["cases"][index]["ns"], "precision/ns identity")
        if case["precision"]:
            _, source_ids["precision"] = t.read(roles["class_precision"], False)
            argv.append(roles["class_precision"]["path"])
        t.need(case["process"]["argv"] == argv, "actual CLASS argv")
        raw, source_ids["class_background.dat"] = t.read(roles[case_id+"/background"])
        digest = source_ids["class_background.dat"]["sha256"]
        if digest not in background_cache:
            background_cache[digest] = theory.parse_table(raw, max_rows=100000)
        del raw
        _, source_ids["class_thermodynamics.dat"] = t.read(roles[case_id+"/thermodynamics"], False)
        log, source_ids["stdout"] = t.read(roles[case_id+"/stdout"])
        prediction = sn.predict({"background": background_cache[digest], "parameters": parameters,
                                 "drag": theory.predicted_drag(log), "source_identities": source_ids}, data, case_id)
        prediction["working_covariance_target_id"] = t.TARGET
        prediction["working_covariance_sha256"] = receipt["S_hat"]["sha256"]
        pin = t.emit(out / (case_id+".prediction.json"), t.encoded(prediction), 4*1024**2)
        mu = [row["distance_modulus_without_M"] for row in prediction["rows"]]
        residual = [row["observed"]-value for row, value in zip(rows, mu)]
        t.need(all(math.isfinite(x) for x in residual), "binary64 residual finite")
        receipt["case_vectors"].append({"id": case_id, "state_sha256": prediction["common_state_sha256"], "prediction": pin,
                                       "mu_hex": [v.hex() for v in mu], "residual_hex": [v.hex() for v in residual]})
        all_mu.append(mu)
    t.need(len(background_cache) == 1 and all([v.hex() for v in mu] == [v.hex() for v in all_mu[0]] for mu in all_mu),
           "same retained background requires identical four mean vectors")
    receipt["same_background_four_means_bitwise_equal"] = True
    wire = out / "native.input"
    with wire.open("xb") as f:
        count = 0
        def put(text):
            nonlocal count
            raw = text.encode("ascii")
            count += len(raw)
            t.need(count <= t.MAX_FILE and f.write(raw) == len(raw), "wire bound/short write")
        put("SN_SHAT_COMMON_M_V1\n1657 4\n" + sn.INPUTS[0]["sha256"]+" "+sn.INPUTS[1]["sha256"]+" "+receipt["S_hat"]["sha256"]+"\n")
        put(" ".join(str(i) for i in indices)+"\n")
        put(" ".join(repr(r["observed"]) for r in rows)+"\n")
        for i in range(t.N):
            put(" ".join(repr(x) for x in shat[i*t.N:(i+1)*t.N])+"\n")
        for vector, mu in zip(receipt["case_vectors"], all_mu):
            put(vector["id"]+" "+vector["state_sha256"]+"\n")
            put(" ".join(repr(v) for v in mu)+"\n")
        f.flush()
        os.fsync(f.fileno())
    wire.chmod(0o444)
    receipt["native_input"] = t.identity(wire)
    receipt["selection"] = selection_pin
    receipt["status"] = "assembled-new-S_hat-and-four-mean-vectors-no-score"
    receipt["sources_after"]=[]
    for pin in pins:receipt["sources_after"].append({"role":pin["role"],**t.read(pin,False)[1]})
    t.need(receipt["sources_before"] == receipt["sources_after"], "terminal ancestry drift")
    _, after_authority = t.read(request["ancestry_authority"], False)
    t.need(after_authority == original_authority, "terminal original authority drift")
    t.need(t.read(request["selection"], False)[1] == selection_pin, "terminal selection drift")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--request", required=True); p.add_argument("--request-sha", required=True); p.add_argument("--output", required=True)
    args = p.parse_args()
    request = None
    out = Path(args.output)
    t.need(out.is_absolute() and not out.exists(), "fresh absolute assembly directory")
    out.mkdir(mode=0o700)
    receipt = {"target_id": t.TARGET, "derivation_id": t.DERIVATION, "request":None, "status": "refused",
               "sources_before": [], "sources_after": [], "normalized_density": None, "joint_target": None}
    start = time.monotonic()
    try:
        request,receipt["request"] = t.admitted_request(args.request, args.request_sha)
        assembly(request, out, receipt)
    except BaseException as exc:
        receipt["status"] = "refused"
        receipt["error"] = t.failure(exc)
        # Retry only terminal identity reads, never assembly/science.
        if request is not None:
            try:
                raw, _ = t.read(request["ancestry_authority"])
                receipt["sources_after"]=[]
                for pin in t.document(raw)["files"]:receipt["sources_after"].append({"role":pin["role"],**t.read(pin,False)[1]})
            except BaseException as terminal:
                receipt["terminal_identity_error"] = t.failure(terminal)
    finally:
        receipt["wall_seconds"] = time.monotonic()-start
        receipt["retained_files"] = [t.identity(p) for p in sorted(out.iterdir()) if p.is_file()]
        t.store_bytes(out)
        t.emit(out / "assembly.json", t.encoded(receipt), 2*1024**2)
        for path in out.iterdir():
            path.chmod(0o444)
        out.chmod(0o555)
    return 0 if receipt["status"] == "assembled-new-S_hat-and-four-mean-vectors-no-score" else 2

if __name__ == "__main__":
    sys.exit(main())
