"""One finite serial SN S_hat experiment, pending source/resource admission."""
import argparse
from array import array
from decimal import Decimal, localcontext
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time
import transport as t
import controls

THREADS = ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","BLIS_NUM_THREADS","NUMEXPR_NUM_THREADS","VECLIB_MAXIMUM_THREADS")

def parsed(path,cap=4*1024**2):
    raw,_ = t.owned_json_bytes(path,cap)
    return json.loads(raw,object_pairs_hook=t.pairs,parse_float=Decimal,parse_constant=lambda _: (_ for _ in ()).throw(ValueError("nonfinite JSON")))

def finite_decimal(value):
    t.need(type(value) in (int,Decimal,str) and type(value) is not bool,"typed decimal")
    x = Decimal(value);t.need(x.is_finite(),"finite emitted decimal")
    return x

def native_shape(value):
    top={"interface","target_id","control","original_B_refusal_only","status","stage","sdk_build_id","n","selected_covariance_sha256","source_indices","input_bytes","maximum_forward_sensitivity","arithmetic_policy","work","preparation","profile_preparation","rows","exception","normalized_density","prior","joint_target"}
    t.need(type(value) is dict and set(value)==top,"closed native record")
    t.need(value["status"] in ("accepted","refused") and all(type(value[k]) is bool for k in ("exception","original_B_refusal_only")),"native status/bools")
    t.need(type(value["n"]) is int and 0<=value["n"]<=1701 and type(value["input_bytes"]) is int and 0<=value["input_bytes"]<=t.MAX_FILE,"native typed work/domain")
    t.need(type(value["source_indices"]) is list and len(value["source_indices"])<=t.N and all(type(x) is int and 0<=x<1701 for x in value["source_indices"]),"native index axis")
    w=value["work"];t.need(set(w)=={"gaussian_prepare_attempts","profile_prepare_attempts","evaluate_attempts","condition_inverse_column_upper_bound"} and all(type(x) is int and x>=0 for x in w.values()) and w["gaussian_prepare_attempts"]<=1 and w["profile_prepare_attempts"]<=1 and w["evaluate_attempts"]<=4 and w["condition_inverse_column_upper_bound"]<=t.N,"native bounded invocation counts")
    t.need(type(value["rows"]) is list and len(value["rows"])==4,"four native case slots")
    if value["profile_preparation"] is not None:
        p=value["profile_preparation"]
        t.need(set(p)=={"status","numerical_status","payload_available","cached_response_backward_residual","cached_response_forward_sensitivity","gram","cached_response_solution","adjusted_solution_vector","native_stationarity"} and type(p["payload_available"]) is bool,"closed prepared response availability")
        if not p["payload_available"]:
            t.need(all(p[k] is None for k in ("cached_response_backward_residual","cached_response_forward_sensitivity","gram","cached_response_solution")),"unearned response diagnostics withheld")
    for case,row in zip(t.CASES,value["rows"]):
        t.need(set(row)=={"case_id","common_state_sha256","input_complete","attempted","accepted","status","numerical_status","payload"} and row["case_id"]==case and all(type(row[k]) is bool for k in ("input_complete","attempted","accepted")),"closed native row/flags")
        t.need(not row["accepted"] or (row["attempted"] and row["input_complete"] and row["payload"] is not None),"native causal scalar availability")
        if row["payload"] is not None:
            p=row["payload"]
            t.need(set(p)=={"M_coefficient","quadratic","relative_profile_score","adjusted_residuals","backward_residual","estimated_forward_sensitivity","coefficient_solve_backward_residual","coefficient_solve_forward_sensitivity","residual_l1","solution_norm_inf","adjusted_residual_l1","adjusted_solution_norm_inf"},"closed native profile payload")
            t.need(type(p["adjusted_residuals"]) is list and len(p["adjusted_residuals"])==t.N,"native full adjusted vector")
            for k,x in p.items():
                if k!="adjusted_residuals":finite_decimal(x)
    t.need(value["normalized_density"] is value["prior"] is value["joint_target"] is None,"relative native target only")
    return value

def reference_shape(value):
    top={"interface","target_id","control","status","request","assembly","runtime_before","runtime_after","eigenvalues","eigenvectors","response_solution","diagnostics","rows","derived_target_sensitivity","error","terminal_errors","work","normalized_density","prior","source_uncertainty","joint_target","error_certificate","wall_seconds"}
    t.need(type(value) is dict and set(value)==top and value["status"] in ("accepted","refused","runtime-only"),"closed reference record/status")
    t.need(type(value["rows"]) is list and len(value["rows"])<=4,"reference earned prefix")
    for case,row in zip(t.CASES,value["rows"]):
        t.need(set(row)=={"case_id","M_coefficient","quadratic","relative_profile_score","eta_adjusted","scaled_gate","stationarity_signed","stationarity","adjusted_residuals","accepted","adjusted_solution_reportcast"} and row["case_id"]==case and type(row["accepted"]) is bool,"closed reference case")
        t.need(type(row["adjusted_residuals"]) is list and len(row["adjusted_residuals"])==t.N,"reference no-drop")
        for key in ("M_coefficient","quadratic","relative_profile_score","eta_adjusted","scaled_gate","stationarity_signed","stationarity"):
            t.need(type(row[key]) is str and len(row[key])<=64,"reference decimal-string representation")
            finite_decimal(row[key])
    t.need(value["normalized_density"] is value["prior"] is value["source_uncertainty"] is value["joint_target"] is value["error_certificate"] is None,"reference target/qualification scope")
    return value

def pin_runtime(record):
    record=t.unpack_runtime(record)
    t.need(record["schema"]=="SN-spectral-runtime/v1" and record["numpy_version"]=="2.2.6","spectral runtime schema/version")
    inventory=record["inventory"]
    t.need(len(inventory)<=8192 and [p["path"] for p in inventory]==sorted({p["path"] for p in inventory}),"canonical runtime inventory")
    t.need(hashlib.sha256(t.encoded(inventory)).hexdigest()==record["inventory_sha256"],"runtime inventory digest")
    for pin in inventory:t.need(t.read(pin,False)[1]==pin,"independent runtime file identity")
    t.need(record["thread_environment"]=={k:"1" for k in THREADS},"runtime thread environment")
    for pool in record["numeric_configuration"]["threadpools"]:
        t.need(type(pool["num_threads"]) is int and pool["num_threads"]==1 and str(Path(pool["filepath"]).resolve()) in {p["path"] for p in inventory},"actual backend one-thread identity")
    t.need(len(record["imported_modules"])<=4096 and len(record["mapped_libraries"])<=256,"actual runtime set bounds")
    return True

def native_admission(record,assembly):
    native_shape(record)
    t.need(record["interface"]=="SN-S_hat-common-M-native/v1" and record["target_id"]==t.TARGET and record["status"]=="accepted","native working target status")
    t.need(record["control"]=="none" and record["original_B_refusal_only"] is False and record["n"]==t.N,"actual native domain")
    t.need(record["sdk_build_id"]=="6d495efb166006c6ce651359a366af5d87686ce516ecd7f0f49eed4e8d8be07e","actual native SDK")
    t.need(record["selected_covariance_sha256"]==assembly["S_hat"]["sha256"],"same stored S_hat")
    selection = t.document(t.read(assembly["selection"])[0])
    t.need(record["source_indices"]==selection["selected_original_indices"],"every original selected occurrence")
    t.need(record["work"]["gaussian_prepare_attempts"]==record["work"]["profile_prepare_attempts"]==1 and record["work"]["evaluate_attempts"]==4,"native exact preparation/evaluation counts")
    prep=record["profile_preparation"]
    t.need(prep["status"]=="finite" and prep["numerical_status"]=="ok" and prep["payload_available"] is True,"cached response status/availability")
    t.need(finite_decimal(prep["gram"])>0 and 0<=finite_decimal(prep["cached_response_forward_sensitivity"])<=Decimal("1e-8"),"cached response sensitivity")
    t.need(prep["adjusted_solution_vector"] is None and prep["native_stationarity"] is None,"native API limitation explicit")
    t.need(len(prep["cached_response_solution"])==t.N,"native response solution all lanes")
    t.need(len(record["rows"])==4,"native four rows")
    for expected,row in zip(assembly["case_vectors"],record["rows"]):
        t.need(row["case_id"]==expected["id"] and row["common_state_sha256"]==expected["state_sha256"] and row["input_complete"] is True and row["attempted"] is True and row["accepted"] is True,"native ordered case availability")
        t.need(row["status"]=="finite" and row["numerical_status"]=="ok","actual profile status")
        p=row["payload"]
        t.need(len(p["adjusted_residuals"])==t.N and finite_decimal(p["quadratic"])>=0,"native all adjusted lanes")
        for key in ("estimated_forward_sensitivity","coefficient_solve_forward_sensitivity"):
            t.need(0<=finite_decimal(p[key])<=Decimal("1e-8"),"native solve sensitivity")
        for value in p["adjusted_residuals"]:finite_decimal(value)
    return True

def reference_admission(record,request):
    reference_shape(record)
    t.need(record["interface"]=="SN-S_hat-spectral-reference/v1" and record["target_id"]==t.TARGET and record["status"]=="accepted","spectral working target status")
    expected_runtime=t.document(t.read(request["reference_runtime"])[0])
    t.need(record["runtime_before"]==record["runtime_after"]==expected_runtime,"same complete runtime before/after")
    pin_runtime(record["runtime_before"])
    t.need(record["work"]=={"eigh_attempts":1,"eigh_completed":1,"dense_check_products_attempted":2,"dense_check_products_completed":2,"profile_attempts":4,"profile_completed":4},"spectral fixed full work")
    d=record["diagnostics"]
    t.need(0<=finite_decimal(d["delta"])<1 and finite_decimal(d["minimum_eigenvalue"])>0 and finite_decimal(d["response_gram"])>0,"full spectral domain")
    t.need(0<=finite_decimal(d["response_scaled_gate"])<=Decimal("1e-8"),"spectral response gate")
    t.need(len(record["rows"])==4,"spectral complete rows")
    for case,row in zip(t.CASES,record["rows"]):
        t.need(row["case_id"]==case and row["accepted"] is True and len(row["adjusted_residuals"])==t.N,"spectral causal all-lane case")
        t.need(finite_decimal(row["quadratic"])>=0 and 0<=finite_decimal(row["scaled_gate"])<=Decimal("1e-8") and 0<=finite_decimal(row["stationarity"])<=Decimal("1e-8"),"adjusted original-matrix solve/stationarity")
    return True

def cached_response_check(native,assembly):
    # Separate native-output query, never an input to independent completion.
    raw,_=t.read(assembly["S_hat"])
    import struct
    c=array("d",(v[0] for v in struct.iter_unpack(">d",raw)))
    x=[float(v) for v in native["profile_preparation"]["cached_response_solution"]]
    norm=error=0.0
    for i in range(t.N):
        row=c[i*t.N:(i+1)*t.N]
        norm=max(norm,math.fsum(abs(v) for v in row))
        error=max(error,abs(math.fsum(a*b for a,b in zip(row,x))-1.0))
    denominator=norm*max(abs(v) for v in x)+1.0
    eta=error/denominator
    t.need(math.isfinite(eta),"native response external residual finite")
    return {"eta":repr(eta),"method":"binary64-products/math.fsum/original-S_hat","gate":None,
            "role":"external native-response diagnostic; independent spectral completion did not consume native output",
            "adjusted_solution_stationarity":None}

class Attempt:
    def __init__(self,request,out):
        self.request,self.out=request,out;self.start=time.monotonic();self.children=[];self.log_bytes=0;self.control_start=None
    def child(self,name,argv,wall,stdin=None,control=False,terminal=False):
        r=self.request["resources"]
        t.need(len(self.children)<16,"direct child count")
        t.store_bytes(self.out)
        remaining=1800-(0 if terminal else 60)-(time.monotonic()-self.start)
        if control:
            if self.control_start is None:self.control_start=time.monotonic()
            remaining=min(remaining,840-(time.monotonic()-self.control_start))
        wall=min(wall,remaining);t.need(wall>3,"remaining phase with TERM grace")
        out=self.out/(name+".stdout");err=self.out/(name+".stderr")
        env=dict(os.environ,**{k:"1" for k in THREADS},PYTHONDONTWRITEBYTECODE="1")
        def limits():
            resource.setrlimit(resource.RLIMIT_AS,(r["address_bytes"],r["address_bytes"]))
            resource.setrlimit(resource.RLIMIT_FSIZE,(r["file_bytes"],r["file_bytes"]))
            resource.setrlimit(resource.RLIMIT_CPU,(max(1,math.ceil(wall)),max(2,math.ceil(wall)+1)))
            resource.setrlimit(resource.RLIMIT_CORE,(0,0))
        entry={"name":name,"argv":argv,"start_utc":datetime.now(timezone.utc).isoformat(),"phase_wall_limit":wall,"address_limit":r["address_bytes"],"file_limit":r["file_bytes"],"thread_environment":{k:"1" for k in THREADS},"attempted":True,"returncode":None,"pid":None,"reaped":False,"timed_out":False,"leader_absent":None,"descendant_cleanup":False,"stdout":None,"stderr":None,"error":None,"actual_mapped_paths":[],"actual_thread_counts":[],"preexec_observations":0}
        self.children.append(entry)
        begun=time.monotonic();leader=None
        try:
            with out.open("xb") as fo,err.open("xb") as fe:
                with (stdin.open("rb") if stdin else open(os.devnull,"rb")) as fi:
                    leader=subprocess.Popen(argv,stdin=fi,stdout=fo,stderr=fe,env=env,start_new_session=True,preexec_fn=limits)
                    entry["pid"]=leader.pid
                    # Poll unreaped leader; bounded output/store checks during run.
                    while leader.poll() is None:
                        if name=="native" or name.startswith("control-native-"):
                            try:
                                exe=Path("/proc/"+str(leader.pid)+"/exe").resolve()
                                if exe!=Path(argv[0]).resolve():
                                    entry["preexec_observations"]+=1
                                else:
                                    pathset=set(entry["actual_mapped_paths"])
                                    for line in Path("/proc/"+str(leader.pid)+"/maps").read_text().splitlines():
                                        fields=line.split(None,5)
                                        if len(fields)==6 and fields[5].startswith("/"):
                                            t.need(not fields[5].endswith(" (deleted)"),"deleted native mapping")
                                            pathset.add(fields[5])
                                    t.need(len(pathset)<=256,"native mapped set bound")
                                    entry["actual_mapped_paths"]=sorted(pathset)
                                    for line in Path("/proc/"+str(leader.pid)+"/status").read_text().splitlines():
                                        if line.startswith("Threads:"):entry["actual_thread_counts"].append(int(line.split()[1]))
                            except FileNotFoundError:
                                pass  # exited between observation and read; no fabricated map
                        elapsed=time.monotonic()-begun
                        output_size=out.stat().st_size+err.stat().st_size
                        if elapsed>=wall-3 or self.log_bytes+output_size>r["logs_family_bytes"] or t.store_bytes(self.out)>r["store_bytes"]:
                            entry["timed_out"]=elapsed>=wall-3
                            os.killpg(leader.pid,signal.SIGTERM)
                            try:leader.wait(timeout=min(2,max(0.01,wall-elapsed)))
                            except subprocess.TimeoutExpired:
                                # Only unreaped live leader permits this signal.
                                if leader.poll() is None:os.killpg(leader.pid,signal.SIGKILL)
                                leader.wait(timeout=1)
                            t.need(False,"phase/output/store bound")
                        time.sleep(0.02)
                    entry["returncode"]=leader.wait();entry["reaped"]=True
        except BaseException as exc:
            entry["error"]=t.failure(exc)
            if leader is not None and leader.poll() is None:
                os.killpg(leader.pid,signal.SIGKILL)
                try:entry["returncode"]=leader.wait(timeout=1);entry["reaped"]=True
                except subprocess.TimeoutExpired:pass
        finally:
            entry["wall_seconds"]=time.monotonic()-begun
            entry["end_utc"]=datetime.now(timezone.utc).isoformat()
            if leader is not None:
                if entry["returncode"] is None and leader.returncode is not None:entry["returncode"]=leader.returncode;entry["reaped"]=True
                entry["leader_absent"]=not Path("/proc"+"/"+str(leader.pid)).exists()
            for path,key in ((out,"stdout"),(err,"stderr")):
                if path.exists():path.chmod(0o444);entry[key]=t.identity(path)
            self.log_bytes+=sum(entry[k]["bytes"] for k in ("stdout","stderr") if entry[k])
        t.need(entry["error"] is None and entry["reaped"] is True and entry["leader_absent"] is True,"child cleanup/execution refusal")
        t.need(self.log_bytes<=r["logs_family_bytes"],"family logs")
        return entry

def main():
    parser=argparse.ArgumentParser();parser.add_argument("--request",required=True);parser.add_argument("--request-sha",required=True);parser.add_argument("--attempt",required=True)
    args=parser.parse_args();request=pin=None
    out=Path(args.attempt);t.need(out.is_absolute() and not out.exists(),"fresh absolute attempt");out.mkdir(mode=0o700)
    attempt=Attempt(request,out)
    receipt={"interface":"SN-S_hat-bounded-controller/v1","target_id":t.TARGET,"request":pin,"status":"refused","children":attempt.children,"source_before":[],"source_after":[],"controls":[],"assembly":None,"native":None,"reference":None,"comparison":None,"error":None,"terminal_errors":[],"normalized_density":None,"joint_target":None,"source_uncertainty":None}
    fingerprint_started=False;sourcepins=[]
    try:
        request,pin=t.admitted_request(args.request,args.request_sha);receipt["request"]=pin;attempt.request=request
        t.need(request["reference_runtime"] is not None and request["reference_python"] is not None,"actual reference runtime/python admission missing")
        ports=request["source_ports"];root=Path(ports["controller.py"]["path"]).parent
        py=request["structural_python"]["path"];refpy=request["reference_python"]["path"]
        sourcepins=list(ports.values())+[request["ancestry_authority"],request["selection"],request["structural_python"],request["reference_python"],request["reference_runtime"]]
        sdk=request["sdk"]
        sourcepins += [sdk["compiler"],sdk["archive"],sdk["admission"],sdk["build_manifest"]]+sdk["headers"]+sdk["libraries"]+[request["native_loader_cache"]]
        for p in sourcepins:receipt["source_before"].append(t.read(p,False)[1])
        receipt["transport_fault_controls"]=controls.transport_fault_controls()
        receipt["projection_controls"]=controls.projection_controls()
        expected_runtime=t.document(t.read(request["reference_runtime"])[0]);pin_runtime(expected_runtime)
        pre=out/"runtime-before";fingerprint_started=True
        process=attempt.child("runtime-before",[refpy,"-B",str(root/"reference.py"),"--runtime-fingerprint","--output",str(pre)],30)
        t.need(process["returncode"]==0,"metadata preflight child")
        pre_record=reference_shape(parsed(pre/"reference.json"))
        t.need(pre_record["status"]=="runtime-only" and pre_record["runtime_before"]==pre_record["runtime_after"]==expected_runtime,"metadata exact pinned runtime")
        assembly_dir=out/"assembly"
        process=attempt.child("assembly",[py,"-B",str(root/"assemble.py"),"--request",args.request,"--request-sha",args.request_sha,"--output",str(assembly_dir)],120)
        if (assembly_dir/"assembly.json").exists():receipt["assembly"]=t.identity(assembly_dir/"assembly.json")
        t.need(process["returncode"]==0,"assembly refusal")
        assembly=t.document(t.owned_json_bytes(assembly_dir/"assembly.json",2*1024**2)[0])
        derivedpins=[assembly[k] for k in ("original_selected_B","S_hat","assembly_ledger","native_input","selection")]+[v["prediction"] for v in assembly["case_vectors"]]+[receipt["assembly"]]
        receipt["derived_before"]=[t.read(p,False)[1] for p in derivedpins]
        binary=out/"sn-working-profile"
        argv=[sdk["compiler"]["path"],"-std=c++20","-O2","-fno-fast-math","-ffp-contract=off","-frounding-math","-I",sdk["include_directory"],'-DIRRED_SN_BUILD_ID="'+sdk["build_id"]+'"',str(root/"consumer.cpp"),sdk["archive"]["path"],"-o",str(binary)]
        process=attempt.child("compile",argv,120);t.need(process["returncode"]==0,"compile refusal")
        receipt["new_executable"]=t.identity(binary);t.need(receipt["new_executable"]["bytes"]<=8*1024**2,"new native executable store reservation")
        binary.chmod(0o555)
        native_analytic=reference_analytic=None
        for name in controls.NATIVE_CONTROLS:
            if name=="original-B-refusal":
                selected=t.document(t.read(request["selection"])[0])["selected_original_indices"]
                raw=("SN_ORIGINAL_B_REFUSAL_ONLY_V1\n1657 4\n1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8 abf806d966485e64afdb359c87bffc0ecc00d05eff0a31ced66f247385df0fdc "+assembly["original_selected_B"]["sha256"]+"\n"+" ".join(str(i) for i in selected)+"\n").encode()
                header=out/"original-B-refusal.header";t.emit(header,raw,32768)
                argv=[str(binary),"--raw-B-refusal",assembly["original_selected_B"]["path"]];stdin=header
            else:argv=[str(binary),"--control",name];stdin=None
            process=attempt.child("control-native-"+name,argv,180 if name in ("analytic","conditioning") else 20,stdin,True)
            value=native_shape(parsed(Path(process["stdout"]["path"])))
            receipt["controls"].append({"kind":"native","id":name,"raw":process["stdout"],"status":value["status"]})
            if name=="analytic":t.need(process["returncode"]==0 and value["status"]=="accepted","native full-n analytic refusal");native_analytic=value
            else:
                t.need(process["returncode"]==2 and value["status"]=="refused" and not any(r["accepted"] or r["payload"] is not None for r in value["rows"]),"required native refusal availability")
                if name=="original-B-refusal":t.need(value["preparation"] is not None and value["preparation"]["status"]!="finite" and value["work"]["profile_prepare_attempts"]==value["work"]["evaluate_attempts"]==0,"original B strict preparation-only refusal")
                if name=="conditioning":
                    p=value["profile_preparation"]
                    t.need(p is not None and p["payload_available"] is False and all(p[k] is None for k in ("cached_response_backward_residual","cached_response_forward_sensitivity","gram","cached_response_solution")),"refused preparation cannot emit default zeros/empty vector")
        for name in controls.REFERENCE_CONTROLS:
            target=out/("control-reference-"+name)
            process=attempt.child("control-reference-"+name,[refpy,"-B",str(root/"reference.py"),"--request",args.request,"--request-sha",args.request_sha,"--control",name,"--output",str(target)],180,control=True)
            value=reference_shape(parsed(target/"reference.json"));receipt["controls"].append({"kind":"reference","id":name,"raw":t.identity(target/"reference.json"),"status":value["status"]})
            if name=="analytic":t.need(process["returncode"]==0,"spectral analytic refusal");reference_admission(value,request);reference_analytic=value
            else:t.need(process["returncode"]==2 and value["status"]=="refused" and value["diagnostics"] is not None and finite_decimal(value["diagnostics"]["response_scaled_gate"])>Decimal("1e-8"),"meaningful full-n spectral conditioning refusal")
        receipt["analytic_control"]=controls.analytic_admission(native_analytic,reference_analytic)
        receipt["asymmetry_algebra_control"]=controls.asymmetry_algebra_control()
        process=attempt.child("native",[str(binary)],300,Path(assembly["native_input"]["path"]))
        receipt["native"]=process["stdout"];native=parsed(Path(receipt["native"]["path"]))
        critical_maps={str(Path(p["path"]).resolve()) for p in sdk["libraries"]}|{str(binary.resolve())}
        actual_maps=set(process["actual_mapped_paths"])
        cache_maps={request["native_loader_cache"]["path"]}
        t.need(actual_maps in (critical_maps,critical_maps|cache_maps) and process["actual_thread_counts"] and set(process["actual_thread_counts"])=={1},"actual native critical maps/only declared transient loader cache/one thread")
        receipt["native_runtime_files"]=[t.identity(Path(p)) for p in sorted(actual_maps)]
        t.need(process["returncode"]==0,"native actual refusal");native_admission(native,assembly)
        target=out/"reference"
        process=attempt.child("reference",[refpy,"-B",str(root/"reference.py"),"--request",args.request,"--request-sha",args.request_sha,"--assembly",receipt["assembly"]["path"],"--output",str(target)],300)
        if (target/"reference.json").exists():receipt["reference"]=t.identity(target/"reference.json")
        t.need(process["returncode"]==0,"spectral actual refusal");reference=parsed(target/"reference.json");reference_admission(reference,request)
        t.need(reference["assembly"]["sha256"]==receipt["assembly"]["sha256"],"same assembly receipt")
        differences=[]
        with localcontext() as ctx:
            ctx.prec=96
            for n,r in zip(native["rows"],reference["rows"]):
                delta=abs(finite_decimal(n["payload"]["relative_profile_score"])-finite_decimal(r["relative_profile_score"]))
                differences.append(str(delta));t.need(delta<=Decimal("1e-6"),"same-S_hat profile difference")
        receipt["comparison"]={"case_order":list(t.CASES),"absolute_relative_profile_differences":differences,"threshold":"1e-6","engineering_only":True}
        receipt["native_response_external_check"]=cached_response_check(native,assembly)
        receipt["status"]="accepted-conditional-working-target-engineering"
    except BaseException as exc:
        receipt["error"]=t.failure(exc)
    finally:
        if fingerprint_started:
            try:
                target=out/"runtime-after"
                process=attempt.child("runtime-after",[refpy,"-B",str(root/"reference.py"),"--runtime-fingerprint","--output",str(target)],30,terminal=True)
                value=reference_shape(parsed(target/"reference.json"))
                t.need(process["returncode"]==0 and value["runtime_before"]==value["runtime_after"]==expected_runtime,"external terminal runtime")
            except BaseException as exc:receipt["status"]="refused";receipt["terminal_errors"].append(t.failure(exc))
        try:
            for p in sourcepins:receipt["source_after"].append(t.read(p,False)[1])
            if receipt["assembly"] is not None and "derivedpins" in locals():
                receipt["derived_after"]=[t.read(p,False)[1] for p in derivedpins]
                t.need(receipt["derived_before"]==receipt["derived_after"],"consumed derived payload terminal drift")
            t.need(receipt["source_before"]==receipt["source_after"],"terminal source identities")
            if pin is not None:t.need(t.read(pin,False)[1]==pin,"terminal request identity")
            t.need(time.monotonic()-attempt.start<=1800,"whole family deadline")
        except BaseException as exc:receipt["status"]="refused";receipt["terminal_errors"].append(t.failure(exc))
        receipt["wall_seconds"]=time.monotonic()-attempt.start
        receipt["retained_files"]=[t.identity(p) for p in sorted(out.rglob("*"),key=str) if p.is_file()]
        receipt["owned_bytes_before_record"]=t.store_bytes(out)
        t.emit(out/"record.json",t.encoded(receipt),4*1024**2)
        t.store_bytes(out)
        for p in out.rglob("*"):
            if p.is_file():p.chmod(0o444 if p!=out/"sn-working-profile" else 0o555)
        out.chmod(0o555)
    print(json.dumps({"status":receipt["status"],"record":t.identity(out/"record.json")},sort_keys=True))
    return 0 if receipt["status"]=="accepted-conditional-working-target-engineering" else 2

if __name__=="__main__":
    sys.exit(main())
