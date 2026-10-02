#!/usr/bin/env python3
"""Pinned-source admission, structural FITS transport and bounded native experiment.

No production numerical fit or uncertainty propagation is implemented in Python.
The rounded-table join is a source-lineage check, not an alternate physical model.
"""
import argparse
import array
import hashlib
import json
import math
from pathlib import Path
import re
import resource
import struct
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
FOLDER = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))
from packet import load, sha256, within

SDK_IDENTITY = {
    "revision": "c9b7febc01ecc85d289056faba8722a8f24a6ad8",
    "manifest_sha256": "bba6ad8d95dbbc141c6adc7c9011b1b74d6f775d69def15f71b33e6ea259136c",
    "archive_sha256": "9fbdf308e1f439cb2dc8ee297bb02997f5743ffde7e928818002488f55e20153",
    "cli_sha256": "61ed95733decd27c39349230ad3af33cfe29c082684aa1d88ed02e424649cfe9",
    "build_id": "a08f62cb32a76097e68beff6a540ccf4da903b5b898ae75921ee53818c913903",
}


FROZEN_BUDGETS = {
    "coefficient_absolute":1e-8,"coefficient_relative":1e-9,
    "quadratic_absolute":1e-7,"quadratic_relative":1e-10,
    "variance_absolute":2e-12,"variance_relative":2e-12,
}


def validate_constrained(target):
    if target["schema_version"]!=1 or target["id"]!="released-fixed44/v1" or target["active_original_indices"]!=[j for j in range(47) if j!=44] or target["fixed_coordinates"]!=[{"index":44,"value":0.0}]:
        raise ValueError("unreviewed fixed-coordinate support/order")
    if target["unchanged_reference_budgets"]!=FROZEN_BUDGETS or target["box"]["halfwidth_multiplier"]!=10 or target["box"]["endpoints"]!="closed":
        raise ValueError("unreviewed fixed target box/budgets")
    if [(v["row"],v["delta_y"]) for v in target["calibration_sensitivities"]]!=[(3211,.10),(3213,.032),(3214,.0263)]:
        raise ValueError("unreviewed source calibration shifts")


def box_support(beta,prior):
    """Structural source-coordinate box admission; no solver/density in Python."""
    if len(beta)!=47 or len(prior)!=47 or prior[44]!=[0.,0.] or beta[44]!=0 or not all(math.isfinite(v) for v in beta):
        raise ValueError("fixed literal-zero finite source support required")
    bounds=[];outside=[]
    for j,(center,width) in enumerate(prior):
        if not math.isfinite(center) or not math.isfinite(width) or (j!=44 and width<=0):
            raise ValueError("positive finite active box width required")
        low,high=center-10*width,center+10*width
        if not math.isfinite(low) or not math.isfinite(high) or (j!=44 and low>=high):
            raise ValueError("finite nondegenerate active box required")
        bounds.append({"original_index":j,"lower":low,"upper":high})
        if not low<=beta[j]<=high:outside.append(j)
    return {"bounds":bounds,"profile_inside_box":not outside,"outside_original_indices":outside,"meaning":"support only; no posterior normalization/inference"}


def constrained_source_audit(store,sources):
    x=array.array("d");x.frombytes((store/"X.f64").read_bytes())
    c=array.array("d");c.frombytes((store/"C.f64").read_bytes())
    expected={37:list(range(2150,2593))+[3213],38:list(range(3130))+[3207,3208],39:list(range(2648,3130))+[3214],40:list(range(2593,2648)),42:list(range(3130,3207))+list(range(3215,3492)),44:[3210],45:list(range(2648,3061))+[3211]}
    for j,indices in expected.items():
        if [i for i in range(3492) if x[i*47+j]!=0]!=indices or any(x[i*47+j]!=1 for i in indices):
            raise ValueError("source coordinate support changed")
    if any(c[3210*3492+i]!=0 or c[i*3492+3210]!=0 for i in range(3492) if i!=3210):
        raise ValueError("fixed-coordinate row not covariance-isolated")
    table=photometry_rows(sources["table2.tex"].decode());groups=[]
    for host,instrument,start,stop in (("N4258","HST",2150,2593),("M31","HST",2593,2648),("LMC","GRND",2648,2918),("SMC","GRND",2918,3061),("LMC","HST",3061,3130)):
        failures=[];locators=set()
        for i in range(start,stop):
            lp,metal=x[i*47+41],x[i*47+43];hp,hm=binary32_half_ulp(lp),binary32_half_ulp(metal)
            matches=[r for r in table if r["host"]==host and r["instrument"]==instrument and lp+hp>=r["period_lower"]-1e-14 and lp-hp<=r["period_upper"]+1e-14 and abs(metal-r["metal"])<=r["metal_half"]+hm]
            if not matches:failures.append({"row":i,"log_period_minus1":lp,"metallicity":metal})
            locators.update(r["line"] for r in matches)
        groups.append({"host":host,"instrument":instrument,"rows":[start,stop-1],"row_count":stop-start,"unmatched":failures,"table2_lines":sorted(locators)})
    return {"original_supports":expected,"fixed44_row_covariance_isolated":True,"anchor_ancillary_join":groups,"qualification":"partial lineage; failed joins retained without adapted tolerance"}


def validate_lineage(m):
    if m["schema_version"] != 1 or m["source_revision"] != "c447f0fea703fcd0fff57de5000947b5ca81286b":
        raise ValueError("unreviewed lineage version/release")
    if [a["index"] for a in m["axes"]] != list(range(47)) or [a["id"] for a in m["axes"]] != [f"released-parameter-{j}" for j in range(47)]:
        raise ValueError("exact ordered full47 axes required")
    if m["decode"]["rows"] != 3492 or m["decode"]["columns"] != 47 or m["decode"]["mask"] != "none":
        raise ValueError("unchanged full design required")
    if len({s["name"] for s in m["sources"]}) != len(m["sources"]):
        raise ValueError("duplicate source identity")
    if m["target"]["reference_budgets"]!=FROZEN_BUDGETS:
        raise ValueError("unreviewed numerical comparison budgets")
    rounded=m["target"]["paper_rounded_coordinate"]
    if rounded["value"]!=9.318 or rounded["half_last_printed_digit"]!=0.0005:
        raise ValueError("unreviewed source rounded-coordinate target")
    if m["constraint_rows"]["indices"]!=list(range(3207,3215)) or m["constraint_rows"]["sensitivity_control"]["delta_y"]!=0.01:
        raise ValueError("unreviewed sensitivity rows/perturbation")
    if m["axes"][44]["physical_identity"] is not None or m["axes"][44]["status"] != "unresolved":
        raise ValueError("column44 has no admitted physical identity")
    if m["target"]["h0_coordinate"]["index"] != 46:
        raise ValueError("source contrast order changed")
    if m["target"]["prior_in_native_target"] != "none; all 47 columns retained":
        raise ValueError("native target measure changed")
    if m["released_mcmc_prior"]["fixed_coordinates"] != [{"index":44,"value":0,"halfwidth":0}]:
        raise ValueError("released zero-width prior distinction changed")


def admitted_sources(directory, m):
    result = {}
    for item in m["sources"]:
        p = within(directory, item["name"])
        if p.stat().st_size != item["bytes"]:
            raise ValueError("source byte count differs: " + item["name"])
        with p.open("rb") as stream:
            b = stream.read(item["bytes"]+1)
            if stream.read(1):raise ValueError("source grew during bounded read")
        if len(b) != item["bytes"] or hashlib.sha256(b).hexdigest() != item["sha256"]:
            raise ValueError("source bytes/hash differ: " + item["name"])
        result[item["name"]] = b
    return result


def fits_image(b, expected_shape):
    """Only the exact pinned unscaled binary32 primary image profile is admitted."""
    header = {}
    end = None
    for offset in range(0, min(len(b), 65536), 80):
        card = b[offset:offset+80].decode("ascii")
        key = card[:8].strip()
        if key == "END":
            end = ((offset + 80 + 2879)//2880)*2880
            break
        if card[8:10] == "= ":
            if key in header:
                raise ValueError("duplicate FITS card")
            header[key] = card[10:].split("/",1)[0].strip()
    if end is None or header.get("SIMPLE") != "T" or header.get("BITPIX") != "-32":
        raise ValueError("unsupported FITS primary profile")
    ndim = int(header["NAXIS"])
    shape = tuple(int(header[f"NAXIS{i}"]) for i in range(1,ndim+1))[::-1]
    if shape != tuple(expected_shape) or "BSCALE" in header or "BZERO" in header:
        raise ValueError("FITS dimension/scaling differs")
    count = math.prod(shape)
    data_end = end + 4*count
    if len(b) != ((data_end+2879)//2880)*2880 or any(b[data_end:]):
        raise ValueError("FITS byte/padding count differs")
    values = array.array("f")
    values.frombytes(b[end:data_end])
    if sys.byteorder == "little":
        values.byteswap()
    if not all(math.isfinite(v) for v in values):
        raise ValueError("nonfinite source coordinate")
    return values


def binary32_half_ulp(value):
    bits = struct.unpack(">I", struct.pack(">f", abs(value)))[0]
    upper = struct.unpack(">f", struct.pack(">I", bits+1))[0]
    return (upper-abs(value))/2


def photometry_rows(text):
    rows = []
    for line_number, line in enumerate(text.splitlines(),1):
        fields = [v.strip() for v in line.split("&")]
        if len(fields)!=11 or not re.fullmatch(r"\d+(?:\.\d+)?",fields[4]):
            continue
        period = float(fields[4])
        decimals = len(fields[4].split(".")[1]) if "." in fields[4] else 0
        half = 0.5*10**(-decimals)
        metal = float(fields[9])
        metal_decimals = len(fields[9].split(".")[1]) if "." in fields[9] else 0
        rows.append({"host":fields[0],"id":fields[3],"line":line_number,"instrument":fields[10].replace("\\","").strip(),
                     "period_lower":math.log10(period-half)-1,
                     "period_upper":math.log10(period+half)-1,
                     "metal":metal,"metal_half":0.5*10**(-metal_decimals)})
    if not rows:
        raise ValueError("no primary table rows")
    return rows


def host_join(x, table, expected_hosts, n=3492, p=47, stop=2150):
    """Intersect all rounded-source host candidates in each original host group.

    Bounds are source printed decimal half-units plus binary32 half-ulp; 1e-14
    in log-period covers library log10 evaluation far below either precision.
    This is a bounded reading/join convention, not a calibrated error budget.
    """
    if len(x)!=n*p or len(expected_hosts)!=37:
        raise ValueError("host join dimensions")
    report = []
    for j in range(37):
        indices = [i for i in range(stop) if x[i*p+j]!=0]
        if not indices or any(x[i*p+j]!=1 for i in indices):
            raise ValueError("host-column support changed")
        candidates = None
        matched_lines = set()
        for i in indices:
            lp, metal = x[i*p+41], x[i*p+43]
            half_lp, half_metal = binary32_half_ulp(lp), binary32_half_ulp(metal)
            matches = [r for r in table if lp+half_lp >= r["period_lower"]-1e-14
                       and lp-half_lp <= r["period_upper"]+1e-14
                       and abs(metal-r["metal"]) <= r["metal_half"]+half_metal]
            hosts = {r["host"] for r in matches}
            candidates = hosts if candidates is None else candidates & hosts
            matched_lines.update(r["line"] for r in matches if r["host"]==expected_hosts[j])
        if candidates != {expected_hosts[j]}:
            raise ValueError(f"host {j} unresolved/changed: {sorted(candidates or [])}")
        report.append({"column":j,"host":expected_hosts[j],"row_count":len(indices),
                       "original_row_indices":indices,"table2_source_lines":sorted(matched_lines)})
    if sum(r["row_count"] for r in report)!=stop:
        raise ValueError("host rows not partitioned exactly once")
    return report


def source_audit(sources,m,store):
    names = {s["name"] for s in m["sources"]}
    def image(prefix,shape):
        name=next(n for n in names if n.startswith(prefix) and n.endswith(".fits"))
        return fits_image(sources[name],shape)
    c=image("allc",(3492,3492)); l=image("alll",(47,3492)); y=image("ally",(3492,))
    x=array.array("d",(l[j*3492+i] for i in range(3492) for j in range(47)))
    for name,values in (("C.f64",c),("X.f64",x),("y.f64",y)):
        a=array.array("d",values)
        if sys.byteorder!="little":a.byteswap()
        b=a.tobytes()
        if hashlib.sha256(b).hexdigest()!=m["decode"]["canonical_sha256"][name]:
            raise ValueError("canonical decode identity differs")
        (store/name).write_bytes(b)
    # Exact source boundary, not diagonal-error repair or a factorization.
    if any(c[i*3492+j]!=c[j*3492+i] for i in range(3492) for j in range(i)):
        raise ValueError("source covariance asymmetric")
    if [i for i in range(3492) if x[i*47+46]!=0] != list(range(3215,3492)) or any(x[i*47+46]!=-1 or x[i*47+42]!=1 for i in range(3215,3492)):
        raise ValueError("released final contrast support differs")
    hosts=[a["physical_identity"].split(":",1)[1] for a in m["axes"][:37]]
    joined=host_join(x,photometry_rows(sources["table2.tex"].decode()),hosts)
    prior=[list(map(float,line.split())) for line in sources["lstsq_results.txt"].decode().splitlines() if line.strip()]
    if len(prior)!=47 or prior[44]!=[0.,0.] or any(v[1]<=0 for j,v in enumerate(prior) if j!=44):
        raise ValueError("released initialization/prior table differs")
    constraints=[]
    for i in range(3207,3215):
        constraints.append({"row":i,"y":y[i],"Cii":c[i*3492+i],"design_nonzero":[{"column":j,"value":x[i*47+j]} for j in range(47) if x[i*47+j]!=0]})
    return {"host_join":joined,"constraints":constraints,"prior_halfwidths":[10*v[1] for v in prior],"fixed_coordinate44":True}


def fingerprint(source,sdk):
    def git(*args):return subprocess.check_output(["git","-C",str(source),*args],text=True,timeout=30).strip()
    # Verify immutable pinned Git blobs; the primary checkout may carry newer work.
    if git("rev-parse",SDK_IDENTITY["revision"])!=SDK_IDENTITY["revision"]:
        raise ValueError("pinned engine commit unavailable")
    manifest=sdk/"build-manifest.json"
    if sha256(manifest)!=SDK_IDENTITY["manifest_sha256"]:raise ValueError("SDK manifest differs")
    build=load(manifest)
    if build["git_head"]!=SDK_IDENTITY["revision"] or build["git_status"] or build["build_id"]!=SDK_IDENTITY["build_id"]:
        raise ValueError("SDK build source differs")
    content={k:v for k,v in build.items() if k not in ("build_id","git_head","git_status")}
    if hashlib.sha256(json.dumps(content,sort_keys=True,separators=(",",":")).encode()).hexdigest()!=build["build_id"]:
        raise ValueError("SDK build identity not reconstructible")
    for name,digest in build["sources"].items():
        blob=subprocess.check_output(["git","-C",str(source),"show",SDK_IDENTITY["revision"]+":"+name],timeout=30)
        if hashlib.sha256(blob).hexdigest()!=digest:raise ValueError("immutable engine source differs: "+name)
    expected={Path(n).name for n in build["sources"] if n.startswith("cpp/include/irred/")}
    actual={p.name for p in (sdk/"include/irred").iterdir()}
    if actual!=expected:raise ValueError("complete SDK header inventory differs")
    for name in expected:
        if sha256(sdk/"include/irred"/name)!=build["sources"]["cpp/include/irred/"+name]:raise ValueError("SDK header differs")
    for path,key in ((sdk/"lib/libirred_core.a","archive_sha256"),(sdk/"bin/irred","cli_sha256")):
        if sha256(path)!=SDK_IDENTITY[key]:raise ValueError("SDK artifact differs")
    for path,digest in ((Path("/usr/bin/c++"),build["compiler_executable_digest"]),(Path(build["standard_library"]),build["standard_library_digest"])):
        if sha256(path)!=digest:raise ValueError("compiler/standard library differs")
    return {str(p):sha256(p) for p in sorted(sdk.rglob("*")) if p.is_file()}


def child(command,store,label,timeout,output_limit=1048576):
    def limits():
        resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
        resource.setrlimit(resource.RLIMIT_FSIZE,(output_limit,output_limit))
    with (store/(label+".out")).open("xb") as out,(store/(label+".err")).open("xb") as err:
        subprocess.run(command,stdout=out,stderr=err,cwd=ROOT,check=True,timeout=timeout,preexec_fn=limits)
    return (store/(label+".out")).read_text()



def validate_native(output,target=None):
    if output["method"]!="retained-whitened-pivoted-householder-qr/v1" or output["variance_method"]!="retained-qr-linear-estimator-variance/v1":
        raise ValueError("unreviewed native method identity")
    if len(output["coefficients"])!=47 or not all(math.isfinite(v) for v in output["coefficients"]):
        raise ValueError("full finite coefficient output required")
    if not all(math.isfinite(output[k]) for k in ("quadratic","relative_log_score","variance46","variance_sensitivity","stationarity")):
        raise ValueError("nonfinite native scalar")
    if output["quadratic"]<0 or output["variance46"]<=0 or output["relative_log_score"]!=-output["quadratic"]/2:
        raise ValueError("native relative-score/variance contract differs")
    if output["variance_sensitivity"]<0 or output["variance_sensitivity"]>1e-10 or output["stationarity"]<0 or output["stationarity"]>1e-10:
        raise ValueError("native numerical screen differs")
    expected=[(v["row"],v["delta_y"]) for v in target["calibration_sensitivities"]] if target else [(i,.01) for i in range(3207,3215)]
    if [(v["row"],v["delta_y"]) for v in output["sensitivities"]]!=expected:
        raise ValueError("native sensitivity order/perturbation differs")
    if not all(math.isfinite(v[k]) for v in output["sensitivities"] for k in ("beta46_plus","beta46_minus","q_plus","q_minus")):
        raise ValueError("nonfinite native sensitivity")
    if target and (output.get("rank")!=46 or output["coefficients"][44]!=0):
        raise ValueError("literal-zero fixed46 native output required")


def execute(args):
    store=within(ROOT,"results/released-ladder/"+args.name)
    store.mkdir(parents=True,exist_ok=False)
    record={"status":"started","engine_identity":SDK_IDENTITY,"gates":{"execution":"unassessed","numerical":"unassessed","inference":"unassessed","interpretation":"conditional; full source qualification blocked"}}
    started=time.monotonic(); before=None
    try:
        m=load(FOLDER/"lineage.json"); validate_lineage(m)
        target=load(FOLDER/"constrained.json") if getattr(args,"target","full47")=="released-fixed44/v1" else None
        if target:validate_constrained(target)
        sources=admitted_sources(args.sources,m)
        record["source_hashes"]={n:hashlib.sha256(b).hexdigest() for n,b in sources.items()}
        record["packet_hashes"]={p.name:sha256(p) for p in FOLDER.iterdir() if p.is_file()}
        record["source_audit"]=source_audit(sources,m,store)
        if target:
            record["target"]=target
            record["constrained_source_audit"]=constrained_source_audit(store,sources)
        before=fingerprint(args.engine_source,args.sdk);record["sdk_before"]=before
        record["engine_tree"]=subprocess.check_output(["git","-C",str(args.engine_source),"rev-parse",SDK_IDENTITY["revision"]+"^{tree}"],text=True).strip()
        # Capture exact adapter sources even during a diagnostic uncommitted run.
        for name in ("controller.py","consumer.cpp","lineage.json","reference.py","constrained.json"):
            (store/name).write_bytes((FOLDER/name).read_bytes())
        record["repository_head"]=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip()
        record["repository_status"]=subprocess.check_output(["git","status","--porcelain"],cwd=ROOT,text=True)
        executable=store/"consumer"
        command=["/usr/bin/c++","-std=c++20","-O2","-Wall","-Wextra","-Wpedantic","-fno-fast-math","-ffp-contract=off",str(store/"consumer.cpp"),"-I",str(args.sdk/"include"),str(args.sdk/"lib/libirred_core.a"),"-o",str(executable)]
        record["compiler_command"]=command
        child(command,store,"compile",120,8*1024**2);record["consumer_sha256"]=sha256(executable)
        record["native_controls"]=json.loads(child([str(executable),"--self-test"],store,"native-controls",30))
        output=json.loads(child([str(executable),str(store)]+([target["id"]] if target else []),store,"native",900))
        record["native"]=output; record["gates"]["execution"]="passed"
        validate_native(output,target)
        if target:
            prior=[list(map(float,line.split())) for line in sources["lstsq_results.txt"].decode().splitlines() if line.strip()]
            record["box_support"]=box_support(output["coefficients"],prior)
        coordinate=m["target"]["paper_rounded_coordinate"]
        if abs(output["coefficients"][46]-coordinate["value"])>coordinate["half_last_printed_digit"]:raise ValueError("paper rounded coordinate mismatch")
        if args.reference_python:
            record["reference_command"]=[str(args.reference_python),str(store/"reference.py"),str(store)]+([target["id"]] if target else [])
            reference=json.loads(child(record["reference_command"],store,"reference",900))
            record["reference"]=reference;record["comparisons"]=compare(output,reference,m,target)
            record["gates"]["numerical"]="passed named fixed46 QR/SVD and source constraint-mean sensitivity comparisons" if target else "passed named full47 QR/SVD and synthetic row-sensitivity comparisons"
        else:
            record["gates"]["numerical"]="native accepted; external comparison not run"
        if fingerprint(args.engine_source,args.sdk)!=before:raise ValueError("SDK changed during execution")
        if admitted_sources(args.sources,m)!=sources:raise ValueError("source changed during execution")
        if {p.name:sha256(p) for p in FOLDER.iterdir() if p.is_file()}!=record["packet_hashes"]:raise ValueError("packet changed during execution")
        record["status"]="completed"
    except Exception as e:
        record["status"]="failed";record["error"]=str(e)
    finally:
        # Retain post-run integrity evidence even when a child/reference fails.
        integrity_errors=[]
        if before is not None:
            try:
                record["sdk_after"]=fingerprint(args.engine_source,args.sdk)
                if record["sdk_after"]!=before:integrity_errors.append("SDK changed")
            except Exception as e:integrity_errors.append("SDK final verification: "+str(e))
        if "source_hashes" in record:
            try:
                final_sources=admitted_sources(args.sources,m)
                record["source_hashes_after"]={n:hashlib.sha256(b).hexdigest() for n,b in final_sources.items()}
                if record["source_hashes_after"]!=record["source_hashes"]:integrity_errors.append("source changed")
            except Exception as e:integrity_errors.append("source final verification: "+str(e))
        if "packet_hashes" in record:
            record["packet_hashes_after"]={p.name:sha256(p) for p in FOLDER.iterdir() if p.is_file()}
            if record["packet_hashes_after"]!=record["packet_hashes"]:integrity_errors.append("packet changed")
        record["integrity_errors"]=integrity_errors
        if integrity_errors:
            record["status"]="failed"
            record.setdefault("error","post-run identity verification failed")
        record["elapsed_seconds"]=time.monotonic()-started
        usage=resource.getrusage(resource.RUSAGE_CHILDREN)
        record["children_resources"]={"user_cpu_seconds":usage.ru_utime,"system_cpu_seconds":usage.ru_stime,"maximum_rss_kib":usage.ru_maxrss,"scope":"all compiler/native/reference children of this controller"}
        (store/"record.json").write_text(json.dumps(record,indent=2)+"\n")
        for p in store.iterdir():
            if p.is_file():p.chmod(0o444)
    print(json.dumps({"status":record["status"],"record":str(store/"record.json"),"gates":record["gates"]}))
    return 0 if record["status"]=="completed" else 1


def compare(native,reference,m,target=None):
    b=m["target"]["reference_budgets"];checks=[]
    if reference.get("rank")!=(46 if target else 47) or set(reference["algorithms"])!={"LAPACK_gesdd_SVD","LAPACK_pivoted_QR"}:
        raise ValueError("complete independent reference algorithms/rank required")
    def scalar(name,a,r,kind):
        if not math.isfinite(a) or not math.isfinite(r):raise ValueError("nonfinite reference/output")
        budget=b[kind+"_absolute"]+b[kind+"_relative"]*abs(r)
        delta=abs(a-r)
        checks.append({"name":name,"absolute_difference":delta,"budget":budget,"fraction":delta/budget})
        if delta>budget:raise ValueError("frozen comparison failed: "+name)
    for algorithm,value in reference["algorithms"].items():
        if len(value["coefficients"])!=47:raise ValueError("reference full47 required")
        if target and (native["coefficients"][44]!=0 or value["coefficients"][44]!=0):
            raise ValueError("reference fixed coordinate must be literal zero")
        for j,(a,r) in enumerate(zip(native["coefficients"],value["coefficients"])):scalar(f"{algorithm}/beta{j}",a,r,"coefficient")
        scalar(algorithm+"/quadratic",native["quadratic"],value["quadratic"],"quadratic")
        scalar(algorithm+"/variance46",native["variance46"],value["variance46"],"variance")
    rows=[v["row"] for v in target["calibration_sensitivities"]] if target else list(range(3207,3215))
    if [v["row"] for v in native["sensitivities"]]!=rows or [v["row"] for v in reference["sensitivities"]]!=rows:
        raise ValueError("sensitivity row order changed")
    for a,r in zip(native["sensitivities"],reference["sensitivities"]):
        delta=next(v["delta_y"] for v in target["calibration_sensitivities"] if v["row"]==a["row"]) if target else .01
        if a["delta_y"]!=delta or (target and r.get("delta_y")!=delta):raise ValueError("sensitivity perturbation changed")
        for field in ("beta46_plus","beta46_minus"):scalar(f"row{a['row']}/{field}",a[field],r[field],"coefficient")
        for field in ("q_plus","q_minus"):scalar(f"row{a['row']}/{field}",a[field],r[field],"quadratic")
    return checks


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sources",type=Path,required=True)
    p.add_argument("--engine-source",type=Path,required=True)
    p.add_argument("--sdk",type=Path,required=True)
    p.add_argument("--name",required=True)
    p.add_argument("--reference-python",type=Path)
    p.add_argument("--target",choices=("full47","released-fixed44/v1"),default="full47")
    args=p.parse_args()
    args.sources=args.sources.resolve();args.engine_source=args.engine_source.resolve();args.sdk=args.sdk.resolve()
    if args.reference_python:args.reference_python=args.reference_python.absolute()
    sys.exit(execute(args))
