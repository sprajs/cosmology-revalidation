"""Independent same-S_hat full-mode spectral profile; no production fallback."""
import argparse
import contextlib
import ctypes
from decimal import Decimal
from fractions import Fraction
import hashlib
import io
import json
import math
import os
from pathlib import Path
import sys
import time
import transport as t
import controls

THREADS = ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","BLIS_NUM_THREADS","NUMEXPR_NUM_THREADS","VECLIB_MAXIMUM_THREADS")

def runtime_fingerprint(np):
    import importlib.metadata
    import numpy.linalg
    import threadpoolctl
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        np.show_config()
    t.need(threadpoolctl.__version__=="3.6.0","pinned private threadpoolctl3.6.0")
    process_symbols=ctypes.CDLL(None)
    threadpoolctl.ThreadpoolController._system_libraries["libc"]=process_symbols
    pools = threadpoolctl.threadpool_info()
    rounding=process_symbols.fegetround
    rounding.argtypes=[];rounding.restype=ctypes.c_int
    rounding_code=rounding()
    t.need(rounding_code==0 and np.finfo(np.longdouble).nmant==63,"pinned Linux RN and 64-significand-bit wide reductions")
    t.need(all(type(p.get("num_threads")) is int and p["num_threads"] == 1 for p in pools) and pools, "actual one-thread BLAS getters")
    t.need(all(os.environ.get(k) == "1" for k in THREADS), "thread environment")
    launch = {k:os.environ.get(k) for k in ("PYTHONPATH","PYTHONNOUSERSITE")}
    t.need(launch=={"PYTHONPATH":"/home/szymon/Projects/reproducible/.venv/lib/python3.12/site-packages:/home/szymon/.codex/worktrees/a4f7/reproducible/.work/sn-threadpoolctl-3.6.0-20261003","PYTHONNOUSERSITE":"1"},"exact inherited reference launch context")
    t.need(len(sys.path)<=128 and all(type(p) is str and len(p.encode())<=4096 for p in sys.path),"bounded actual sys.path")
    # Metadata/helper imports finish before recording actual sets. No eigensolve,
    # GL, profile, matrix, or scientific prewarm is performed here.
    distribution = importlib.metadata.distribution("numpy")
    files = set()
    def add_file(path):
        path=Path(path).resolve()
        if path not in files:
            t.need(len(files)<8192,"runtime file cardinality before admission")
            files.add(path)
    package = Path(np.__file__).resolve().parent
    for p in package.rglob("*"):
        if p.is_file() and (p.suffix in (".py",".pyc",".so") or ".so." in p.name):
            add_file(p)
    for sibling in package.parent.glob("numpy.libs/*"):
        if sibling.is_file():
            add_file(sibling)
    for p in distribution.files or ():
        path = Path(distribution.locate_file(p)).resolve()
        if path.is_file() and (".dist-info" in str(path) or path.suffix == ".pyc"):
            add_file(path)
    modules = []
    for name,module in sys.modules.items():
        if module is None:
            continue
        t.need(len(modules)<4096,"runtime module cardinality before admission")
        origin = getattr(module,"__file__",None)
        locations = getattr(getattr(module,"__spec__",None),"submodule_search_locations",None)
        entry = {"name": name, "origin": None, "cached": None, "namespace_paths": None}
        if origin and Path(origin).is_file():
            path = Path(origin).resolve();add_file(path);entry["origin"] = str(path)
            cached = getattr(module,"__cached__",None)
            if cached and Path(cached).is_file():
                cached = Path(cached).resolve();add_file(cached);entry["cached"] = str(cached)
        elif locations is not None:
            entry["namespace_paths"] = [str(Path(p).resolve()) for p in locations]
        modules.append(entry)
    modules.sort(key=lambda p:p["name"])
    maps = set()
    with Path("/proc/self/maps").open() as map_stream:
        for line in map_stream:
            fields = line.split(None,5)
            if len(fields)==6 and fields[5].startswith("/"):
                lexical = fields[5].rstrip("\n")
                t.need(not lexical.endswith(" (deleted)"), "deleted mapped file")
                resolved = Path(lexical).resolve()
                t.need(resolved.is_file(), "mapped regular file")
                add_file(resolved)
                pair=(lexical,str(resolved))
                if pair not in maps:
                    t.need(len(maps)<256,"runtime map cardinality before admission")
                    maps.add(pair)
    maps = [dict(lexical_path=a,resolved_path=b) for a,b in sorted(maps)]
    executable = Path(sys.executable).resolve();add_file(executable)
    inventory = [];total=0
    for path in sorted(files,key=str):
        # One bounded consumed stream: descriptor stat enforces leaf and remaining
        # aggregate allowance before reading; failures retain actual read evidence.
        _,sealed = t.read({"path":str(path),"bytes":None,"sha256":None},False,
            limit=min(256*1024**2,1024**3-total),allow_unpinned=True,unknown_size=True)
        total+=sealed["bytes"];inventory.append(sealed)
    record = {"schema":"SN-spectral-runtime/v1", "python_version":sys.version, "python_executable":str(executable),
              "sys_prefix":str(Path(sys.prefix).resolve()), "sys_base_prefix":str(Path(sys.base_prefix).resolve()),
              "numpy_version":np.__version__, "backend":"numpy.linalg.eigh-UPLO-L-LAPACK-syevd",
              "thread_environment":{k:os.environ.get(k) for k in THREADS},"launch_environment":launch,"sys_path":list(sys.path), "numeric_configuration":{"show_config":buf.getvalue(),"threadpools":pools,
              "threadpoolctl_libc_binding":"3.6.0/ThreadpoolController._system_libraries[libc]=ctypes.CDLL(None)/already-loaded-process-symbols/no-discovery-children",
              "fegetround_linux_code":rounding_code,"rounding":"FE_TONEAREST/Linux-code-0",
              "longdouble":{k:str(getattr(np.finfo(np.longdouble),k)) for k in ("bits","nmant","iexp","eps","tiny","max")}},
              "imported_modules":modules,"mapped_libraries":maps,"inventory":inventory}
    record["inventory_sha256"] = hashlib.sha256(t.encoded(inventory)).hexdigest()
    return t.pack_runtime(record)

def scalar(np,x):
    t.need(np.isfinite(x), "nonfinite diagnostic/scalar")
    return np.format_float_scientific(np.longdouble(x),unique=False,precision=36,trim="k")

def norm_inf(np,a):
    if a.ndim==2:
        return np.max(np.sum(np.abs(a),axis=1,dtype=np.longdouble))
    return np.max(np.abs(a))

def reportcast(np,path,x):
    cast=x.astype(np.float64)
    t.need(np.isfinite(cast).all() and not np.any((x!=0)&(cast==0)),"unrepresentable reportcast solution")
    return t.f64_payload(path,cast.tolist())

def normalized_residual(np,c,x,b,cnorm):
    error = c @ x - b
    numerator = norm_inf(np,error)
    denominator = cnorm*norm_inf(np,x)+norm_inf(np,b)
    if denominator==0:
        t.need(np.array_equal(x,np.zeros_like(x)) and np.array_equal(b,np.zeros_like(b)) and numerator==0, "zero residual denominator")
        return np.longdouble(0)
    return numerator/denominator

def evaluate(np,c,residuals,out,receipt):
    t.need(np.__version__=="2.2.6", "exact numpy2.2.6")
    t.need(c.shape==(t.N,t.N) and np.isfinite(c).all() and np.array_equal(c,c.T), "full original derived matrix admission")
    t.need(len(residuals)==4, "four ordered residuals")
    receipt["work"]["eigh_attempts"] += 1
    lam,q = np.linalg.eigh(c,UPLO="L")
    receipt["work"]["eigh_completed"] += 1
    # Earned eigenpairs are retained even if positivity or later diagnostics fail.
    receipt["eigenvalues"] = t.f64_payload(out/"eigenvalues.f64be",lam.tolist())
    receipt["eigenvectors"] = t.f64_payload(out/"eigenvectors-rowmajor.f64be",q.ravel().tolist())
    t.need(lam.shape==(t.N,) and q.shape==(t.N,t.N) and np.isfinite(lam).all() and np.isfinite(q).all() and (lam>0).all(), "all eigenpairs finite/positive; no mode clipping")
    w = np.longdouble
    cw,qw,lw = c.astype(w),q.astype(w),lam.astype(w)
    one = np.ones(t.N,dtype=w)
    cnorm = norm_inf(np,cw)
    t.need(np.isfinite(cnorm) and cnorm>0, "positive matrix norm")
    receipt["work"]["dense_check_products_attempted"] += 1
    reconstruct = (q*lam)@q.T
    receipt["work"]["dense_check_products_completed"] += 1
    rho = norm_inf(np,reconstruct.astype(w)-cw)/cnorm
    del reconstruct
    receipt["work"]["dense_check_products_attempted"] += 1
    orthogonal = q.T@q
    receipt["work"]["dense_check_products_completed"] += 1
    delta = norm_inf(np,orthogonal.astype(w)-np.eye(t.N,dtype=w))
    del orthogonal
    kappa = lw[-1]/lw[0]
    b = qw.T@one
    g = np.sum(b*b/lw,dtype=w)
    t.need(np.isfinite(g) and g>0, "finite positive response gram")
    response = qw@(b/lw)
    receipt["response_solution"] = reportcast(np,out/"response-solution-reportcast.f64be",response)
    eta_response = normalized_residual(np,cw,response,one,cnorm)
    receipt["diagnostics"] = {"rho":scalar(np,rho),"delta":scalar(np,delta),"kappa_lambda":scalar(np,kappa),
        "matrix_norm_inf":scalar(np,cnorm),"response_gram":scalar(np,g),"eta_response":scalar(np,eta_response),
        "minimum_eigenvalue":scalar(np,lw[0]),"maximum_eigenvalue":scalar(np,lw[-1]),"stationarity_threshold":"1e-8",
        "dimension_roundoff_floor":scalar(np,w(t.N)*w(2)**-52),"engineering_only":True}
    base = kappa*(rho+delta+eta_response+w(t.N)*w(2)**-52)
    receipt["diagnostics"]["response_scaled_gate"] = scalar(np,base)
    t.need(delta<1 and np.isfinite(base) and base<=w("1e-8"), "spectral response engineering gate")
    for case_id,raw in zip(t.CASES,residuals):
        receipt["work"]["profile_attempts"] += 1
        r = np.asarray(raw,dtype=np.float64)
        t.need(r.shape==(t.N,) and np.isfinite(r).all(), "exact ordered residual vector")
        rw = r.astype(w)
        a = qw.T@rw
        m = np.sum(a*b/lw,dtype=w)/g
        modes = a-m*b
        quadratic = np.sum(modes*modes/lw,dtype=w)
        adjusted = rw-m*one
        x = qw@(modes/lw)
        eta = normalized_residual(np,cw,x,adjusted,cnorm)
        stationarity_signed = np.sum(x,dtype=w)
        denom = np.sqrt(w(t.N))*np.sqrt(np.sum(x*x,dtype=w))
        if denom==0:
            t.need(np.array_equal(x,np.zeros_like(x)) and stationarity_signed==0, "zero stationarity vector")
            stationarity = w(0)
        else:
            stationarity = abs(stationarity_signed)/denom
        scaled = kappa*(rho+delta+max(eta,eta_response)+w(t.N)*w(2)**-52)
        row = {"case_id":case_id,"M_coefficient":scalar(np,m),"quadratic":scalar(np,quadratic),
               "relative_profile_score":scalar(np,-quadratic/2),"eta_adjusted":scalar(np,eta),"scaled_gate":scalar(np,scaled),
               "stationarity_signed":scalar(np,stationarity_signed),"stationarity":scalar(np,stationarity),
               "adjusted_residuals":[scalar(np,v) for v in adjusted],"accepted":False}
        row["adjusted_solution_reportcast"] = None
        receipt["rows"].append(row)
        row["adjusted_solution_reportcast"] = reportcast(np,out/(case_id+".adjusted-solution-reportcast.f64be"),x)
        t.need(quadratic>=0 and scaled<=w("1e-8") and stationarity<=w("1e-8"), "spectral adjusted solve/stationarity screen")
        row["accepted"] = True
        receipt["work"]["profile_completed"] += 1
    receipt["status"] = "accepted"

def sensitivity(np,ledger,receipt):
    minimum = np.longdouble(receipt["diagnostics"]["minimum_eigenvalue"])
    def value(key):
        p = ledger[key]
        # Longdouble estimate of the exact rational norm; not outward/validated.
        return np.longdouble(p["numerator"])/np.longdouble(p["denominator"])
    rnorm = value("R_norm_inf_exact")
    values = {}
    for key in ("R_norm_inf_exact","lower_completion_minus_S_hat_norm_inf_exact","upper_completion_minus_S_hat_norm_inf_exact"):
        xi = value(key)/minimum
        values[key] = {"xi_estimate":scalar(np,xi),"profile_change_bounds":None}
        if xi<1:
            values[key]["profile_change_bounds"] = [scalar(np,np.longdouble(row["quadratic"])*xi/(2*(1-xi))) for row in receipt["rows"]]
    denominator = minimum-rnorm
    values["A_t_estimate"] = scalar(np,value("A_norm_inf_exact")/denominator) if denominator>0 else None
    receipt["derived_target_sensitivity"] = {"engineering_estimates":values,"source_uncertainty":None,
        "triangle_variants_executed":False,"raw_B_score_executed":False,"validated_eigenvalue_interval":False}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--runtime-fingerprint",action="store_true")
    parser.add_argument("--request");parser.add_argument("--request-sha");parser.add_argument("--assembly");parser.add_argument("--control",choices=controls.REFERENCE_CONTROLS)
    parser.add_argument("--output",required=True)
    args = parser.parse_args()
    out = Path(args.output);t.need(out.is_absolute() and not out.exists(), "fresh absolute reference directory")
    out.mkdir(mode=0o700)
    receipt = {"interface":"SN-S_hat-spectral-reference/v1","target_id":t.TARGET,"control":args.control,"status":"refused",
        "request":None,"assembly":None,"runtime_before":None,"runtime_after":None,"eigenvalues":None,"eigenvectors":None,
        "response_solution":None,
        "diagnostics":None,"rows":[],"derived_target_sensitivity":None,"error":None,"terminal_errors":[],
        "work":{"eigh_attempts":0,"eigh_completed":0,"dense_check_products_attempted":0,"dense_check_products_completed":0,"profile_attempts":0,"profile_completed":0},
        "normalized_density":None,"prior":None,"source_uncertainty":None,"joint_target":None,"error_certificate":None}
    start = time.monotonic();request = np = expected_runtime = None
    try:
        if args.runtime_fingerprint:
            t.need(not any((args.request,args.request_sha,args.assembly,args.control)), "exclusive metadata CLI")
        else:
            t.need(bool(args.control) != bool(args.assembly),"one scientific input/control branch")
            request,receipt["request"] = t.admitted_request(args.request,args.request_sha)
            t.need(request["reference_runtime"] is not None,"actual reference runtime pin missing")
            expected_runtime=t.document(t.read(request["reference_runtime"])[0])
        import numpy as np
        t.need(np.__version__=="2.2.6", "exact NumPy2.2.6")
        receipt["runtime_before"] = runtime_fingerprint(np)
        if args.runtime_fingerprint:
            receipt["status"] = "runtime-only"
        else:
            t.need(receipt["runtime_before"] == expected_runtime, "pinned complete runtime preflight")
            if args.control:
                c,residuals = controls.synthetic(args.control)
                c = np.asarray(c,dtype=np.float64)
            else:
                assembly_path = Path(args.assembly)
                raw,receipt["assembly"] = t.owned_json_bytes(assembly_path,2*1024**2)
                assembly = t.document(raw)
                t.need(assembly["status"]=="assembled-new-S_hat-and-four-mean-vectors-no-score" and assembly["request"]["sha256"]==args.request_sha, "new target assembly request")
                raw,_ = t.read(assembly["S_hat"])
                t.need(len(raw)==8*t.N*t.N,"S_hat exact byte shape")
                c = np.frombuffer(raw,dtype=">f8").astype(np.float64).reshape(t.N,t.N)
                residuals = [[float.fromhex(v) for v in row["residual_hex"]] for row in assembly["case_vectors"]]
            evaluate(np,c,residuals,out,receipt)
            if not args.control:
                ledger_raw,_ = t.read(assembly["assembly_ledger"])
                sensitivity(np,t.document(ledger_raw),receipt)
    except BaseException as exc:
        receipt["status"] = "refused";receipt["error"]=t.failure(exc)
    finally:
        if np is not None and receipt["runtime_before"] is not None:
            try:
                receipt["runtime_after"] = runtime_fingerprint(np)
                t.need(receipt["runtime_before"]==receipt["runtime_after"],"actual imported/mapped/runtime terminal drift")
                if request is not None:
                    t.need(t.read(receipt["request"],False)[1]==receipt["request"],"terminal request drift")
                    for p in request["source_ports"].values():t.read(p,False)
                    t.read(request["reference_runtime"],False)
            except BaseException as exc:
                receipt["status"]="refused";receipt["terminal_errors"].append(t.failure(exc))
        receipt["wall_seconds"] = time.monotonic()-start
        t.store_bytes(out)
        t.emit(out/"reference.json",t.encoded(receipt),1024**2)
        for p in out.iterdir():p.chmod(0o444)
        out.chmod(0o555)
    return 0 if receipt["status"] in ("accepted","runtime-only") else 2

if __name__=="__main__":
    sys.exit(main())
