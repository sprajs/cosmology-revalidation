"""Reference-only full-C LAPACK SVD/QR; not a production inference route."""
import os
for name in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","BLIS_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[name]="1"
import contextlib
import hashlib
import io
import json
from pathlib import Path
import sys
import numpy as np
import scipy
import scipy.linalg as la

def library_identity():
    files=set()
    for package in (np,scipy):
        root=Path(package.__file__).parent
        files.update(root.rglob("*.so"))
        sibling=root.parent/(root.name+".libs")
        if sibling.is_dir():files.update(p for p in sibling.iterdir() if p.is_file())
    return {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}

libraries_before=library_identity()

folder=Path(sys.argv[1]);n,p=3492,47
fixed=len(sys.argv)==3
if fixed and sys.argv[2]!="released-fixed44/v1":raise ValueError("unknown reference target")
target=json.loads((folder/"constrained.json").read_text()) if fixed else None
m=json.loads((folder/"lineage.json").read_text())
for name,digest in m["decode"]["canonical_sha256"].items():
    if hashlib.sha256((folder/name).read_bytes()).hexdigest()!=digest:raise ValueError("reference input hash differs")
c=np.fromfile(folder/"C.f64",dtype="<f8").reshape(n,n)
x=np.fromfile(folder/"X.f64",dtype="<f8").reshape(n,p)
active=[j for j in range(47) if not fixed or j!=44]
if fixed:x=x[:,active];p=46
y=np.fromfile(folder/"y.f64",dtype="<f8")
if not np.array_equal(c,c.T):raise ValueError("asymmetric source")
factor=la.cholesky(c,lower=True)
a=la.solve_triangular(factor,x,lower=True)
z=la.solve_triangular(factor,y,lower=True)
u,s,vt=la.svd(a,full_matrices=False,lapack_driver="gesdd")
if np.count_nonzero(s>np.finfo(float).eps*max(a.shape)*s[0])!=p:raise ValueError("full47 SVD rank not admitted")
svd=(vt.T/s)@(u.T@z)
q,r,pivot=la.qr(a,mode="economic",pivoting=True)
qr_p=la.solve_triangular(r,q.T@z)
qr=np.empty(p);qr[pivot]=qr_p
w=np.zeros(p);w[-1]=1
sv=(vt@w)/s
qv=la.solve_triangular(r,w[pivot],trans="T")
def quadratic(beta,observed):
    residual=observed-x@beta
    white=la.solve_triangular(factor,residual,lower=True)
    return float(white@white)
report={"software":{"python":sys.version,"numpy":np.__version__,"scipy":scipy.__version__},
        "ancestry":"SVD and QR share exact inputs and LAPACK Cholesky whitening; native portable wide factor/QR is separate implementation ancestry. No raw-data or physical independence implied.",
        "rank":p,"condition2":float(s[0]/s[-1]),"algorithms":{
          "LAPACK_gesdd_SVD":{"coefficients":np.insert(svd,44,0).tolist() if fixed else svd.tolist(),"quadratic":quadratic(svd,y),"variance46":float(sv@sv)},
          "LAPACK_pivoted_QR":{"coefficients":np.insert(qr,44,0).tolist() if fixed else qr.tolist(),"quadratic":quadratic(qr,y),"variance46":float(qv@qv)}},"sensitivities":[]}
# Full47 retained SVD law with synthetic perturbed observations only.
variations=target["calibration_sensitivities"] if fixed else [{"row":i,"delta_y":.01} for i in range(3207,3215)]
for variation in variations:
    row,delta=variation["row"],variation["delta_y"]
    item={"row":row,"delta_y":delta}
    for sign,label in ((1,"plus"),(-1,"minus")):
        observed=y.copy();observed[row]+=sign*delta
        whitened=la.solve_triangular(factor,observed,lower=True)
        beta=(vt.T/s)@(u.T@whitened)
        item["beta46_"+label]=float(beta[-1]);item["q_"+label]=quadratic(beta,observed)
    report["sensitivities"].append(item)
config=io.StringIO()
with contextlib.redirect_stdout(config):
    np.show_config();scipy.show_config()
report["numeric_configuration"]=config.getvalue()
report["reference_executable_sha256"]=hashlib.sha256(Path(sys.executable).read_bytes()).hexdigest()
report["thread_environment"]={name:os.environ[name] for name in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","BLIS_NUM_THREADS","NUMEXPR_NUM_THREADS")}
report["libraries_before"]=libraries_before;report["libraries_after"]=library_identity()
if report["libraries_after"]!=libraries_before:raise ValueError("reference libraries changed during execution")
print(json.dumps(report))
