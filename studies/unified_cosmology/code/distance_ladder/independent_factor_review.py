"""Independent SN-free host-factor review using orthogonal nuisance projection.

No producer modules imported; no new cosmological prior or target is exported.
"""
import argparse
from collections import defaultdict
from decimal import Decimal
import hashlib
import importlib.metadata
import json
from pathlib import Path
import numpy as np
from astropy.io import fits
from scipy import linalg


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(root):
    work=root/'.work/unified-cosmology/distance-ladder'
    code=root/'studies/unified_cosmology/code/distance_ladder'
    out=root/'studies/unified_cosmology/results/distance_ladder'
    paths=[work/('all'+c+'_shoes_ceph_topantheonwt6.0_112221.fits')for c in 'ylc']
    deps=paths+[work/'table2.tex',work/'cepheid-host-factor.npz',out/'cepheid-factor.json',
        code/'cepheid_factor.py',code/'cepheid-design.json',code/'solve.py',out/'acquisition.json',out/'reconstruction.json']
    hashes={str(p.relative_to(root)):sha(p)for p in deps}
    acquisition=json.loads((out/'acquisition.json').read_text())
    for row in acquisition['files']:
        p=root/row['path']
        assert p.stat().st_size==row['bytes'] and sha(p)==row['sha256']
    producer=json.loads((out/'cepheid-factor.json').read_text())
    for key,digest in producer['source_sha256'].items():assert sha(root/key)==digest
    for key,digest in producer['dependencies_sha256'].items():assert sha(root/key)==digest
    for key,digest in producer['output_sha256'].items():assert sha(root/key)==digest
    y,A,C=[np.asarray(fits.getdata(p),dtype=float)for p in paths];A=A.T
    rows=np.flatnonzero(A[:,42]==0)
    cols=np.r_[np.arange(42),43,44,45]
    discarded=np.flatnonzero(A[:,42]!=0)
    assert len(rows)==3138 and len(discarded)==354 and len(cols)==45
    assert np.all(A[rows,42]==0) and np.all(A[rows,46]==0)
    assert np.count_nonzero(C[np.ix_(rows,discarded)])==0
    # Marginalize nuisance coordinates in data space with an orthogonal projection,
    # instead of the producer's 45D fit/covariance slice or normal-matrix Schur solve.
    L=linalg.cholesky(C[np.ix_(rows,rows)],lower=True)
    T=linalg.solve_triangular(L,A[np.ix_(rows,cols)],lower=True)
    d=linalg.solve_triangular(L,y[rows],lower=True)
    H,N=T[:,:37],T[:,37:]
    Q,R=linalg.qr(N,mode='economic')
    assert np.linalg.matrix_rank(R)==8
    Hp=H-Q@(Q.T@H);dp=d-Q@(Q.T@d)
    U,s,Vh=linalg.svd(Hp,full_matrices=False,lapack_driver='gesvd')
    assert np.all(s>np.finfo(float).eps*max(Hp.shape)*s[0])
    mean=Vh.T@((U.T@dp)/s)
    cov=(Vh.T/s**2)@Vh
    released=np.load(work/'cepheid-host-factor.npz')
    meanerr=float(np.max(abs(mean-released['mean_mu'])))
    coverr=float(np.max(abs(cov-released['covariance_mu'])))
    assert meanerr<1e-8 and coverr<1e-10
    assert np.max(abs(mean-np.array(producer['mean_distance_modulus'])))<1e-8
    assert np.max(abs(cov-np.array(producer['covariance_distance_modulus'])))<1e-10
    residual=dp-Hp@mean
    chi2=float(residual@residual)
    assert abs(chi2-producer['validation']['chi2'])<1e-8
    # Deterministic density-shape checks, independent of producer's random probes.
    factors=linalg.cholesky(cov,lower=True)
    curvature=[]
    for k in range(8):
        v=np.cos(np.arange(37)*(k+1)/7)
        v=v/np.linalg.norm(v)*(k+1)/2
        delta=factors@v
        r=dp-Hp@(mean+delta)
        actual=float(r@r-chi2)
        curvature.append({'offset_norm':float(np.linalg.norm(v)),
            'projected_delta_chi2':actual,'expected':float(v@v)})
    curve_error=max(abs(v['projected_delta_chi2']-v['expected'])for v in curvature)
    assert curve_error<1e-8
    # Independent public-table parsing and candidate search, with decimal bounds.
    groups=defaultdict(list)
    for line in (work/'table2.tex').read_text().splitlines():
        if '&' in line and not line.lstrip().startswith('\\'):
            parts=line.split('&')
            groups[parts[0].strip()].append([x.strip().removesuffix('\\\\').strip()for x in parts])
    assert sum(len(x)for x in groups.values())==3130
    labels=[]
    for j in range(37):
        ix=np.flatnonzero((A[:,j]!=0)&(A[:,38]!=0))
        matching=[]
        for name,table in groups.items():
            if len(table)!=len(ix):continue
            valid=True
            for row,item in zip(ix,table):
                for field,column,log_period in [(4,41,True),(9,43,False)]:
                    val=Decimal(item[field]);half=Decimal(5).scaleb(val.as_tuple().exponent-1)
                    lo,hi=float(val-half),float(val+half)
                    if log_period:lo,hi=np.log10([lo,hi])-1
                    native=A[row,column]
                    tolerance=2*abs(float(np.spacing(np.float32(native))))
                    if not lo-tolerance<=native<=hi+tolerance:valid=False;break
                if not valid:break
            if valid:matching.append(name)
        assert len(matching)==1,(j,matching)
        label=matching[0]
        assert label==str(released['host'][j])==producer['host_order'][j]
        assert producer['host_mapping'][j]['released_row_indices_zero_based']==ix.tolist()
        labels.append({'column':j,'host':label,'Cepheid_rows':len(ix),'unique_candidates':1})
    assert len({x['host']for x in labels})==37
    # Quantify the necessary distinction from hosts fitted using SN data.
    full=json.loads((out/'reconstruction.json').read_text())
    full_means=np.array([x['fit']for x in full['parameters_by_index'][:37]])
    shift=float(np.max(abs(mean-full_means)))
    assert abs(shift-producer['maximum_host_mean_change_after_removing_SN_mag'])<1e-8
    final={str(p.relative_to(root)):sha(p)for p in deps};assert final==hashes
    return {'status':'passed_independent_SN_free_factor_review','source_sha256':sha(Path(__file__)),
        'dependency_sha256':hashes,'versions':{p:importlib.metadata.version(p)for p in ['numpy','scipy','astropy']},
        'method':'Full covariance whitening, QR orthogonal nuisance projection, independent GESVD of projected 37-host design. No producer module imports.',
        'retained_rows':len(rows),'removed_SN_rows':len(discarded),'nuisance_rank':8,'host_rank':37,
        'mean_max_absolute_difference_mag':meanerr,'covariance_max_absolute_difference_mag2':coverr,
        'chi2':chi2,'chi2_difference':float(chi2-producer['validation']['chi2']),
        'profile_curvature_checks':curvature,'profile_max_absolute_error':curve_error,
        'host_labels':labels,'removed_retained_covariance_nonzero_entries':0,
        'maximum_host_change_from_full_ladder_mag':shift,
        'qualification':'Correct conditional SN-free compression under the released Gaussian model with fixed design/covariance and flat linear nuisance measure. All eight nuisance coordinates are integrated; the host covariance is full. Not an independent reconstruction of photometry, a progenitor-age result, a Gaussian H0 prior, or a ready-to-multiply factor for unmatched/restandardized SN.',
        'integration_requirements':['One physical host/SN linkage, preserving duplicate light-curve rows and aliases.',
            'Consistent calibrated SN absolute-magnitude standardization and its covariance; shared cross-release calibration response must not be ignored.',
            'Preserve correlated37-host covariance and compressed anchor assumptions; do not add anchor priors a second time.']}


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path.cwd());p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();result=audit(args.root.resolve());args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k]for k in ['status','mean_max_absolute_difference_mag','covariance_max_absolute_difference_mag2','profile_max_absolute_error']}))
if __name__=='__main__':main()
