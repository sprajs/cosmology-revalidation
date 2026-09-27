"""Independent released-matrix audit; general LU, Schur profile and measure checks.

Does not import the first-party ladder solver, change any likelihood, or export an H0 prior.
Run with OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import numpy as np
from astropy.io import fits
from scipy import integrate, linalg, stats


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def audit(root):
    work=root/'.work/unified-cosmology/distance-ladder'
    out=root/'studies/unified_cosmology/results/distance_ladder'
    code=root/'studies/unified_cosmology/code/distance_ladder'
    inputs=[work/('all'+c+'_shoes_ceph_topantheonwt6.0_112221.fits')for c in 'ylc']
    dependencies=inputs+[work/'README.md',work/'MCMC_utils.py',work/'run_mcmc.py',
        work/'read_chains_example.py',work/'lstsq_results.txt',work/'Riess2022v3.pdf',
        code/'solve.py',code/'common.py',code/'design.json',out/'acquisition.json',out/'reconstruction.json',
        work/'gls-result.npz']
    initial={str(p.relative_to(root)):sha(p)for p in dependencies}
    release=json.loads((out/'acquisition.json').read_text())
    for item in release['files']:
        p=root/item['path']
        assert p.stat().st_size==item['bytes'] and sha(p)==item['sha256']
    y,A,C=[np.asarray(fits.getdata(p),dtype=np.float64)for p in inputs]
    A=A.T
    assert y.shape==(3492,) and A.shape==(3492,47) and C.shape==(3492,3492)
    assert np.array_equal(C,C.T) and all(np.isfinite(x).all()for x in [y,A,C])
    # Deliberately do not Cholesky-whiten C or call the producer's solve/gls.
    lu=linalg.lu_factor(C)
    solved=linalg.lu_solve(lu,np.column_stack((y,A)))
    F=A.T@solved[:,1:]; b=A.T@solved[:,0]
    q=linalg.solve(F,b,assume_a='gen')
    V=linalg.inv(F)
    r=y-A@q
    chi2=float(r@linalg.lu_solve(lu,r))
    reference=json.loads((out/'reconstruction.json').read_text())
    published_q=np.array([v['fit']for v in reference['parameters_by_index']])
    published_V=np.load(work/'gls-result.npz')['covariance']
    fiterr=float(np.max(abs(q-published_q)))
    coverr=float(np.max(abs(V-published_V)))
    assert fiterr<1e-8 and coverr<1e-10 and abs(chi2-reference['validation']['chi2'])<1e-8
    # Marginalizing all other linear coordinates gives the same H0-coordinate variance.
    Fn=F[:-1,:-1]; cross=F[:-1,-1]
    inverse_cross=linalg.solve(Fn,cross,assume_a='gen')
    k=float(F[-1,-1]-F[-1,:-1]@inverse_cross)
    rhs=float(b[-1]-F[-1,:-1]@linalg.solve(Fn,b[:-1],assume_a='gen'))
    profile_q=rhs/k;profile_variance=1/k
    assert abs(profile_q-q[-1])<1e-9 and abs(profile_variance-V[-1,-1])<1e-12
    curvature=[]
    for z in [-3.,-2.,-1.,0.,1.,2.,3.]:
        h=q[-1]+z*np.sqrt(V[-1,-1])
        nuisance=linalg.solve(Fn,b[:-1]-cross*h,assume_a='gen')
        residual=y-A@np.r_[nuisance,h]
        delta=float(residual@linalg.lu_solve(lu,residual)-chi2)
        curvature.append({'standard_deviation_offset':z,'delta_chi2':delta,'expected':z*z})
    assert max(abs(v['delta_chi2']-v['expected'])for v in curvature)<1e-8
    # Author README: native slope is a delta from -3.285.
    y_physical=y-3.285*A[:,41]
    q_physical=linalg.solve(F,A.T@linalg.lu_solve(lu,y_physical),assume_a='gen')
    expected=q.copy();expected[41]-=3.285
    offset_error=float(np.max(abs(q_physical-expected)))
    assert offset_error<1e-8
    # HF sign is MB - 5log10(H0). Calibrator rows are mu_host + MB.
    hf=A[:,46]!=0;sn=A[:,42]!=0;cal=sn&~hf;cep=~sn
    pattern=np.zeros(47);pattern[42]=1;pattern[46]=-1
    assert np.array_equal(A[hf],np.repeat(pattern[None,:],int(hf.sum()),axis=0))
    assert np.all(A[cal,42]==1) and np.all(np.count_nonzero(A[cal],axis=1)==2)
    assert np.all(A[cal,:37].sum(axis=1)==1)
    assert np.all(A[cal,:37]>=0) and np.all(A[cal,:37]<=1)
    # The zero-width author box applies to one algebraically isolated auxiliary only.
    auxiliary_rows=np.flatnonzero(A[:,44])
    assert auxiliary_rows.tolist()==[3210]
    row=3210
    assert y[row]==0 and np.count_nonzero(A[row])==1 and A[row,44]==1
    assert np.count_nonzero(C[row])==1
    mask=np.arange(3492)!=row;cols=np.arange(47)!=44
    Cm=C[np.ix_(mask,mask)];Am=A[np.ix_(mask,cols)];ym=y[mask]
    reduced=linalg.solve(Cm,np.column_stack((ym,Am)),assume_a='gen')
    Fm=Am.T@reduced[:,1:]
    qm=linalg.solve(Fm,Am.T@reduced[:,0],assume_a='gen')
    auxerr=float(np.max(abs(qm-q[cols])))
    assert auxerr<1e-8 and abs(V[44,44]-C[row,row])<1e-20
    # Positive-transform Jacobian and moments, independent numerical integration in H0.
    mean_q=q[-1]; sd_q=np.sqrt(V[-1,-1]);a=np.log(10)/5
    lo,hi=np.exp(a*(mean_q+sd_q*np.array([-10,10])))
    def density(h):
        return stats.norm.pdf((5*np.log10(h)-mean_q)/sd_q)/sd_q*5/(np.log(10)*h)
    numeric=[]
    for power in [0,1,2]:
        numeric.append(integrate.quad(lambda h:h**power*density(h),lo,hi,epsabs=1e-10,epsrel=1e-12)[0])
    analytic=[1.,np.exp(a*mean_q+.5*a*a*sd_q*sd_q),np.exp(2*a*mean_q+2*a*a*sd_q*sd_q)]
    moment_error=float(np.max(abs(np.array(numeric)/analytic-1)))
    assert moment_error<1e-10
    probabilities=np.array([.025,.16,.5,.84,.975]);quantiles=np.exp(a*(mean_q+sd_q*stats.norm.ppf(probabilities)))
    integrated_cdf=[integrate.quad(density,lo,x,epsabs=1e-12,epsrel=1e-12)[0]for x in quantiles]
    cdferr=float(np.max(abs(np.array(integrated_cdf)-probabilities)))
    assert cdferr<1e-10
    means,uncertainties=np.loadtxt(work/'lstsq_results.txt',unpack=True)
    j=np.arange(47)!=44
    lower=(means[j]-10*uncertainties[j]-q[j])/np.sqrt(np.diag(V)[j])
    upper=(means[j]+10*uncertainties[j]-q[j])/np.sqrt(np.diag(V)[j])
    discarded_bound=float(np.sum(stats.norm.cdf(lower)+stats.norm.sf(upper)))
    assert discarded_bound<1e-12
    cross_C=C[np.ix_(np.flatnonzero(cal),np.flatnonzero(hf))]
    corr=cross_C/np.sqrt(np.diag(C)[cal,None]*np.diag(C)[None,hf])
    assert np.count_nonzero(cross_C)==77*277
    assert np.count_nonzero(C[np.ix_(np.flatnonzero(cep),np.flatnonzero(sn))])==0
    singleton=[]
    for row in np.flatnonzero(np.count_nonzero(A,axis=1)==1):
        singleton.append({'row':int(row),'column':int(np.flatnonzero(A[row])[0]),'y':float(y[row]),
            'sigma':float(np.sqrt(C[row,row])),'other_covariance_entries':int(np.count_nonzero(C[row])-1)})
    assert {str(p.relative_to(root)):sha(p)for p in dependencies}==initial
    return {'status':'passed_independent_released_matrix_review','calls_to_producer_solver':0,
        'source_sha256':sha(__file__),'dependency_sha256':initial,
        'versions':{n:importlib.metadata.version(n)for n in ['numpy','scipy','astropy']},
        'method':'General LU solve C against [y,A]; general linear solve/inverse of normal system; scalar Schur elimination and direct residual profiling. No producer code imported.',
        'rows':3492,'parameters':47,'df':3445,'full_covariance_retained':True,
        'H0_median':float(10**(q[-1]/5)),'H0_propagated_sd':float(a*10**(q[-1]/5)*sd_q),
        'H0_mean':float(analytic[1]),'H0_equal_tail_quantiles':quantiles.tolist(),
        'H0_quantile_probabilities':probabilities.tolist(),'H0_density_measure':'Flat native linear q; H0 density includes 5/(ln(10)*H0). Not a flat-H0 posterior or an exportable independent H0 prior.',
        'chi2':chi2,'parameter_max_absolute_difference':fiterr,'covariance_max_absolute_difference':coverr,
        'Schur_mean_difference':float(profile_q-q[-1]),'Schur_variance_difference':float(profile_variance-V[-1,-1]),
        'profile_curvature':curvature,'physical_slope':float(q_physical[41]),'slope_reparameterization_error':offset_error,
        'auxiliary_removed_max_physical_parameter_change':auxerr,'H0_numerical_moment_relative_error':moment_error,
        'H0_numerical_CDF_max_error':cdferr,'author_broad_box_discarded_probability_union_bound_excluding_isolated_auxiliary':discarded_bound,
        'minimum_author_bound_distance_in_marginal_sd':float(min(np.min(-lower),np.min(upper))),
        'singleton_embedded_constraints':singleton,
        'row_counts':{'Cepheid_and_external':int(cep.sum()),'SN_calibrator':int(cal.sum()),'SN_Hubble_flow':int(hf.sum())},
        'calibrator_HF_covariance_nonzero_entries':int(np.count_nonzero(cross_C)),
        'calibrator_HF_max_absolute_correlation':float(np.max(abs(corr))),
        'conditional_release_Cepheid_to_SN_covariance_zero':True,
        'physical_review':[
          'HF matrix uses MB - 5log10H0; calibrators use mu_host + MB. Follow actual matrix and paper Eq3; paper Eq16 reverses the MB sign in its printed expression.',
          'Author README explicitly supplies the -3.285 slope offset; all external Gaussian constraints are already rows in released y/A/C and are not added again.',
          'Author MCMC zero-width coordinate44 fixes an isolated auxiliary. Removing it leaves physical fit unchanged; flat-q analytic reconstruction is not a rerun of author MCMC.',
          'Released baseline Hubble-flow y includes the adopted low-redshift expansion, q0=-0.55 per paper section5.2. This is not an unrestricted expansion-history-independent H0 measurement.',
          'Table5 prints metallicity-slope sigma0.046 whereas matrix yields '+str(float(np.sqrt(V[43,43])))+'; preserve this small release/table discrepancy rather than claiming every uncertainty rounds identically.',
          'The covariance includes within-release shared calibrator/Hubble-flow uncertainties. Existing Dovekie/Pantheon events and cross-release calibration cannot be made independent by renaming this H0 summary.',
          'Cepheid-only factor extraction may exploit released zero Cepheid-to-SN covariance, conditional on that release model. It still needs consistent host mapping and common SN standardization before joint cosmology.'
        ],'source_references':[
          'https://arxiv.org/abs/2112.04510v3',
          'https://github.com/PantheonPlusSH0ES/DataRelease/tree/c447f0fea703fcd0fff57de5000947b5ca81286b/SH0ES_Data'],
        'qualification':'Independent numerical/physical audit of released compressed products only; not a raw-photometry reconstruction, new cosmology result, or permission to multiply overlapping likelihoods.'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path.cwd());p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();result=audit(args.root.resolve());args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k]for k in ['status','H0_median','H0_propagated_sd','parameter_max_absolute_difference','covariance_max_absolute_difference']}))
if __name__=='__main__':main()
