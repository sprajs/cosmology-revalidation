"""Independent multistart maximum-likelihood comparison of alternative SN releases.

Each SN compilation is fitted separately with the same DESI DR2 BAO likelihood.
No posterior sampling, CMB, age correction, or inter-release multiplication.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import datetime, hashlib, json, time
import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import differential_evolution, minimize
from late_geometry import Geometry, ROOT, BAO, sample_path

OUT=ROOT/'studies/unified_cosmology/results/inference/release-comparison.json'
WORK=ROOT/'.work/unified-cosmology/inference/release-comparison'
SEEDS=[927701,927702,927703]
SAMPLES=['dovekie','pantheon','des3yr']
MODELS=['lcdm','wcdm','cpl']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def pantheon_audit():
    directory=sample_path('pantheon').parent
    table=pd.read_csv(directory/'Pantheon+SH0ES.dat',sep=r'\s+',dtype={'CID':str})
    raw=np.loadtxt(directory/'Pantheon+SH0ES_STAT+SYS.cov')
    assert int(raw[0])==len(table) and len(raw)==1+len(table)**2
    original=raw[1:].reshape(len(table),len(table));cov=(original+original.T)/2
    keep=table.zHD.gt(.01).to_numpy();selected=table.loc[keep].reset_index(drop=True)
    with np.load(sample_path('pantheon')) as f:
        assert np.array_equal(f['CID'],selected.CID.to_numpy())
        for field,column in [('zHD','zHD'),('zHEL','zHEL'),('MU','m_b_corr'),('IDSURVEY','IDSURVEY')]:
            assert np.array_equal(f[field],selected[column].to_numpy()),field
        assert np.array_equal(f['covariance'],cov[np.ix_(keep,keep)])
        inverse_error=float(np.max(abs(f['precision']@f['covariance']-np.eye(len(selected)))))
        assert inverse_error<1e-10
        sub=f['covariance'];counts=selected.CID.value_counts()
        corr=[]
        for name,count in counts.items():
            if count<2:continue
            ix=np.flatnonzero(selected.CID.to_numpy()==name)
            for a,j in enumerate(ix):
                for k in ix[a+1:]:corr.append(sub[j,k]/np.sqrt(sub[j,j]*sub[k,k]))
    return {'raw_rows':len(table),'retained_zHD_gt_0_01_rows':len(selected),'distinct_released_CID_names':int(counts.size),'excess_repeat_measurement_rows':int(len(selected)-len(counts)),'repeated_CID_objects':int((counts>1).sum()),'measurement_multiplicity_distribution':{str(k):int(v)for k,v in counts.value_counts().items()},'same_CID_covariance_correlation_min_median_max':[float(x)for x in [min(corr),np.median(corr),max(corr)]],'retained_IS_CALIBRATOR_rows':int(selected.IS_CALIBRATOR.sum()),'row_order_redshift_magnitude_survey_and_full_covariance_exactly_match_raw_release':True,'inverse_identity_max_error':inverse_error,'original_covariance_max_asymmetry':float(np.max(abs(original-original.T))),'semantics':'Rows are released light-curve measurements, not all independent physical SNe. Repeated CID rows and their full released covariance are retained. m_b_corr with one free global offset is used; neither MU_SH0ES nor CEPH_DIST is a likelihood observation. No absolute Cepheid anchor.'}


def independent_components(g,theta):
    om,scale,w,wa,_=[float(x[0])for x in g.physical(theta)]
    def inv_e(z):return 1/np.sqrt(om*(1+z)**3+(1-om)*(1+z)**(3*(1+w+wa))*np.exp(-3*wa*z/(1+z)))
    ratio=np.array([quad(inv_e,0,z,epsabs=2e-12,epsrel=2e-12)[0]/z for z in g.sz])
    shape=5*np.log10(g.sz*(1+g.zhel)*ratio)
    raw=g.observed-shape
    chol=cho_factor(g.cov,lower=True);p1=cho_solve(chol,np.ones(len(raw)));pr=cho_solve(chol,raw)
    offset=float(pr.sum()/p1.sum());res=raw-offset
    schi=float(res@cho_solve(chol,res))
    ib=np.array([quad(inv_e,0,z,epsabs=2e-12,epsrel=2e-12)[0]for z in g.bz])
    dh=299792.458/scale*np.array([inv_e(z)for z in g.bz]);dm=299792.458/scale*ib
    predictions=np.where(g.btype=='DH_over_rs',dh,np.where(g.btype=='DM_over_rs',dm,np.cbrt(g.bz*dh*dm**2)))
    residual=predictions-g.bmean;bchi=float(residual@g.bprecision@residual)
    return schi,bchi,offset


def fit(g):
    low,high=np.array(g.bounds).T;span=high-low
    def convert(x):return low+span*np.asarray(x)
    def objective(x):return sum(float(v[0])for v in g.components(convert(x)))
    bounds=[(1e-10,1-1e-10)]*len(low)
    starts=[];t=time.monotonic()
    for seed in SEEDS:
        d=differential_evolution(objective,bounds,seed=seed,popsize=16,maxiter=900,tol=1e-10,atol=1e-8,polish=False,workers=1)
        p=minimize(objective,d.x,method='Nelder-Mead',bounds=bounds,options={'maxiter':5000,'xatol':2e-10,'fatol':1e-9})
        best=p if p.fun<=d.fun else d
        starts.append({'method':'differential_evolution_then_Nelder_Mead','seed':seed,'success':bool(d.success and p.success),'de_success':bool(d.success),'local_success':bool(p.success),'de_message':str(d.message),'local_message':str(p.message),'theta':convert(best.x).tolist(),'chi2':float(best.fun),'de_evaluations':int(d.nfev),'local_evaluations':int(p.nfev)})
    # Additional local searches begin away from the global starts, with scaled coordinates.
    rng=np.random.default_rng(927704)
    for i in range(8):
        initial=rng.uniform(.05,.95,len(low))
        p=minimize(objective,initial,method='Nelder-Mead',bounds=bounds,options={'maxiter':5000,'xatol':2e-9,'fatol':1e-8})
        starts.append({'method':'independent_Nelder_Mead','start':i,'success':bool(p.success),'initial_theta':convert(initial).tolist(),'theta':convert(p.x).tolist(),'chi2':float(p.fun),'evaluations':int(p.nfev),'message':str(p.message)})
    best=min(starts,key=lambda r:r['chi2']);theta=np.array(best['theta']);native=[float(x[0])for x in g.components(theta)]
    direct=independent_components(g,theta)
    assert max(abs(np.array(native)-np.array(direct[:2])))<1e-4
    de=[s['chi2']for s in starts if s['method'].startswith('differential')]
    assert max(de)-min(de)<1e-4,(g.model,de)
    om,_,w,_,_=[float(x[0])for x in g.physical(theta)]
    boundary=np.minimum((theta-low)/span,(high-theta)/span)
    return {'parameter_names':g.names,'bounds':g.bounds,'theta':theta.tolist(),'chi2_total':sum(native),'chi2_SN':native[0],'chi2_BAO':native[1],'q0':.5*om+.5*(1+3*w)*(1-om),'profiled_SN_global_offset_in_adapter_units':direct[2],'independent_SN_BAO_chi2':list(direct[:2]),'independent_max_abs_chi2_difference':float(max(abs(np.array(native)-np.array(direct[:2])))),'de_multistart_chi2_spread':max(de)-min(de),'parameter_bound_min_fraction':boundary.tolist(),'parameters_within_1e_5_fraction_of_boundary':[n for n,b in zip(g.names,boundary)if b<1e-5],'starts_within_1e_4_of_best':sum(s['chi2']-best['chi2']<1e-4 for s in starts),'all_starts':starts,'elapsed_seconds':time.monotonic()-t}


def main():
    WORK.mkdir(parents=True,exist_ok=True)
    frozen=Path(__file__).with_name('late_geometry.py');before=sha(frozen)
    design={'samples_separately':SAMPLES,'models':MODELS,'evolution':'none','CMB':False,'spatial_curvature':0,'radiation':'omitted, matching frozen late-time geometry','early_dark_energy_or_w0_plus_wa_cut':False,'seeds':SEEDS,'global_starts_per_fit':3,'additional_local_starts_per_fit':8,'bounds':'unchanged Geometry bounds; normalized optimizer coordinates, no posterior prior interpretation','direct_validation':'Uncompressed per-row adaptive distance integration and explicit scalar-offset profiling; independent BAO integration','estimand':'Maximum likelihood conditional on fixed released covariance and corrected distances; optimizer convergence is not posterior convergence.'}
    (WORK/'execution-design.json').write_text(json.dumps(design,indent=2)+'\n')
    audit=pantheon_audit();fits={};interfaces={}
    for sample in SAMPLES:
        fits[sample]={}
        for model in MODELS:
            g=Geometry(model=model,evolution='none',sample=sample)
            result=fit(g);fits[sample][model]=result
            print(json.dumps({'sample':sample,'model':model,'theta':result['theta'],'chi2':result['chi2_total'],'global_spread':result['de_multistart_chi2_spread'],'near_boundary':result['parameters_within_1e_5_fraction_of_boundary']}),flush=True)
            if model=='lcdm':
                interfaces[sample]={'likelihood_rows':len(g.sz),'released_data_variable':'m_b_corr' if sample=='pantheon' else 'MU adapter convention','redshift_range':[float(min(g.sz)),float(max(g.sz))],'SN_logdet_covariance':float(np.linalg.slogdet(g.cov)[1]),'BAO_rows':len(g.bz),'BAO_redshift_range':[float(min(g.bz)),float(max(g.bz))],'covariance_positive_definite':True,'physical_events_semantics':'20 bin measurements represent329SNe;18nonempty,2released999mag-error sentinels' if sample=='des3yr' else ('1590measurements of1473releasedCID names' if sample=='pantheon' else '1820distance rows; use physical-event registry for inter-release overlap')}
        l=fits[sample]['lcdm']['chi2_total'];w=fits[sample]['wcdm']['chi2_total'];c=fits[sample]['cpl']['chi2_total']
        assert c<=w+1e-4 and w<=l+1e-4
        fits[sample]['within_sample_improvements']={'LCDM_minus_wCDM':l-w,'LCDM_minus_CPL':l-c,'wCDM_minus_CPL':w-c,'meaning':'Differences of minimum chi-square on the same data/covariance. Not posterior model probability, Bayes factor, or calibrated discovery significance.'}
        (WORK/f'{sample}.json').write_text(json.dumps(fits[sample],indent=2)+'\n')
    assert sha(frozen)==before,'Frozen geometry changed during computation.'
    inputs=[sample_path(s)for s in SAMPLES]+[BAO/'desi_gaussian_bao_ALL_GCcomb_mean.txt',BAO/'desi_gaussian_bao_ALL_GCcomb_cov.txt']
    inputs += [sample_path('pantheon').parent/n for n in ['Pantheon+SH0ES.dat','Pantheon+SH0ES_STAT+SYS.cov','README']]
    result={'schema':'late-geometry-release-comparison-v1','status':'passed','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'design':design,'fits':fits,'interfaces':interfaces,'pantheon_row_covariance_audit':audit,'likelihood_normalization':'Reported chi2 = profiled SN residual quadratic + BAO quadratic. Fixed logdet(C), Gaussian2pi terms and offset-integration constants omitted. They cancel for parameter/model comparisons on one unchanged dataset, but raw chi2 across different releases or binned/unbinned data is not a relative model likelihood or goodness comparison.','shared_BAO_covariance_used_once_per_alternative':True,'SN_releases_combined':False,'sample_overlap_warning':'These alternatives share supernova events/calibration and the identical BAO data; best-fit differences are not independent tension measurements.','code_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),frozen]},'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in inputs},'design_sha256':sha(WORK/'execution-design.json')}
    OUT.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'output':str(OUT.relative_to(ROOT))},indent=2))


if __name__=='__main__':main()
