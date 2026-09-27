#!/usr/bin/env python3
"""Conditional observational tests of the age-standardization dispute.

Inputs remain read-only. Large row tables are regenerated in --work. This is
not a replacement for a selection-aware flux likelihood or joint host posterior.
"""
from pathlib import Path
import argparse, datetime, hashlib, json, platform
import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.linalg import cho_factor, cho_solve, cholesky, solve_triangular
from scipy.optimize import minimize
from astropy.io import ascii

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent
DEFAULT_ARCHIVE = Path.home()/'.local/share/cosmology-revalidation/archive/17487bf659fcbdeeea072221492bac14b04a0a85'

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def serial(x):
    if isinstance(x, np.ndarray): return x.tolist()
    if isinstance(x, np.generic): return x.item()
    raise TypeError(type(x).__name__)

def save(path, item):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(item, indent=2, default=serial, allow_nan=False)+'\n')

def read_age(path, source):
    rows=[]
    for row, line in enumerate(path.read_text().splitlines()[1:], 2):
        c=[s.replace('\\','').strip() for s in line.split('&')]
        if len(c)!=5: continue
        hr=c[3].split('(')
        rows.append([str(int(c[0])),source,row,float(c[1]),float(c[2]),float(hr[0]),float(hr[-1].strip(') ')),float(c[4])])
    return pd.DataFrame(rows,columns=['CID','age_source','age_source_line','age','age_err','hr_C25','hr_original','hr_C25_err'])

def mu(z, zhel, om=.3):
    radial=np.array([quad(lambda t:1/np.sqrt(om*(1+t)**3+1-om),0,zz,epsabs=1e-11)[0] for zz in z])
    return 5*np.log10((1+np.asarray(zhel))*299792.458/70*radial)+25

def gls(X,y,C):
    ci=cho_factor(C, lower=True,check_finite=False)
    cix=cho_solve(ci,X,check_finite=False)
    V=np.linalg.inv(X.T@cix)
    L=V@cix.T
    beta=L@y
    r=y-X@beta
    return beta,V,L,float(r@cho_solve(ci,r,check_finite=False))

def summary_gls(X,y,C,labels):
    b,v,_,chi=gls(X,y,C)
    return {'coefficients':dict(zip(labels,b)), 'standard_errors':dict(zip(labels,np.sqrt(np.diag(v)))), 'age_interval_normal_95': [b[1]-1.96*np.sqrt(v[1,1]),b[1]+1.96*np.sqrt(v[1,1])], 'chi2':chi,'dof':len(y)-len(b)}

def young_operator(age,z,error,folds=None,counts=None, weighted=True):
    n=len(age); counts=np.ones(n) if counts is None else counts
    Z=np.c_[np.ones(n),z]
    weight=counts/(error**2 if weighted else np.ones(n))
    H=np.eye(n); coefs=[]
    masks=[np.ones(n,dtype=bool)] if folds is None else [folds==f for f in range(5)]
    for test in masks:
        train=(age<4)&(counts>0)
        if folds is not None:train &= ~test
        ix=np.flatnonzero(train); zz=Z[ix]; w=weight[ix]
        B=np.linalg.solve(zz.T@(w[:,None]*zz),zz.T*w)
        H[np.ix_(np.flatnonzero(test),ix)] -= Z[test]@B
        coefs.append((ix,B))
    return H,coefs

def young_test(d, rng, work, name):
    age=d.age.to_numpy(); z=d.z.to_numpy(); y=d.hr_original.to_numpy(); e=d.hr_C25_err.to_numpy(); n=len(d)
    fold=np.array([int(hashlib.sha256(('2026092701:'+str(i)).encode()).hexdigest()[:8],16)%5 for i in d.CID])
    X=np.c_[np.ones(n),age]; C=np.diag(e*e)
    b,v,L,_=gls(X,y,C); H,co=young_operator(age,z,e); Hcf,cf=young_operator(age,z,e,fold)
    # The intercept is singular after full fitted trend removal. Use the original
    # age estimator and propagate H C H^T; do not invert that singular covariance.
    new=L@(H@y); cross=L@(Hcf@y)
    vin=L@H@C@H.T@L.T; vcf=L@Hcf@C@Hcf.T@L.T
    joint=summary_gls(np.c_[X,z],y,C,['intercept','age','z'])
    young=age<4
    published_diff=d.hr_C25.to_numpy()-y
    coefficient=np.linalg.lstsq(np.c_[np.ones(n),z],published_diff,rcond=None)[0]
    boot=[]
    for _ in range(2000):
        counts=np.bincount(rng.integers(0,n,n),minlength=n)
        w=counts/e**2
        LL=np.linalg.solve(X.T@(w[:,None]*X),X.T*w)
        hh,_=young_operator(age,z,e,counts=counts)
        hc,_=young_operator(age,z,e,fold,counts)
        xj=np.c_[X,z]; bj=np.linalg.solve(xj.T@(w[:,None]*xj),xj.T@(w*y))
        ix,bb=young_operator(age,z,e,counts=counts)[1][0]
        boot.append([*(bb@y[ix]),(LL@y)[1],(LL@(hh@y))[1],(LL@(hc@y))[1],bj[1]])
    draws=np.asarray(boot); np.savez_compressed(work/f'{name}-young-bootstrap.npz',draws=draws)
    labels=['young_intercept','young_z','original_age','same_sample_adjusted_age','crossfit_adjusted_age','joint_age']
    results={'n':n,'young_n':int(young.sum()),'young_regression_WLS':gls(np.c_[np.ones(young.sum()),z[young]],y[young],C[np.ix_(young,young)])[0], 'young_regression_OLS':np.linalg.lstsq(np.c_[np.ones(young.sum()),z[young]],y[young],rcond=None)[0], 'table_correction_linear_coefficients':coefficient,'table_correction_max_linear_remainder':float(np.max(abs(published_diff-np.c_[np.ones(n),z]@coefficient))), 'original_age':{'slope':b[1],'se':np.sqrt(v[1,1])},'same_sample_adjustment':{'slope':new[1],'naive_se':np.sqrt(v[1,1]),'propagated_se':np.sqrt(vin[1,1])},'crossfit_adjustment':{'slope':cross[1],'naive_se':np.sqrt(v[1,1]),'propagated_se':np.sqrt(vcf[1,1])}, 'joint_age_z':joint, 'bootstrap':{lab:{'sd':draws[:,i].std(ddof=1),'interval_95':np.quantile(draws[:,i],[.025,.975])} for i,lab in enumerate(labels)},'bootstrap_covariance_labels':labels,'bootstrap_covariance':np.cov(draws,rowvar=False),'same_minus_original':{'mean':np.mean(draws[:,3]-draws[:,2]),'sd':np.std(draws[:,3]-draws[:,2],ddof=1),'interval_95':np.quantile(draws[:,3]-draws[:,2],[.025,.975])}, 'fold_sizes':np.bincount(fold,minlength=5),'limits':'Fixed observed ages/young membership; original reported diagonal HR errors only. Shared correction covariance propagated conditional on that input; unreported inter-SN calibration and host-age covariance cannot be restored. Bootstrap resamples physical SNe and refits young trend; duplicate draws keep one fold.'}
    # Direct noise propagation independently verifies linear-operator variance.
    noise=rng.normal(size=(n,20000))*e[:,None]
    sim=L@H@noise
    results['linear_covariance_check']={'replicates':20000,'analytic_age_variance':vin[1,1],'monte_carlo_age_variance':np.var(sim[1],ddof=1),'relative_error':np.var(sim[1],ddof=1)/vin[1,1]-1}
    assert abs(results['linear_covariance_check']['relative_error'])<.05
    return results

def heldout(d,C,seed,work,policy):
    n=len(d); age=d.age.to_numpy(); zz=np.c_[np.ones(n),d.zHD,d.HOST_LOGMASS,d.c,d.x1]; Xs=[zz,np.c_[zz,age]]
    rng=np.random.default_rng(seed); sums={k:np.zeros((n,2)) for k in ['hr_corrected','hr_no_bias']}; score={k:[] for k in sums}
    Q=np.zeros((n,n))
    for repeat in range(20):
        order=rng.permutation(n); folds=np.empty(n,int); folds[order]=np.arange(n)%5
        row={k:np.zeros(2) for k in sums}
        for fold in range(5):
            tr=np.flatnonzero(folds!=fold); te=np.flatnonzero(folds==fold)
            ct=C[np.ix_(tr,tr)]; cte=C[np.ix_(te,tr)]; ctt=C[np.ix_(te,te)]
            ci=cho_factor(ct,lower=True); A=cho_solve(ci,cte.T).T; conditional=ctt-A@cte.T
            operators=[]
            for model,X in enumerate(Xs):
                _,v,linear,_=gls(X[tr],d.hr_corrected.to_numpy()[tr],ct)
                R=X[te]-A@X[tr]
                H=np.zeros((len(te),n)); H[:,te]=np.eye(len(te)); H[:,tr]=-A-R@linear
                operators.append(H)
                for ycol in sums:
                    y=d[ycol].to_numpy(); b,v,_,_=gls(X[tr],y[tr],ct)
                    R=X[te]-A@X[tr]; pred=X[te]@b+A@(y[tr]-X[tr]@b)
                    vv=conditional+R@v@R.T; ll=cholesky(vv,lower=True); residual=y[te]-pred; wr=solve_triangular(ll,residual,lower=True)
                    logp=-.5*len(te)*np.log(2*np.pi)-np.log(np.diag(ll)).sum()-.5*(wr@wr)
                    row[ycol][model]+=logp; sums[ycol][te,model]+=residual**2/20
            Q+=(operators[0].T@operators[0]-operators[1].T@operators[1])/(20*n)
        for k in sums: score[k].append(row[k])
    # Calibrate repeated-fold MSE gain with shared errors and repeated training
    # intact. This is an exactly specified Gaussian conditional-design null,
    # not an empirical age/dust null or a survey-selection calibration.
    noise=cholesky(C,lower=True)@rng.normal(size=(n,5000))
    null=np.sum(noise*(Q@noise),axis=0); threshold=np.quantile(null,.95)
    calibration={'null_replicates':5000,'null_MSE_gain_threshold_95':threshold,
                 'scope':'Full correlated Gaussian released covariance; fixed measured ages, nuisance design and folds. No age inference or selection uncertainty.',
                 'power_for_measured_age_slopes':{}}
    for b in [0.,-.01,-.03]:
        draw=noise+b*age[:,None]
        power=np.mean(np.sum(draw*(Q@draw),axis=0)>threshold)
        calibration['power_for_measured_age_slopes'][str(b)]={'power':power,'binomial_MC_se':np.sqrt(power*(1-power)/5000)}
    assert np.max(abs(Q@zz))<1e-9
    out={}
    for k in sums:
        loss=sums[k]; gain=loss[:,0]-loss[:,1]; bs=[]
        for _ in range(2000):bs.append(np.mean(gain[rng.integers(0,n,n)]))
        sc=np.asarray(score[k]);out[k]={'rmse_without_age':np.sqrt(loss[:,0].mean()),'rmse_with_age':np.sqrt(loss[:,1].mean()),'MSE_gain_positive_favours_age':gain.mean(),'paired_SN_bootstrap_MSE_gain_95':np.quantile(bs,[.025,.975]),'joint_conditional_logscore_gain_mean':np.mean(sc[:,1]-sc[:,0]),'joint_logscore_gain_repeat_range':[np.min(sc[:,1]-sc[:,0]),np.max(sc[:,1]-sc[:,0])],'age_better_MSE_repeats':'Not treated as independent replication; point scores are averaged per physical SN.','conditional_correlated_null':dict(calibration, p_gain_ge_observed=(1+np.sum(null>=gain.mean()))/(len(null)+1)),'limits':'Held-out folds use same survey and released train/test covariance. Nuisance parameter uncertainty included; fixed ages and fixed covariance. Paired SN bootstrap is descriptive and does not calibrate shared covariance, repeated training, or age inference uncertainty.'}
        assert abs(float(d[k].to_numpy()@Q@d[k].to_numpy())-gain.mean())<1e-12
        pd.DataFrame({'CID':d.CID,'loss_no_age':loss[:,0],'loss_age':loss[:,1]}).to_csv(work/f'{policy}-{k}-prediction.csv',index=False)
    return out

def eiv(d,C,ycol,work,policy):
    # Analytic integration of independent normal latent ages conditional on
    # fixed z/mass/c/x1. This intentionally exposes posterior-summary assumptions.
    age=d.age.to_numpy(); sa=d.age_err.to_numpy(); y=d[ycol].to_numpy()
    raw=np.c_[d.zHD,d.HOST_LOGMASS,d.c,d.x1]; means=raw.mean(0); scales=raw.std(0); Z=(raw-means)/scales; F=np.c_[np.ones(len(d)),Z]
    # [latent age mean coefficients(5), log tau, HR intercept/nuisance(5), age slope, log scatter]
    bounds=[(0,14)]+[(-8,8)]*4+[(-4,2.6)]+[(-3,3)]*5+[(-.2,.15),(-9,0)]
    def nll(t):
        m=F@t[:5];tau2=np.exp(2*t[5]);b=t[11];s2=np.exp(2*t[12]);av=tau2+sa**2;am=m+tau2/av*(age-m);pv=tau2*sa**2/av
        cc=C+np.diag(s2+b*b*pv); L=cholesky(cc,lower=True,check_finite=False); rr=solve_triangular(L,y-F@t[6:11]-b*am,lower=True,check_finite=False)
        return float(.5*np.sum(np.log(av)+(age-m)**2/av)+np.log(np.diag(L)).sum()+.5*rr@rr)
    starts=[]
    af=np.linalg.lstsq(F,age,rcond=None)[0]
    for tau,b in [(1,-.03),(2,-.01),(.3,-.06)]:
        st=np.r_[af,np.log(tau),np.mean(y)-b*np.mean(age),[0]*4,b,np.log(.04)];starts.append(st)
    fits=[minimize(nll,s,method='L-BFGS-B',bounds=bounds,options={'maxiter':1800,'ftol':1e-11,'gtol':1e-5}) for s in starts];best=min(fits,key=lambda f:f.fun)
    grid=np.arange(-.08,.040001,.004); keep=np.array([i for i in range(13) if i!=11]); curves=[]
    for b in grid:
        def ff(p): t=best.x.copy();t[keep]=p;t[11]=b;return nll(t)
        rr=minimize(ff,best.x[keep],method='L-BFGS-B',bounds=[bounds[i] for i in keep],options={'maxiter':1200,'ftol':1e-10,'gtol':1e-5})
        curves.append([b,rr.fun,bool(rr.success)])
    a=np.asarray(curves); base=min(best.fun,a[:,1].min()); allowed=grid[2*(a[:,1]-base)<=3.84]
    pd.DataFrame(curves,columns=['slope','nll','success']).to_csv(work/f'{policy}-{ycol}-eiv-profile.csv',index=False)
    return {'slope_mle':best.x[11],'profile_interval_delta_3p84_grid_004': [allowed.min(),allowed.max()] if len(allowed) else None,'tau_age_mle':np.exp(best.x[5]),'intrinsic_scatter_mle':np.exp(best.x[12]),'optimizer_success':bool(best.success),'multistart_nll': [f.fun for f in fits],'profile_successes':int(a[:,2].sum()),'profile_fits':len(grid),'nll':best.fun,'theta':best.x,'nuisance_predictor_means':means,'nuisance_predictor_scales':scales,'limits':'Normal classical-error approximation to posterior age summaries; unknown joint host likelihood, host/SALT measurement covariance and selection omitted. Nominal profile intervals uncalibrated, not Bayesian credible intervals.'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--archive',type=Path,default=DEFAULT_ARCHIVE);ap.add_argument('--work',type=Path,default=ROOT/'.work/age-reconciliation');ap.add_argument('--skip-eiv',action='store_true');args=ap.parse_args(); archive=args.archive; work=args.work; work.mkdir(parents=True,exist_ok=True)
    out=ROOT/'studies/host_ages/results/age_reconciliation';out.mkdir(parents=True,exist_ok=True)
    inputs=[archive/f'data/host_ages/chung2025/table{i}.dat' for i in [1,2]]
    ages=pd.concat([read_age(p,l) for p,l in zip(inputs,['G11','R19'])],ignore_index=True)
    ppdir=archive/'sources/repos/PantheonPlusSH0ES__DataRelease@7fc6805/Pantheon+_Data/4_DISTANCES_AND_COVAR'
    inputs += [ppdir/'Pantheon+SH0ES.dat',ppdir/'Pantheon+SH0ES_STAT+SYS.cov']
    pp=pd.read_csv(inputs[-2],sep=r'\s+',dtype={'CID':str});pp['pp_row']=np.arange(len(pp));a=np.loadtxt(inputs[-1]);n=int(a[0]);C=a[1:].reshape(n,n)
    source_covariance_asymmetry=np.max(abs(C-C.T));C=(C+C.T)/2
    joined=ages.merge(pp[pp.IDSURVEY==1],on='CID',how='left',validate='many_to_one',indicator=True)
    joined.to_csv(work/'all-source-crosswalk.csv',index=False)
    matched=joined[joined._merge=='both'].copy();matched['mu_ref']=mu(matched.zHD,matched.zHEL);matched['hr_corrected']=matched.MU_SH0ES-matched.mu_ref;matched['hr_no_bias']=matched.hr_corrected+matched.biasCor_m_b;matched['hr_tripp']=matched.mB+.148*matched.x1-3.112*matched.c+19.253-matched.mu_ref;matched['raw_magnitude_residual']=matched.mB+19.253-matched.mu_ref
    # Named signed components close exactly, including a recorded release remainder.
    matched['tripp_width']=.148*matched.x1;matched['tripp_colour']=-3.112*matched.c;matched['signed_released_bias']=-matched.biasCor_m_b
    matched['release_remainder']=matched.hr_corrected-(matched.raw_magnitude_residual+matched.tripp_width+matched.tripp_colour+matched.signed_released_bias)
    closure=matched.hr_corrected-(matched.raw_magnitude_residual+matched.tripp_width+matched.tripp_colour+matched.signed_released_bias+matched.release_remainder)
    assert np.max(abs(closure))<1e-12
    matched.to_csv(work/'signed-correction-ledger.csv',index=False)
    result={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'counts':{'source_age_rows':len(ages),'unique_source_SNe':ages.CID.nunique(),'duplicate_G11_R19_SNe':int(ages.CID.duplicated().sum()),'matched_age_rows':len(matched),'unmatched_age_rows':int((joined._merge!='both').sum())},'ledger':{'max_absolute_closure_mag':np.max(abs(closure)),'remainder_range_mag':[matched.release_remainder.min(),matched.release_remainder.max()],'bias_correction_scope':'Exported combined simulation bias; not a separately identifiable pure dust or mass step.','rows_file':'signed-correction-ledger.csv in --work directory'},'covariance':{'source_full_matrix_max_asymmetry':source_covariance_asymmetry,'treatment':'Symmetrized (C+C.T)/2 before row selection; the age-matched 196-row subset is exactly symmetric in the release.'},'samples':{}}
    for policy,keep in [('G11_first','first'),('R19_first','last')]:
        d=matched.drop_duplicates('CID',keep=keep);d=d[(d.zHD>.06)&(d.zHD<.42)].copy().reset_index(drop=True); assert d.CID.is_unique
        ii=d.pp_row.to_numpy(int);cc=C[np.ix_(ii,ii)];age=d.age.to_numpy();z=d.zHD.to_numpy();F=np.c_[np.ones(len(d)),z,d.HOST_LOGMASS,d.c,d.x1];X=np.c_[F[:,0],age,F[:,1:]]
        binid=np.digitize(z,[.06,.15,.25,.35,.42])-1;bincols=np.column_stack([binid==i for i in range(1,4)]).astype(float);Xbin=np.c_[np.ones(len(d)),age,bincols]
        fits={}
        for ycol in ['hr_C25','hr_original','hr_corrected','hr_no_bias','hr_tripp','raw_magnitude_residual']:
            metric=np.diag(d.hr_C25_err.to_numpy()**2) if ycol in ['hr_C25','hr_original'] else cc
            fits[ycol]={name:summary_gls(xx,d[ycol].to_numpy(),metric,ll) for name,xx,ll in [('age',X[:,:2],['intercept','age']),('age_z',X[:,:3],['intercept','age','z']),('age_z_mass_c_x1',X,['intercept','age','z','mass','c','x1']),('age_free_z_bins',Xbin,['intercept','age','bin1','bin2','bin3'])]}
        # Nuisance projection of measured age in the full covariance metric.
        aa,_,_,_=gls(F,age,cc);res=age-F@aa;one=np.ones((len(d),1));a0=age-one[:,0]*gls(one,age,cc)[0][0];W=cho_factor(cc)
        info=float(res@cho_solve(W,res));info0=float(a0@cho_solve(W,a0));within=[]
        for k in range(4):
            q=d[binid==k];within.append({'z_range':[[.06,.15],[.15,.25],[.25,.35],[.35,.42]][k],'n':len(q),'young_lt4':int((q.age<4).sum()),'old_ge4':int((q.age>=4).sum()),'age_minmax':[q.age.min(),q.age.max()],'age_q10_q90':np.quantile(q.age,[.1,.9])})
        # Correlated parametric noise calibrates paired slope contrasts conditional
        # on measured design and frozen covariance, not unknown correction errors.
        La=gls(X[:,:2],d.hr_corrected.to_numpy(),cc)[2][1];Lj=gls(X,d.hr_corrected.to_numpy(),cc)[2][1]
        items={'n':len(d),'fits':fits,'measured_age_information':{'fraction_left_after_z_mass_c_x1':info/info0,'conditional_slope_se':1/np.sqrt(info),'age_variance_minus_mean_error_squared':np.var(age,ddof=1)-np.mean(d.age_err**2),'bins':within},'heldout':heldout(d,cc,2026092701,work,policy),'corrected_joint_minus_marginal_conditional':{'difference':(Lj-La)@d.hr_corrected.to_numpy(),'se_full_covariance':np.sqrt((Lj-La)@cc@(Lj-La))}}
        der=(mu(d.zHD,d.zHEL,.301)-mu(d.zHD,d.zHEL,.299))/.002
        geometry={}
        for label,metric in [('plotting_error_diagonal',np.diag(d.MU_SH0ES_ERR_DIAG.to_numpy()**2)),('released_covariance_diagonal',np.diag(np.diag(cc))),('released_full_covariance',cc)]:
            oo=np.ones((len(d),1));aa=age-gls(oo,age,metric)[0][0];gg=der-gls(oo,der,metric)[0][0];ci=cho_factor(metric)
            frac=(aa@cho_solve(ci,gg))**2/(aa@cho_solve(ci,aa))/(gg@cho_solve(ci,gg))
            geometry[label]={'fraction_removed_by_one_local_Omega_m_derivative':frac,'retained_fraction':1-frac}
        items['cosmology_projection_weight_reconciliation']=geometry
        if not args.skip_eiv:
            # Prioritize common G11-first sample; the fixed-age analysis above runs both.
            if policy=='G11_first':items['joint_eiv']={y:eiv(d,cc,y,work,policy) for y in ['hr_corrected','hr_no_bias']}
        result['samples'][policy]=items;save(out/'summary.json',result)
    gp=archive/'data/host_ages/gupta2011';inputs += [gp/'table2.dat',gp/'ReadMe']
    old=ascii.read(gp/'table2.dat',format='cds',readme=str(gp/'ReadMe')).to_pandas();old['CID']=old.SNID.astype(str)
    gd=ages[ages.age_source=='G11'].merge(old,on='CID',validate='one_to_one');quality=(abs(gd.x1)<=3)&(gd.e_x1<=1)&(abs(gd.c)<=.3)&(gd.e_c<=.1)
    gd.to_csv(work/'G11-original-lightcurve-crosswalk.csv',index=False)
    result['original_quality_reconstruction']={'original_G11_rows':len(old),'C25_crossmatched':len(gd),'quality_rows':int(quality.sum()),'quality_young_lt4':int(((gd.age<4)&quality).sum()),'quality_young_le4':int(((gd.age<=4)&quality).sum()),'age_equal4_count':int((gd.age==4).sum()),'original_residual_table_rounding_max_mag':np.max(abs(gd.HR-gd.hr_original)),'status':'Original Gupta2011 columns; no substitution of modern Pantheon SALT parameters.'}
    result['young_host']={'C25_G11_all':young_test(gd,np.random.default_rng(2026092702),work,'C25'), 'Park_original_quality':young_test(gd[quality].reset_index(drop=True),np.random.default_rng(2026092703),work,'Park')}
    result['same_G11_objects_original_vs_updated_age']={}
    for label,q in [('all_199',gd),('quality_175',gd[quality])]:
        cc=np.diag(q.hr_C25_err.to_numpy()**2)
        result['same_G11_objects_original_vs_updated_age'][label]={}
        for age_col in ['Age','age']:
            xx=np.c_[np.ones(len(q)),q[age_col],q.z,q.M,q.c,q.x1]
            result['same_G11_objects_original_vs_updated_age'][label][age_col]={yc:summary_gls(xx,q[yc].to_numpy(),cc,['intercept','age','z','mass','c','x1']) for yc in ['hr_original','hr_C25']}
    # Fit the released algebra independently, rather than treating a named
    # remainder as proof of a physical decomposition.
    q=matched.drop_duplicates('CID');xx=np.c_[np.ones(len(q)),q.x1,q.c,(q.HOST_LOGMASS>10).astype(int)];target=q.m_b_corr-q.mB+q.biasCor_m_b
    coef=np.linalg.lstsq(xx,target,rcond=None)[0]
    result['ledger']['empirical_release_coefficients']={'labels':['constant','alpha','minus_beta','mass_step'],'values':coef,'max_absolute_remainder_mag':np.max(abs(target-xx@coef))}
    # Independent GLS implementation checks using explicit inverse on a subset.
    d=matched.drop_duplicates('CID').iloc[:20];ii=d.pp_row.to_numpy(int);co=C[np.ix_(ii,ii)];xx=np.c_[np.ones(len(d)),d.age]; yy=d.hr_corrected.to_numpy();direct=np.linalg.solve(xx.T@np.linalg.solve(co,xx),xx.T@np.linalg.solve(co,yy));b=gls(xx,yy,co)[0]; assert np.max(abs(direct-b))<1e-10
    # Test the integrated multivariate conditional-age model against a single
    # block Gaussian, with nonconstant age means and nonconstant HR nuisance.
    rng=np.random.default_rng(2026092704);ncheck=17
    q=rng.normal(size=(ncheck,ncheck));co=q@q.T*.0002+np.eye(ncheck)*.01
    age_sigma=rng.uniform(.3,2,ncheck);mean_age=rng.uniform(3,7,ncheck);tau=1.3;slope=-.03;scatter=.025
    age=mean_age+rng.normal(size=ncheck);extra=rng.normal(0,.03,ncheck);yy=rng.normal(0,.1,ncheck)
    av=tau*tau+age_sigma**2;am=mean_age+tau*tau/av*(age-mean_age);pv=tau*tau*age_sigma**2/av
    def gnll(r,co):
        ll=cholesky(co,lower=True);vv=solve_triangular(ll,r,lower=True)
        return np.log(np.diag(ll)).sum()+.5*vv@vv
    conditional=.5*np.sum(np.log(av)+(age-mean_age)**2/av)+gnll(yy-extra-slope*am,co+np.diag(scatter**2+slope*slope*pv))
    block=np.block([[np.diag(av),np.eye(ncheck)*slope*tau*tau],[np.eye(ncheck)*slope*tau*tau,co+np.eye(ncheck)*(slope*slope*tau*tau+scatter**2)]])
    joint=gnll(np.r_[age-mean_age,yy-extra-slope*mean_age],block)
    assert abs(conditional-joint)<1e-10
    result['validation']={'gls_independent_max_abs_parameter_difference':np.max(abs(direct-b)),'eiv_block_vs_conditional_nll_difference':abs(conditional-joint),'closure_pass':True,'physical_SN_unique_primary_samples':True}
    save(out/'summary.json',result)
    # The live redirected execution log is still being written after this
    # manifest; only completed numerical products belong in output identities.
    generated=[p for p in work.iterdir() if p.is_file() and p.suffix in {'.csv','.npz'}]
    manifest={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'python':platform.python_version(),'versions':{n:__import__(n).__version__ for n in ['numpy','scipy','pandas','astropy']},'archive_snapshot':archive.name,'inputs':{str(p.relative_to(archive)):digest(p) for p in inputs},'code_sha256':digest(Path(__file__)),'design_sha256':digest(HERE/'design.json'),'summary_sha256':digest(out/'summary.json'),'outputs':{str(p.relative_to(work)):digest(p) for p in generated},'command':'.venv/bin/python studies/host_ages/code/age_reconciliation/run.py --archive <restored-snapshot-directory>','scope':'New conditional observational analyses; no new LINMIX chain, survey refit, causal identification or cosmological measurement.'}
    save(out/'manifest.json',manifest);print(json.dumps({'samples':{k:v['n'] for k,v in result['samples'].items()},'quality':result['original_quality_reconstruction'],'output':str(out)},indent=2,default=serial))

if __name__=='__main__':main()
