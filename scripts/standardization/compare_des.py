#!/usr/bin/env python3
"""STD-01: corrected-distance release comparison, not photometry reconstruction."""
from pathlib import Path
import hashlib, json, platform, subprocess
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve, eigvalsh
from scipy.optimize import minimize_scalar, brentq
from scipy.integrate import quad
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'runs/standardization/des_comparison'
DER = ROOT/'data/derived/standardization'
OUT.mkdir(parents=True, exist_ok=True)
DER.mkdir(parents=True, exist_ok=True)
OLD = ROOT/'sources/repos/des-science__DES-SN5YR@1.3/4_DISTANCES_COVMAT'
NEW = ROOT/'sources/repos/des-science__DES-SN5YR/4_DISTANCES_COVMAT'
PP = ROOT/'sources/repos/PantheonPlusSH0ES__DataRelease@7fc6805/Pantheon+_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES.dat'
inputs = []

def sn_table(path):
    inputs.append(path)
    if path.read_text().lstrip().startswith('CID,'):
        d = pd.read_csv(path, dtype={'CID':str})
    else:
        d = pd.read_csv(path, sep=r'\s+', comment='#', dtype={'CID':str})
        if 'VARNAMES:' in d:
            d = d.drop(columns='VARNAMES:')
    d['CID'] = d['CID'].astype(str)
    return d

def key(d):
    return pd.MultiIndex.from_frame(d[['CID','IDSURVEY']])

def old_cov(name, d):
    path=OLD/name; inputs.append(path)
    flat=np.loadtxt(path)
    assert int(flat[0])==len(d)
    c=flat[1:].reshape(len(d),len(d))
    return c+np.diag(d.MUERR_FINAL.to_numpy()**2)

def new_precision(name):
    path=NEW/name; inputs.append(path)
    d=np.load(path)
    n=int(d['nsn'][0]); p=np.zeros((n,n))
    p[np.triu_indices(n)] = d['cov']
    p+=np.triu(p,1).T
    return p

def inverse(c):
    return cho_solve(cho_factor(c,lower=True),np.eye(len(c)))

def diagnostics(c):
    cf=cho_factor(c,lower=True)
    return {'n':len(c),'symmetric_max_abs':float(np.max(abs(c-c.T))),
            'cholesky':True,'min_eigenvalue':float(eigvalsh(c,subset_by_index=[0,0])[0]),
            'max_eigenvalue':float(eigvalsh(c,subset_by_index=[len(c)-1,len(c)-1])[0]),
            'median_sigma':float(np.median(np.sqrt(np.diag(c)))),
            'logdet':float(2*np.log(np.diag(cf[0])).sum())}

gx,gw=np.polynomial.legendre.leggauss(64)
def mu_lcdm(om,z,zhel):
    x=z[:,None]*(gx+1)/2
    integral=z/2*np.sum(gw/np.sqrt(om*(1+x)**3+1-om),axis=1)
    return 5*np.log10(299792.458/70*(1+zhel)*integral)+25

def fit(d,p,label):
    z=d.zHD.to_numpy(); zh=d.zHEL.to_numpy(); obs=d.MU.to_numpy()
    one=np.ones(len(d)); po=p@one; den=one@po
    def objective(om,ret=False):
        r=obs-mu_lcdm(om,z,zh)
        offset=(po@r)/den
        r=r-offset
        chi=float(r@p@r)
        return (chi,float(offset)) if ret else chi
    result=minimize_scalar(objective,bounds=(0.01,0.99),method='bounded',options={'xatol':1e-10})
    om=result.x; ch,offset=objective(om,True)
    lo=brentq(lambda x:objective(x)-ch-1,.001,om)
    hi=brentq(lambda x:objective(x)-ch-1,om,1.5)
    # Independent quadrature spot check, including largest redshift.
    ix=np.unique([0,len(d)//2,int(np.argmax(z))])
    q=np.array([5*np.log10(299792.458/70*(1+zh[i])*quad(lambda x:1/np.sqrt(om*(1+x)**3+1-om),0,z[i],epsabs=1e-12)[0])+25 for i in ix])
    return {'label':label,'n':len(d),'Omega_m':float(om),'profile_delta_chi2_1_low':float(lo),
            'profile_delta_chi2_1_high':float(hi),'chi2_profile':ch,'chi2_with_M_marginalization_term':float(ch+np.log(den/(2*np.pi))),'dof':len(d)-2,
            'offset_mag':offset,'quad_check_max_mag':float(np.max(abs(q-mu_lcdm(om,z[ix],zh[ix])))),
            'assumptions':'Flat LCDM; free common magnitude intercept; no Cepheid/CMB/BAO input; published fixed distance covariance.'}

def summary(d,group):
    records=[]
    for g,s in d.groupby(group,observed=True):
        row={'group':str(g),'n':len(s)}
        for col in ['delta_mu','delta_amp','delta_stretch','delta_colour','delta_host','delta_bias','delta_intercept','delta_closure']:
            row[col+'_mean']=float(s[col].mean())
            row[col+'_median']=float(s[col].median())
        row['delta_mu_sd']=float(s.delta_mu.std())
        records.append(row)
    return records

def main():
    a=sn_table(OLD/'DES-SN5YR_HD.csv'); b=sn_table(NEW/'DES-Dovekie_HD.csv')
    am=sn_table(OLD/'DES-SN5YR_HD+MetaData.csv'); bm=sn_table(NEW/'DES-Dovekie_Metadata.csv')
    pp=sn_table(PP)
    audit={}
    pp_x=np.column_stack([np.ones(len(pp)),pp.x1,-pp.c,np.where(pp.HOST_LOGMASS>10,1,-1)/2])
    pp_y=pp.m_b_corr+pp.biasCor_m_b-pp.mB
    pp_coef=np.linalg.lstsq(pp_x,pp_y,rcond=None)[0]
    pp_err=pp_y-pp_x@pp_coef
    audit['Pantheon_algebraic_standardization_closure']={'constant_alpha_beta_residual_mass_step':pp_coef.tolist(),
        'rms_mag':float(np.sqrt(np.mean(pp_err**2))),'max_abs_mag':float(np.max(abs(pp_err))),
        'interpretation':'Algebraic diagnostic only; mass dependence remains inside biasCor_m_b even when residual explicit mass step is near zero.'}
    for label,d in [('original_HD',a),('Dovekie_HD',b),('original_metadata',am),('Dovekie_metadata',bm),('Pantheon_W26',pp)]:
        audit[label]={'rows':len(d),'unique_CID':int(d.CID.nunique()),'unique_CID_survey':len(key(d).unique()),'surveys':{str(k):int(v) for k,v in d.IDSURVEY.value_counts().items()}}
        if label!='Pantheon_W26': assert not key(d).duplicated().any(),label
    am=am.set_index(['CID','IDSURVEY']).loc[key(a)].reset_index()
    bm_orig=key(bm)
    bm=bm.set_index(['CID','IDSURVEY']).loc[key(b)].reset_index()
    audit['Dovekie_metadata_order_equal']=bool(bm_orig.equals(key(b)))
    audit['Dovekie_HD_vs_metadata_mu_max']=float(np.max(abs(b.MU-bm.MU)))
    audit['original_HD_vs_metadata_mu_max']=float(np.max(abs(a.MU-am.MU)))
    common=key(a).intersection(key(b),sort=False)
    ai=key(a).get_indexer(common); bi=key(b).get_indexer(common)
    aa=a.iloc[ai].reset_index(drop=True); bb=b.iloc[bi].reset_index(drop=True)
    ama=am.iloc[ai].reset_index(drop=True); bmb=bm.iloc[bi].reset_index(drop=True)
    audit['matched']=len(common)
    audit['original_only']=[list(x) for x in key(a).difference(key(b))]
    audit['Dovekie_only']=[list(x) for x in key(b).difference(key(a))]
    audit['matched_max_delta_zHD']=float(np.max(abs(aa.zHD-bb.zHD)))
    audit['matched_max_delta_zHEL']=float(np.max(abs(aa.zHEL-bb.zHEL)))
    audit['Pantheon_overlap']={}
    for name,d in [('original',a),('Dovekie',b)]:
        shared=set(d.CID)&set(pp.CID)
        audit['Pantheon_overlap'][name]={'literal_unique_CID':len(shared),'literal_CID_survey':len(key(d).intersection(key(pp))),
                                      'by_DES_survey':{str(k):int(v) for k,v in d.loc[d.CID.isin(shared),'IDSURVEY'].value_counts().items()}}
    oldstat=old_cov('STATONLY.txt.gz',a); oldtot=old_cov('STAT+SYS.txt.gz',a)
    newp=new_precision('STAT+SYS.npz'); newps=new_precision('STATONLY.npz')
    newtot=inverse(newp); newstat=inverse(newps)
    mats={}
    for label,c in [('original_stat',oldstat),('original_total',oldtot),('Dovekie_stat',newstat),('Dovekie_total',newtot)]: mats[label]=diagnostics(c)
    mats['Dovekie_stat_sigma_vs_HD_max_abs']=float(np.max(abs(np.sqrt(np.diag(newstat))-b.MUERR)))
    mats['Dovekie_total_inverse_residual_max']=float(np.max(abs(newp@newtot-np.eye(len(b)))))
    fits=[fit(a,inverse(oldtot),'native_original_total'),fit(b,newp,'native_Dovekie_total'),fit(a,inverse(oldstat),'native_original_stat'),fit(b,newps,'native_Dovekie_stat')]
    ca=oldtot[np.ix_(ai,ai)]; cb=newtot[np.ix_(bi,bi)]
    pa=inverse(ca); pb=inverse(cb)
    for name,d in [('original',aa),('Dovekie',bb)]:
        for cname,p in [('original',pa),('Dovekie',pb)]: fits.append(fit(d,p,f'matched_{name}_distances_{cname}_covariance'))
    # Documentation-based decomposition (rounded published coefficients).
    out=aa[['CID','IDSURVEY','zHD','zHEL']].copy()
    out['delta_mu']=bb.MU-aa.MU
    out['delta_amp']=-2.5*np.log10(bmb.x0/ama.x0)
    out['delta_stretch']=.169*bmb.x1-.16087*ama.x1
    out['delta_colour']=-3.14*bmb.c+3.11780*ama.c
    out['delta_host']=.033/2*np.where(bmb.HOST_LOGMASS>10,1,-1)-.03754/2*np.where(ama.HOST_LOGMASS>10,1,-1)
    out['delta_bias']=-bmb.biasCor_mu+ama.biasCor_mu
    out['delta_intercept']=29.96210-29.95821
    components=['delta_amp','delta_stretch','delta_colour','delta_host','delta_bias','delta_intercept']
    out['delta_closure']=out.delta_mu-out[components].sum(axis=1)
    out['survey_group']=np.where(out.IDSURVEY==10,'DES',np.where(out.IDSURVEY==150,'Foundation','Other low-z'))
    out['redshift_bin']=pd.cut(out.zHD,[0,.1,.3,.5,.7,1,2],right=False)
    out.to_csv(DER/'des_matched_comparison.csv',index=False)
    # Linear algebra closure diagnostic: recover documented coefficients from release alone,
    # not infer physical causes. Account for smooth host step at tau=.001.
    closure={}
    from scipy.special import expit
    for label,hd,meta in [('original',a,am),('Dovekie',b,bm)]:
        step=expit((meta.HOST_LOGMASS.to_numpy()-10)/.001)-.5
        X=np.column_stack([np.ones(len(meta)),meta.x1,-meta.c,step])
        y=hd.MU.to_numpy()+2.5*np.log10(meta.x0.to_numpy())+meta.biasCor_mu.to_numpy()
        coef=np.linalg.lstsq(X,y,rcond=None)[0]; err=y-X@coef
        closure[label]={'inferred_neg_M0_alpha_beta_gamma':coef.tolist(),'rms_mag':float(np.sqrt(np.mean(err**2))),
                        'max_abs_mag':float(np.max(abs(err))), 'host_step_width_logmass':.001}
        away=abs(meta.HOST_LOGMASS.to_numpy()-10)>.02
        coef_away=np.linalg.lstsq(X[away],y[away],rcond=None)[0]
        error_away=y-X@coef_away
        closure[label]['away_from_rounded_mass_boundary']={'n':int(away.sum()),'coefficients':coef_away.tolist(),
            'rms_mag':float(np.sqrt(np.mean(error_away[away]**2))),'max_abs_mag':float(np.max(abs(error_away[away])))}
        closure[label]['objects_with_abs_closure_over_1_mmag']=[{'CID':str(meta.CID.iloc[i]),'mass':float(meta.HOST_LOGMASS.iloc[i]),'error_mag':float(error_away[i])} for i in np.flatnonzero(abs(error_away)>.001)]
    res={'audit':audit,'matrix_diagnostics':mats,'fits':fits,'decomposition_by_survey':summary(out,'survey_group'),
         'decomposition_by_redshift':summary(out,'redshift_bin'),'documented_decomposition_closure_rms_mag':float(np.sqrt(np.mean(out.delta_closure**2))),
         'documented_decomposition_closure_max_abs_mag':float(np.max(abs(out.delta_closure))), 'linear_algebra_closure_diagnostic':closure,
         'difference_uncertainty_warning':'Shared observations/calibration; cross-release covariance unavailable. No independent-error significance for release differences.'}
    (OUT/'results.json').write_text(json.dumps(res,indent=2)+'\n')
    pd.DataFrame(fits).to_csv(OUT/'flat_lcdm_fits.csv',index=False)
    fig,axs=plt.subplots(1,2,figsize=(12,4))
    for label,g in out.groupby('survey_group'):
        axs[0].scatter(g.zHD,g.delta_mu,s=6,alpha=.35,label=label)
    axs[0].set(xlabel='Hubble-diagram redshift',ylabel='Dovekie − original distance (mag)',ylim=(-.28,.28))
    axs[0].legend(fontsize=8); axs[0].axhline(0,c='k',lw=.5)
    rows=res['decomposition_by_survey']; x=np.arange(len(rows))
    bottom=np.zeros(len(x))
    for col in components:
        vals=np.array([r[col+'_mean'] for r in rows])
        axs[1].plot(x,vals,marker='o',label=col.replace('delta_',''))
    axs[1].set(xticks=x,xticklabels=[r['group'] for r in rows],ylabel='Mean change of term (mag)')
    axs[1].legend(fontsize=8,ncol=2); axs[1].axhline(0,c='k',lw=.5)
    fig.suptitle('Same supernovae, jointly changed distance pipeline (no cross-release significance)')
    fig.tight_layout(); fig.savefig(OUT/'matched_distance_changes.png',dpi=170); plt.close(fig)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    manifest={'experiment':'STD-01','utc':datetime.now(timezone.utc).isoformat(),'code_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
              'script_sha256':sha(Path(__file__)),'plan_sha256':sha(ROOT/'docs/experiments/standardization-plan.md'),
              'inputs':{str(p.relative_to(ROOT)):sha(p) for p in sorted(set(inputs))},
              'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'pandas':pd.__version__},
              'configuration':{'Omega_m_bounds':[.01,.99],'profile_delta_chi2':1,'magnitude_intercept':'analytic free','quadrature':'64 point Gauss-Legendre checked against scipy quad',
                               'matrix_subset_rule':'Invert precision to covariance; subset covariance; invert subset','seed':None},
              'outputs':[str(p.relative_to(ROOT)) for p in sorted(OUT.iterdir()) if p.name!='manifest.json']+['data/derived/standardization/des_matched_comparison.csv']}
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'audit':audit,'fits':fits,'closure':closure},indent=2))

if __name__=='__main__': main()
