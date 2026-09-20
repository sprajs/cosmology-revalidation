#!/usr/bin/env python3
"""Public-table age-HR reconstruction; never raw-photometry reconstruction."""
from pathlib import Path
import argparse, hashlib, json, platform, subprocess, sys, datetime
import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.linalg import cho_factor, cho_solve
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/age_signal'; DER=ROOT/'data/derived/age_signal'
PP=ROOT/'sources/repos/PantheonPlusSH0ES__DataRelease@7fc6805/Pantheon+_Data/4_DISTANCES_AND_COVAR'
AGES=[ROOT/f'data/host_ages/chung2025/table{i}.dat' for i in [1,2]]
CUTS=[.20,.25,.30,.35,.42]
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def save_manifest(label, inputs, config, outputs):
    manifest={'time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'purpose':label,'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'python':sys.version,'platform':platform.platform(),'versions':{k:__import__(k).__version__ for k in ['numpy','scipy','pandas','matplotlib']},'inputs':{str(p.relative_to(ROOT)):digest(p) for p in inputs},'code_sha256':digest(Path(__file__)),'plan_sha256':digest(ROOT/'docs/experiments/age_signal-plan.md'),'config':config,'outputs':{str(p.relative_to(ROOT)):digest(p) for p in outputs}}
    (OUT/f'{label}-manifest.json').write_text(json.dumps(manifest,indent=2))
def read_age(p, label):
    rows=[]
    for line in p.read_text().splitlines()[1:]:
        c=[v.strip().replace('\\','').strip() for v in line.split('&')]
        if len(c)!=5: continue
        hr=c[3].split('(')
        rows.append([str(int(c[0])),label,float(c[1]),float(c[2]),float(hr[0]),float(hr[-1].strip(') ')),float(c[4])])
    return pd.DataFrame(rows,columns=['CID','age_source','age','age_err','hr_C25','hr_original','hr_C25_err'])
def mu(z, zhel=None, om=.3):
    z=np.asarray(z); zhel=z if zhel is None else np.asarray(zhel)
    d=np.array([quad(lambda t: 1/np.sqrt(om*(1+t)**3+1-om),0,zz,epsabs=1e-11)[0] for zz in z])
    return 5*np.log10((1+zhel)*299792.458/70*d)+25

def prepare():
    pp=pd.read_csv(PP/'Pantheon+SH0ES.dat',sep=r'\s+',dtype={'CID':str}); pp['pp_row']=np.arange(len(pp))
    a,b=[read_age(p,l) for p,l in zip(AGES,['G11','R19'])]; ages=pd.concat([a,b],ignore_index=True)
    ages.to_csv(DER/'all_age_rows.csv',index=False)
    overlap=a.merge(b,on='CID',suffixes=('_G11','_R19')); overlap.to_csv(DER/'age_overlaps.csv',index=False)
    join=ages.merge(pp[pp.IDSURVEY==1],on='CID',how='left',indicator=True,validate='many_to_one')
    join.to_csv(DER/'crosswalk_all.csv',index=False)
    full=join[join._merge=='both'].drop(columns='_merge').copy()
    full['mu_model']=mu(full.zHD,full.zHEL)
    full['hr_corrected']=full.MU_SH0ES-full.mu_model
    full['hr_no_bias']=full.hr_corrected+full.biasCor_m_b
    # A pure Tripp residual need not equal merely reversing biasCor_m_b:
    full['hr_tripp']=full.mB+.148*full.x1-3.112*full.c+19.253-full.mu_model
    full['error_diag']=full.MU_SH0ES_ERR_DIAG
    full['error_no_covadd']=np.sqrt(np.maximum(0,full.error_diag**2-full.biasCor_m_b_COVADD))
    count={}
    for pol in ['all_rows','G11_first','R19_first']:
        d=full if pol=='all_rows' else full.drop_duplicates('CID',keep=('first' if pol=='G11_first' else 'last'))
        d.to_csv(DER/f'matched_{pol}.csv',index=False)
        count[pol]={'matched_rows':len(d),'unique_SN':d.CID.nunique(),'cut_counts_zHD':[int(((d.zHD>.06)&(d.zHD<q)).sum()) for q in CUTS], 'cut_counts_zHEL':[int(((d.zHEL>.06)&(d.zHEL<q)).sum()) for q in CUTS]}
    report={'source_rows':{'G11':len(a),'R19':len(b)},'shared_SN':len(overlap),'sdss_distance_rows':len(pp[pp.IDSURVEY==1]),'matched_age_rows':len(full),'unmatched_age_rows':int((join._merge!='both').sum()),'counts':count,'covadd_min':float(full.biasCor_m_b_COVADD.min()),'covadd_max':float(full.biasCor_m_b_COVADD.max()),'no_covadd_nonpositive':int((full.error_no_covadd<=0).sum()),'tripp_minus_no_bias_range':[float((full.hr_tripp-full.hr_no_bias).min()),float((full.hr_tripp-full.hr_no_bias).max())]}
    (OUT/'sample-audit.json').write_text(json.dumps(report,indent=2)); print(json.dumps(report,indent=2))
    save_manifest('A1-sample-audit',AGES+[PP/'Pantheon+SH0ES.dat'],{'cuts':CUTS,'zmin':.06,'survey':1,'overlap_policies':list(count)},list(DER.glob('*.csv'))+[OUT/'sample-audit.json'])

def wls(x,y,e,X=None,cov=None):
    if X is None: X=np.column_stack([np.ones(len(x)),x])
    if cov is None:
        wx=X/e[:,None]; wy=y/e
        c=np.linalg.inv(wx.T@wx); b=c@wx.T@wy; chi=float(np.sum(((y-X@b)/e)**2))
    else:
        cf=cho_factor(cov); c=np.linalg.inv(X.T@cho_solve(cf,X)); b=c@X.T@cho_solve(cf,y); chi=float((y-X@b)@cho_solve(cf,y-X@b))
    return b,c,chi

def diagnostics():
    rows=[]; bootstrap={}; rng=np.random.default_rng(20260920)
    for pol in ['all_rows','G11_first','R19_first']:
        full=pd.read_csv(DER/f'matched_{pol}.csv'); full=full[(full.zHD>.06)&(full.zHD<.42)]
        for correction in ['hr_corrected','hr_no_bias','hr_tripp']:
            for cut in CUTS:
                d=full[full.zHD<cut]; x=d.age.to_numpy(); y=d[correction].to_numpy(); e=d.error_diag.to_numpy()
                for formula in ['age','age_z','age_z_mass_c']:
                    X={'age':np.column_stack([np.ones(len(d)),x]),'age_z':np.column_stack([np.ones(len(d)),x,d.zHD]),'age_z_mass_c':np.column_stack([np.ones(len(d)),x,d.zHD,d.HOST_LOGMASS,d.c])}[formula]
                    b,c,chi=wls(x,y,e,X)
                    rows.append({'policy':pol,'correction':correction,'zcut':cut,'n':len(d),'formula':formula,'slope':b[1],'slope_se':np.sqrt(c[1,1]),'chi2':chi,'dof':len(d)-X.shape[1],'age_z_r':np.corrcoef(x,d.zHD)[0,1]})
            if pol=='all_rows': continue # primary paired inference uses unique-SN samples
            slopes=[]
            for i in range(2000):
                d=full.iloc[rng.integers(0,len(full),len(full))]; bs=[]
                for cut in CUTS:
                    q=d[d.zHD<cut]; bs.append(wls(q.age.to_numpy(),q[correction].to_numpy(),q.error_diag.to_numpy())[0][1])
                slopes.append(bs)
            slopes=np.array(slopes); contrast=slopes[:,0]-slopes[:,-1]
            bootstrap[f'{pol}:{correction}']={'covariance':np.cov(slopes,rowvar=False).tolist(),'narrow_minus_full_mean':float(contrast.mean()),'narrow_minus_full_sd':float(contrast.std()),'contrast_95pct':np.quantile(contrast,[.025,.975]).tolist()}
    pd.DataFrame(rows).to_csv(OUT/'wls-factorial.csv',index=False)
    (OUT/'cut-bootstrap.json').write_text(json.dumps(bootstrap,indent=2))
    save_manifest('A2-A3-wls',[DER/f'matched_{p}.csv' for p in ['all_rows','G11_first','R19_first']],{'seed':20260920,'bootstrap':2000,'note':'WLS ignores age uncertainty; diagnostic only'},[OUT/'wls-factorial.csv',OUT/'cut-bootstrap.json'])

def rhat_four(values):
    v=np.asarray(values).reshape(4,-1); n2=v.shape[1]//2; v=np.concatenate([v[:,:n2],v[:,-n2:]],axis=0); n=v.shape[1]; W=np.mean(np.var(v,axis=1,ddof=1)); B=n*np.var(v.mean(axis=1),ddof=1); return float(np.sqrt(((n-1)/n*W+B/n)/W))

def linmix_run():
    vendor=ROOT/'sources/external/age_signal/linmix-933dbb1359dcb5404cf881bfdd5cf433b0195152'; sys.path.insert(0,str(vendor)); import linmix
    tasks=[]
    # No overlap choice is declared an exact author table. R19-first chosen because R19 uses higher-quality nearby photometry.
    for pol in ['G11_first','R19_first']:
        d=pd.read_csv(DER/f'matched_{pol}.csv'); d=d[(d.zHD>.06)&(d.zHD<.42)]
        variants=[('hr_corrected','error_diag',.42),('hr_no_bias','error_diag',.42),('hr_tripp','error_diag',.42),('hr_tripp','error_no_covadd',.42),('hr_tripp','error_no_covadd',.2)]
        if pol!='all_rows': variants += [('hr_tripp','error_no_covadd',q) for q in [.25,.30,.35]]
        for y,e,cut in variants: tasks.append((f'{pol}_{y}_{e}_{cut}',d[d.zHD<cut],y,e))
    for path,label in zip(AGES,['G11','R19']):
        d=read_age(path,label)
        for y in (['hr_C25','hr_original'] if label=='G11' else ['hr_C25']): tasks.append((f'C25_{label}_{y}',d,y,'hr_C25_err'))
    results=[]
    for i,(name,d,y,e) in enumerate(tasks):
        out=OUT/'linmix'/f'{name}.json'; out.parent.mkdir(exist_ok=True)
        if out.exists(): results.append(json.loads(out.read_text())); continue
        print('START',name,len(d),flush=True)
        lm=linmix.LinMix(d.age.to_numpy(),d[y].to_numpy(),xsig=d.age_err.to_numpy(),ysig=d[e].to_numpy(),K=3,nchains=4,parallelize=True,seed=20260920+i)
        lm.run_mcmc(miniter=20000,maxiter=100000,silent=False)
        chain=lm.chain; np.savez_compressed(out.with_suffix('.npz'),chain=chain)
        beta=chain['beta']; item={'name':name,'n':len(d),'mean':float(beta.mean()),'sd':float(beta.std()),'median':float(np.median(beta)),'q16_84':np.quantile(beta,[.16,.84]).tolist(),'q025_975':np.quantile(beta,[.025,.975]).tolist(),'p_negative':float(np.mean(beta<0)),'scatter_mean':float(np.sqrt(chain['sigsqr']).mean()),'retained_draws':len(chain),'seed':20260920+i,'rhat_retained':{k:rhat_four(chain[k]) for k in ['alpha','beta','sigsqr']}}
        out.write_text(json.dumps(item,indent=2)); results.append(item); print(json.dumps(item),flush=True)
    (OUT/'linmix-summary.json').write_text(json.dumps(results,indent=2))
    save_manifest('A2-A3-A5-linmix',AGES+[DER/f'matched_{p}.csv' for p in ['all_rows','G11_first','R19_first']]+[vendor/'linmix/linmix.py'],{'K':3,'nchains':4,'miniter':20000,'maxiter':100000,'seed_base':20260920,'intrinsic_variance':'free nonnegative'},[OUT/'linmix-summary.json']+list((OUT/'linmix').glob('*.json')))
if __name__=='__main__':
    OUT.mkdir(exist_ok=True,parents=True); DER.mkdir(exist_ok=True,parents=True)
    cmd=argparse.ArgumentParser();cmd.add_argument('action',choices=['prepare','diagnostics','linmix']);args=cmd.parse_args()
    {'prepare':prepare,'diagnostics':diagnostics,'linmix':linmix_run}[args.action]()
