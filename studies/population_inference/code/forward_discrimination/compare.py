"""Finite-ensemble conditional distribution checks, not a population likelihood.

All inferential choices are in the hashed preregistration. The forecast is the
weighted empirical distribution itself, so energy scores use its V statistic.
"""
from pathlib import Path
import argparse, datetime, hashlib, json, os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'phase2/forward_discrimination'
FWD = ROOT / 'phase2/hierarchy/forward'
MODELS = ['P21','BS21','G10','P21_dmplus010','P21_dmminus010','P21_dmz020','P21_rho000','P21_rho090','P21_noisetrue120']
KNOTS = np.array([.025,.2,.35,.5,.65,.8,1.2])
INPUTS = set()

def source(path):
    path = Path(path); INPUTS.add(path); return path

def mu(z, zh, q=0):
    r = np.log1p(z) if q == 0 else -np.expm1(-q*np.log1p(z))/q
    return 5*np.log10((299792.458/70)*(1+zh)*r)+25

def basis(z):
    return np.column_stack([np.interp(z,KNOTS,np.eye(len(KNOTS))[j]) for j in range(len(KNOTS))])

def quality(f, strict=False):
    ans = ((f.x1.abs()<3)&(f.c.abs()<.3)&(f.x1ERR<1)&(f.PKMJDERR<2)&(f.cERR<1.5)&(f.FITPROB>.001)&(f.zHD>.025)&(f.zHD<1.2))
    if strict: ans &= (f.c.abs()<.2)&(f.x1ERR<.7)&(f.PKMJDERR<1)&(f.FITPROB>.01)
    return ans

def real_data(arm, all_des=False, strict=False, q=0):
    root = ROOT/'phase2/hierarchy/data'/arm
    z = np.load(source(root/'data.npz')); f = pd.read_csv(source(root/'rows.csv'), dtype={'CID':str})
    prefix = {'conditioned-multistart-best':'conditioned_multistart_best_refit','recovered-mask':'recovered_refit'}[arm]
    ob = pd.read_csv(source(ROOT/f'phase2/official/results/{prefix}_observables.csv'),dtype={'CID':str})
    ix=ob.set_index(['IDSURVEY','CID']).loc[list(zip(f.IDSURVEY,f.CID))].reset_index()
    use = (z['survey']==10)&quality(ix,strict).to_numpy()
    if not all_des: use &= z['pIa']>.999
    f['z']=z['z'];f['mass']=z['host_mass'];f['masserr']=z['host_mass_error'];f['host']=(f.mass>=10).astype(int)
    f['u']=z['y'][:,0]-mu(z['z'],z['zhel'],q);f['x']=z['y'][:,1];f['c']=z['y'][:,2]
    f['sigma_m']=np.sqrt(z['cov'][:,0,0]);f['sigma_x']=np.sqrt(z['cov'][:,1,1]);f['sigma_c']=np.sqrt(z['cov'][:,2,2]);f['fitprob']=ix.FITPROB.to_numpy()
    assert (np.linalg.eigvalsh(z['cov'][use]).min(axis=1)>0).all()
    f['same_basic_quality']=quality(ix,strict).to_numpy()
    return f.loc[use].reset_index(drop=True), {'n_des_before_cuts':int((z['survey']==10).sum()),'n_pass':int(use.sum()),'n_test_before_support':int(((f.fold==0)&use).sum())}

def sim_data(model, strict=False):
    f=pd.read_csv(source(FWD/(model+'-fitted.csv.gz')))
    cv=np.load(source(FWD/(model+'-covariance.npz')))
    np.testing.assert_array_equal(f.generated_attempt_index,cv['generated_attempt_index'])
    f['u']=f.mB_double-f.MU;f['x']=f.x1_double;f['c']=f.c_double
    f['z']=f.zHD;f['field']=f.FIELD.astype(str);f['mass']=f.HOST_LOGMASS;f['masserr']=f.HOST_LOGMASS_ERR;f['host']=(f.mass>=10).astype(int)
    f['sigma_m']=np.sqrt(cv['mag_covariance'][:,0,0]);f['sigma_x']=np.sqrt(cv['mag_covariance'][:,1,1]);f['sigma_c']=np.sqrt(cv['mag_covariance'][:,2,2]);f['fitprob']=f.FITPROB
    good=quality(f,strict)&(np.linalg.eigvalsh(cv['mag_covariance']).min(axis=1)>0)&(f.mass>0)
    return f.loc[good].reset_index(drop=True)

def raw_groups(real, sim, bw):
    groups=[];support=np.zeros((len(real),2))
    for (field,host), rows in real.groupby(['field','host']).groups.items():
        ri=np.array(rows);si=np.flatnonzero((sim.field==field)&(sim.host==host))
        if not len(si):continue
        dz=(real.z.to_numpy()[ri,None]-sim.z.to_numpy()[si][None,:])/bw
        w=np.exp(-.5*dz*dz)*(abs(dz)<=3)
        sw=w.sum(1);nw=np.count_nonzero(w,axis=1);ess=sw*sw/np.maximum((w*w).sum(1),1e-300)
        support[ri]=np.column_stack([nw,ess]);groups.append((ri,si,w))
    return groups,support

def quantile(v,w,p):
    idx=np.argsort(v); c=np.cumsum(w[idx]);return float(v[idx[np.searchsorted(c,p,side='left').clip(0,len(v)-1)]])

class Forecast:
    def __init__(self, real, sim, bw):
        self.r=real;self.s=sim;self.groups,self.support=raw_groups(real,sim,bw)
        self.B=basis(real.z);self.test=np.flatnonzero(real.fold==0);self.train=np.flatnonzero(real.fold!=0)
        self.testmap={i:k for k,i in enumerate(self.test)}
        self.values=sim[['u','c','x']].to_numpy();self.values2=self.values[:,1:]/[.1,1];self.values3=self.values/[.3,.1,1]
        self.cache=[]
        for ri,si,w in self.groups:
            ts=np.array([k for k,i in enumerate(ri) if i in self.testmap],dtype=int)
            oi=np.array([self.testmap[ri[k]] for k in ts],dtype=int)
            self.cache.append((ri,si,w,ts,oi,cdist(self.values2[si],self.values2[si]),cdist(self.values3[si],self.values3[si]),[cdist(self.values[si,j,None],self.values[si,j,None]) for j in range(3)]))

    def run(self, sim_mult=None, train_mult=None, details=False):
        mean=np.empty(len(self.r));ww=[]
        for ri,si,w,ts,oi,*dist in self.cache:
            w=w.copy()
            if sim_mult is not None:w*=sim_mult[self.s.generated_attempt_index.to_numpy()[si]-1][None,:]
            sw=w.sum(1)
            if np.any(sw==0):raise ValueError('Bootstrap removes every support member')
            w/=sw[:,None];mean[ri]=w@self.values[si,0];ww.append(w)
        tr=self.train;tw=np.ones(len(tr)) if train_mult is None else train_mult
        fit=np.linalg.lstsq(self.B[tr]*np.sqrt(tw[:,None]),(self.r.u.to_numpy()[tr]-mean[tr])*np.sqrt(tw),rcond=None)
        offset=self.B@fit[0]
        metrics={k:np.zeros(len(self.test)) for k in ['energy_cx','energy_mcx','crps_m','crps_c','crps_x','delta_c','delta_x','delta_c2','delta_x2','delta_red_fraction','delta_sigma_m','delta_sigma_x','delta_sigma_c','delta_fitprob']}
        if details:
            for dim in ['m','c','x']:
                for suffix in ['pit','coverage68','coverage95']:metrics[dim+'_'+suffix]=np.zeros(len(self.test))
        for cached,w in zip(self.cache,ww):
            ri,si,_,ts,oi,d2,d3,ds=cached
            if not len(ts):continue
            wt=w[ts];rv=self.r.iloc[ri[ts]][['u','c','x']].to_numpy().copy();rv[:,0]-=offset[ri[ts]];sv=self.values[si]
            metrics['energy_cx'][oi]=(wt*cdist(rv[:,1:]/[.1,1],sv[:,1:]/[.1,1])).sum(1)-.5*((wt@d2)*wt).sum(1)
            metrics['energy_mcx'][oi]=(wt*cdist(rv/[.3,.1,1],sv/[.3,.1,1])).sum(1)-.5*((wt@d3)*wt).sum(1)
            for j,dim in enumerate(['m','c','x']):
                metrics['crps_'+dim][oi]=(wt*np.abs(rv[:,j,None]-sv[None,:,j])).sum(1)-.5*((wt@ds[j])*wt).sum(1)
                if details:
                    metrics[dim+'_pit'][oi]=(wt*(sv[None,:,j]<=rv[:,j,None])).sum(1)
                    for kk in range(len(ts)):
                        for level,lo,hi in [(68,.16,.84),(95,.025,.975)]:
                            metrics[dim+f'_coverage{level}'][oi[kk]]=quantile(sv[:,j],wt[kk],lo)<=rv[kk,j]<=quantile(sv[:,j],wt[kk],hi)
            metrics['delta_c'][oi]=rv[:,1]-wt@sv[:,1];metrics['delta_x'][oi]=rv[:,2]-wt@sv[:,2]
            metrics['delta_c2'][oi]=rv[:,1]**2-wt@(sv[:,1]**2);metrics['delta_x2'][oi]=rv[:,2]**2-wt@(sv[:,2]**2)
            metrics['delta_red_fraction'][oi]=(rv[:,1]>.1)-wt@(sv[:,1]>.1)
            for col in ['sigma_m','sigma_x','sigma_c','fitprob']:
                metrics['delta_'+col][oi]=self.r[col].to_numpy()[ri[ts]]-wt@self.s[col].to_numpy()[si]
        return metrics,{'knots':KNOTS.tolist(),'offsets':fit[0].tolist(),'rank':int(fit[2]),'training_residual_rms':float(np.sqrt(np.average((self.r.u.to_numpy()[tr]-mean[tr]-offset[tr])**2,weights=tw)))}

def field_interval(values,fields):
    # Equal field weights expose one-field dependence; not the object-weighted estimand.
    v=pd.DataFrame({'v':values,'f':fields}).groupby('f').v.mean().to_numpy();se=v.std(ddof=1)/np.sqrt(len(v));h=t.ppf(.975,len(v)-1)*se
    return {'mean_equal_field':float(v.mean()),'ci95':[float(v.mean()-h),float(v.mean()+h)],'n_fields':len(v)}

def run(name, arm='conditioned-multistart-best', all_des=False, strict=False,bw=.10,q=0,nboot=300):
    dest=OUT/name;dest.mkdir(parents=True,exist_ok=True)
    real,counts=real_data(arm,all_des,strict,q);sims={m:sim_data(m,strict) for m in MODELS}
    good=np.ones(len(real),dtype=bool);support=real[['CID','fold','field','z','mass']].copy()
    for m,s in sims.items():
        _,n=raw_groups(real,s,bw);support[m+'_n']=n[:,0];support[m+'_ess']=n[:,1];good&=(n[:,0]>=10)&(n[:,1]>=8)
    support['common_support']=good;support.to_csv(dest/'support.csv',index=False)
    real=real.loc[good].reset_index(drop=True);real.to_csv(dest/'real-cohort.csv',index=False)
    test=real.loc[real.fold==0].reset_index(drop=True);ntrain=int((real.fold!=0).sum());ntest=len(test)
    assert ntrain>=20 and ntest>=10
    forecasts={m:Forecast(real,s,bw) for m,s in sims.items()};base={};align={}
    for m,f in forecasts.items():
        v,a=f.run(details=True);base[m]=v;align[m]=a;pd.DataFrame({'CID':test.CID,**v}).to_csv(dest/(m+'-test-scores.csv'),index=False)
    keys=list(forecasts['P21'].run()[0]);boot={m:np.empty((nboot,len(keys))) for m in MODELS};rng=np.random.default_rng(2026092141)
    for b in range(nboot):
        sm=rng.poisson(1,26518);tr=rng.multinomial(ntrain,np.full(ntrain,1/ntrain));te=rng.multinomial(ntest,np.full(ntest,1/ntest))
        for m,f in forecasts.items():
            v,_=f.run(sm,tr);boot[m][b]=[np.average(v[k],weights=te) for k in keys]
        if b%50==0:print(f'{name}: bootstrap {b}/{nboot}',flush=True)
    reports={}
    for m in MODELS:
        scores={}
        for i,k in enumerate(keys):
            val=float(np.mean(base[m][k]));contrast=base[m][k]-base['P21'][k]
            scores[k]={'mean':val,'delta_vs_P21':float(contrast.mean()),'field_contrast':field_interval(contrast,test.field)}
            if nboot:
                scores[k].update(ci95=np.quantile(boot[m][:,i],[.025,.975]).tolist(),delta_vs_P21_ci95=np.quantile(boot[m][:,i]-boot['P21'][:,i],[.025,.975]).tolist())
        reports[m]={'n_sim_quality':len(sims[m]),'alignment':align[m],'metrics':scores,'coverage_and_pit':{k:float(v.mean()) for k,v in base[m].items() if k not in keys}}
    np.savez_compressed(dest/'bootstrap.npz',metric_names=np.array(keys),**boot)
    summary={'configuration':{'arm':arm,'all_des':all_des,'strict':strict,'bandwidth':bw,'reference_q':q,'bootstrap_replicates':nboot},'cohort':{**counts,'n_train_common':ntrain,'n_test_common':ntest,'n_fields_test':int(test.field.nunique())},'models':reports,'scope':'Retrospective finite-ensemble conditional forecast comparison; incomplete original selection; no physical-model posterior or independent heldout population validation.'}
    (dest/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({'name':name,**summary['cohort']}),flush=True)
    return summary

def manifest():
    files=[source(Path(__file__)),source(OUT/'preregistration.json'),source(ROOT/'phase2/hierarchy/folds.csv'),source(FWD/'manifest.json')]
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    out={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'inputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(INPUTS)},'outputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json'}}
    (OUT/'manifest.json').write_text(json.dumps(out,indent=2)+'\n')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--bootstrap',type=int,default=300);ap.add_argument('--sensitivities',action='store_true');a=ap.parse_args()
    summary=json.loads(source(FWD/'summary.json').read_text());assert all(m['fit_complete'] for m in summary['models'])
    run('primary',nboot=a.bootstrap)
    if a.sensitivities:
        for name,kw in [('all-des',{'all_des':True}),('recovered-mask',{'arm':'recovered-mask'}),('strict-cuts',{'strict':True}),('bandwidth075',{'bw':.075}),('bandwidth150',{'bw':.15}),('reference-qplus05',{'q':.5}),('reference-qminus1',{'q':-1})]:run(name,nboot=a.bootstrap,**kw)
    manifest()

if __name__=='__main__':main()
