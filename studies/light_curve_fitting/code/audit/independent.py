"""Independent implementation audit; no import of investigator calculation modules."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
from pathlib import Path
import json, hashlib, subprocess, datetime
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp, quad
from scipy.interpolate import PchipInterpolator
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import minimize, brentq

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/audit'; OUT.mkdir(exist_ok=True, parents=True)
INPUTS=[]
def track(p):
    INPUTS.append(p)
    return p

def ode_distance(theta, z):
    """Evolve rho_DE/rho_DE0 and chi, avoiding explicit CPL density formula."""
    om,w0,wa=theta
    def rhs(v, y):
        h=np.sqrt(om*(1+v)**3+(1-om)*y[0])
        return [3*(1+w0+wa*v/(1+v))/(1+v)*y[0], 1/h]
    sol=solve_ivp(rhs, [0,float(np.max(z))], [1.,0.],method='DOP853',rtol=2e-11,atol=1e-12,dense_output=True)
    assert sol.success
    rho,chi=sol.sol(z)
    return chi, np.sqrt(om*(1+z)**3+(1-om)*rho)

def raw_data():
    b=ROOT/'sources/repos/CobayaSampler__sn_data/PantheonPlus'
    df=pd.read_csv(track(b/'Pantheon+SH0ES.dat'),sep=r'\s+',dtype={'CID':str})
    v=np.loadtxt(track(b/'Pantheon+SH0ES_STAT+SYS.cov')); n=int(v[0]);assert len(df)==n
    cov=v[1:].reshape(n,n); mask=df.zHD.to_numpy()>.01
    cov=(cov+cov.T)/2
    sub=df[mask].copy();csub=cov[np.ix_(mask,mask)];L=cholesky(csub,lower=True)
    u=solve_triangular(L,np.ones(len(sub)),lower=True)
    bw=ROOT/'sources/repos/CobayaSampler__bao_data/desi_bao_dr2'
    bao=pd.read_csv(track(bw/'desi_gaussian_bao_ALL_GCcomb_mean.txt'),sep=r'\s+',comment='#',names=['z','value','kind'])
    LB=cholesky(np.loadtxt(track(bw/'desi_gaussian_bao_ALL_GCcomb_cov.txt')),lower=True)
    return df,cov,sub,L,u,bao,LB

def main():
    df,C,sn,L,u,bao,LB=raw_data();results={}
    zz=np.r_[sn.zHD,bao.z];ns=len(sn)
    def likelihood(t, correction):
        chi,h=ode_distance(t[:3],zz)
        model=5*np.log10((1+sn.zHEL.to_numpy())*chi[:ns])
        obs=sn.m_b_corr.to_numpy()-correction
        whitened=solve_triangular(L,obs-model,lower=True,check_finite=False)
        offset=whitened@u/(u@u);r=whitened-offset*u
        dm=299792.458/(10000*t[3])*chi[ns:]
        dh=299792.458/(10000*t[3])/h[ns:]
        dv=(bao.z.to_numpy()*dm**2*dh)**(1/3)
        pred=np.select([bao.kind=='DM_over_rs',bao.kind=='DH_over_rs',bao.kind=='DV_over_rs'],[dm,dh,dv],default=np.nan)
        assert np.isfinite(pred).all()
        br=solve_triangular(LB,bao.value.to_numpy()-pred,lower=True)
        return float(r@r),float(br@br),float(offset)
    fits=[]
    for name in ['pantheon-bao-cpl','pantheon-bao-cpl-c14fixed','pantheon-bao-cpl-c14mean']:
        old=json.loads(track(ROOT/'runs/cosmology'/name/'summary.json').read_text())
        conf=json.loads(track(ROOT/'runs/cosmology'/name/'configuration.json').read_text())
        correction=np.zeros(ns)
        if conf['correction']:
            tab=pd.read_csv(track(ROOT/conf['correction']))
            correction=PchipInterpolator(tab.z,tab.delta_mu)(sn.zHD.to_numpy())
        orig=np.array(min(old['optimizer_starts'],key=lambda x:x['objective'])['x']);start=orig.copy();start[3]/=10000
        comp=likelihood(start,correction)
        bounds=[(.01,.99),(-3,1),(-3,2),(.5,1.5)]
        fit=minimize(lambda t:sum(likelihood(t,correction)[:2]) if t[1]+t[2]<0 else 1e20,
                     start+np.array([.001,.004,-.008,.001]),method='Nelder-Mead',bounds=bounds,
                     options={'maxiter':1600,'xatol':1e-8,'fatol':1e-8})
        restored=fit.x.copy();restored[3]*=10000
        q=lambda t: .5+1.5*t[1]*(1-t[0])
        rec={'experiment':name,'independent_at_original_mode_sn_bao_offset':comp,
             'original_mode_sn_bao':[old['mode_sn_chisq'],old['mode_bao_chisq']],
             'total_chisq_difference_at_original_mode':sum(comp[:2])-old['posterior_mode_chisq'],
             'independent_optimum':restored.tolist(),'original_optimum':orig.tolist(),'independent_optimum_chisq':float(fit.fun),
             'independent_optimum_q0':float(q(fit.x)),'original_optimum_q0':float(q(orig)),
             'q0_optimum_difference':float(q(fit.x)-q(orig)),'optimizer_success':bool(fit.success)}
        fits.append(rec);print(name,rec,flush=True)
    results['cosmology_likelihood']=fits
    # Direct ODE-versus-stored-template distance check at arbitrary extreme parameters.
    results['distance_ode_analytic_limits']={}
    z=np.array([.01,.1,1,2.3]); chi,_=ode_distance([1.,-1.,0.],z)
    results['distance_ode_analytic_limits']['matter_only_max_dimensionless_error']=float(abs(chi-2*(1-(1+z)**-.5)).max())
    chi,_=ode_distance([0.,-1.,0.],z)
    results['distance_ode_analytic_limits']['de_sitter_max_dimensionless_error']=float(abs(chi-z).max())
    results['mapping']=mapping()
    results['age']=age_audit(df,C)
    results['mock']=mock()
    results['posterior_prior_sensitivities']=posterior()
    path=OUT/'independent-results.json';path.write_text(json.dumps(results,indent=2)+'\n')
    source_audit_files=[
        'papers/text/chung2026-published.txt','papers/text/son2025-published.txt','papers/text/wiseman2026-published.txt',
        'runs/mapping/sources-2026-09-20/behroozi2013.txt',
        'sources/updates/2026-09-20-standardization/2601.19424.txt',
        'sources/updates/2026-09-20-standardization/2607.24443.txt',
        'scripts/cosmology/core.py','scripts/cosmology/run.py','scripts/cosmology/diagnostics.py',
        'scripts/age_signal/sensitivity.py','scripts/age_signal/fullcov_eiv.py',
        'scripts/mapping/csfh_dtd.py','scripts/standardization/compare_des.py',
        'scripts/standardization/pantheon_mass_update.py','runs/cosmology/diagnostics/results.json']
    inputs=set(INPUTS+[Path(__file__),ROOT/'docs/experiments/independent-audit-plan.md',ROOT/'uv.lock']+[ROOT/p for p in source_audit_files])
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    manifest={'time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'purpose':'Independent falsification and numerical audit',
              'git_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
              'inputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(inputs)},
              'outputs_sha256':{str(path.relative_to(ROOT)):sha(path)},'seed':38192,
              'implementations':'CPL rhoDE+chi ODE; explicit whitened intercept; formation-redshift DTD quadrature; raw ID parsing'}
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

def mapping():
    # Analytic flat matter+Lambda cosmic clock and integration over formation z.
    H0=70.;om=.3;clock=lambda z:2*977.7922216807892/(3*H0*np.sqrt(1-om))*np.arcsinh(np.sqrt((1-om)/om)/(1+z)**1.5)
    sfr=lambda z:.180/(10**(-.997*(z-1.243))+10**(.241*(z-1.243)))
    dtds=['C14_smooth','W26_cut40'];res={}
    published=pd.read_csv(track(ROOT/'runs/mapping/csfh-dtd-curves.csv'))
    for dtd in dtds:
        rec=[]
        for z in [0.,1.]:
            t=clock(z)
            lo=z if dtd=='C14_smooth' else brentq(lambda x:t-clock(x)-.04,z,100.)
            def pdf(zf):
                delay=t-clock(zf)
                if delay<=0:return 0.
                if dtd=='C14_smooth':
                    x=np.log(delay/.3);phi=np.exp(20*x-np.logaddexp(21*x,0))
                else:phi=delay**-1.13
                jac=977.7922216807892/H0/((1+zf)*np.sqrt(om*(1+zf)**3+1-om))
                return sfr(zf)*phi*jac
            norm=quad(pdf,lo,100,epsabs=1e-10,epsrel=1e-10,limit=300)[0]
            mean=quad(lambda v:(t-clock(v))*pdf(v),lo,100,epsabs=1e-10,limit=300)[0]/norm
            med_z=brentq(lambda v:quad(pdf,lo,v,epsabs=1e-10,limit=300)[0]-.5*norm,lo,100,xtol=1e-12)
            median=t-clock(med_z)
            p=published[(published.cosmology=='lcdm_H70')&(published.csfh=='B13')&(published.dtd==dtd)&(published.z==z)].iloc[0]
            rec.append({'z':z,'mean':mean,'median':median,'mean_difference_gyr':mean-p.mean_delay_gyr,'median_difference_gyr':median-p.median_delay_gyr})
        res[dtd]={'values':rec,'mean_evolution_gyr':rec[0]['mean']-rec[1]['mean'],'median_evolution_gyr':rec[0]['median']-rec[1]['median'],
                  'median_correction_mag':.03*(rec[0]['median']-rec[1]['median'])}
    return res

def age_audit(df,C):
    ages=[]
    for n in [1,2]:
        p=track(ROOT/f'data/host_ages/chung2025/table{n}.dat');rows=[]
        for line in p.read_text().splitlines():
            fields=line.replace('\\','').split('&')
            try: rows.append({'CID':str(int(fields[0])),'age':float(fields[1]),'age_error':float(fields[2])})
            except (ValueError,IndexError):pass
        ages.append(pd.DataFrame(rows))
    overlap=set(ages[0].CID)&set(ages[1].CID)
    df=df.copy();df['row']=np.arange(len(df))
    res={'age_counts':[len(d) for d in ages],'overlap':len(overlap),'variants':[]}
    for policy in ['G11_first','R19_first']:
        d=pd.concat(ages).drop_duplicates('CID',keep='first' if policy=='G11_first' else 'last').merge(df[df.IDSURVEY==1],on='CID',validate='one_to_one')
        d=d[(d.zHD>.06)&(d.zHD<.42)];idx=d.row.to_numpy();c=C[np.ix_(idx,idx)];L=cholesky(c,lower=True)
        a=d.age.to_numpy();z=d.zHD.to_numpy();zh=d.zHEL.to_numpy();X=np.column_stack([np.ones(len(a)),a]);wx=solve_triangular(L,X,lower=True)
        chi,_=ode_distance([.3,-1,0],z);mu=5*np.log10((1+zh)*chi)
        hr=d.m_b_corr.to_numpy()-mu
        fits=[]
        for correction,y in [('corrected',hr),('bias_reversed',hr+d.biasCor_m_b.to_numpy())]:
            wy=solve_triangular(L,y,lower=True);coef=np.linalg.lstsq(wx,wy,rcond=None)[0];r=wy-wx@coef
            fits.append({'correction':correction,'fixed_age_slope':float(coef[1]),'chi2':float(r@r)})
        # Within the stated fixed-age diagnostic, flexible cosmology can absorb more than a single Om derivative.
        center=lambda v:v-u*(u@v)/(u@u)
        u=solve_triangular(L,np.ones(len(a)),lower=True);v=center(solve_triangular(L,a,lower=True))
        derivatives=[]
        theta=np.array([.3,-1,0.]);eps=1e-4
        for k in range(3):
            tp=theta.copy();tm=theta.copy();tp[k]+=eps;tm[k]-=eps
            cp,_=ode_distance(tp,z);cm,_=ode_distance(tm,z)
            derivatives.append(center(solve_triangular(L,5*np.log10(cp/cm)/(2*eps),lower=True)))
        J=np.column_stack(derivatives)
        fractions={}
        for n in [1,3]:
            fitted=J[:,:n]@np.linalg.lstsq(J[:,:n],v,rcond=None)[0]
            fractions[str(n)+'_cosmology_derivatives']=float(fitted@fitted/(v@v))
        res['variants'].append({'policy':policy,'n':len(d),'cut_counts':[int((z<cut).sum()) for cut in [.2,.25,.3,.35,.42]],
              'table_covdiag_rms_mag':float(np.sqrt(np.mean((d.MU_SH0ES_ERR_DIAG.to_numpy()-np.sqrt(np.diag(c)))**2))),
              'fits':fits,'observed_age_variation_projected_fraction':fractions})
    return res

def mock():
    rng=np.random.default_rng(38192);z=np.repeat(np.arange(.05,.401,.05),200);age=7-10*z+rng.normal(0,1.5,len(z));b=-.034;y=b*age
    centred=y.copy();ac=age.copy()
    for zz in np.unique(z):
        m=z==zz;centred[m]-=y[m].mean();ac[m]-=age[m].mean()
    slope=lambda x,y:float(np.cov(x,y,ddof=0)[0,1]/np.var(x))
    # Project the luminosity evolution on a smooth fitted baseline. Finite-sample ratio still gives attenuation.
    Z=np.column_stack([np.ones(len(z)),z]);fit=Z@np.linalg.lstsq(Z,y,rcond=None)[0]
    return {'b_true':b,'fixed_correct_baseline_slope':slope(age,y),
            'bin_centred_only_slope':slope(age,centred),'within_total_prediction':float(b*np.var(ac)/np.var(age)),
            'both_centred_slope':slope(ac,centred),'fitted_linear_z_baseline_slope':slope(age,y-fit),
            'n':len(age),'noise':0,'purpose':'Noiseless mathematical counterexample; not an empirical mock reproduction'}

def posterior():
    res=[]
    for name in ['bao-cpl','pantheon-cpl','pantheon-bao-cpl','pantheon-bao-cpl-c14fixed','pantheon-bao-cpl-c14slope']:
        p=track(ROOT/'runs/cosmology'/name/'chains.npz');f=np.load(p);t=f['chain'].reshape(-1,f['chain'].shape[-1]);names=f['names'].tolist()
        q=.5+1.5*t[:,1]*(1-t[:,0]);case={'experiment':name,'draws':len(t),'sample_fraction_q0_negative':float((q<0).mean()),
             'q0_mean':float(q.mean()),'q0_sd':float(q.std()),'fraction_w0_plus_wa_above_minus_point05':float((t[:,1]+t[:,2]>-.05).mean()),
             'fraction_wa_below_minus2point95':float((t[:,2]<-2.95).mean()),'reweighted':{}}
        weights={'flat_q0_instead_of_flat_w0':1-t[:,0]}
        if 'H0_rd' in names:weights['log_H0rd_instead_of_flat_H0rd']=1/t[:,names.index('H0_rd')]
        for label,w in weights.items():
            w=w/w.sum();m=w@q
            case['reweighted'][label]={'mean_q0':float(m),'sd_q0':float(np.sqrt(w@((q-m)**2))),'fraction_q0_negative':float(w@(q<0)),
                'importance_ESS_fraction':float(1/(w@w)/len(w)), 'same_original_parameter_support':True}
        res.append(case)
    return res

if __name__=='__main__':main()
