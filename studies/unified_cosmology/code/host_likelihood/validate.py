"""Independent arithmetic, observational identities and explicit failed gates."""
from __future__ import annotations
import ast,hashlib,json,re
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from astropy.cosmology import WMAP9
from scipy.optimize import minimize
from scipy.special import logsumexp
from acquire import ROOT,WORK,OUT,sha
from likelihood import host_flux_loglike,conditional_brightness_loglike,grey_compensator
from ozdes_bands import ratio_profile
from ozdes_repeat import common_ratio


def main():
    code=Path(__file__).parent;identities={};checks={}
    def identity(path,digest):
        path=ROOT/path if not Path(path).is_absolute() else Path(path)
        key=str(path.relative_to(ROOT));actual=identities.get(key)
        if actual is None:actual=sha(path);identities[key]=actual
        assert actual==digest,(key,actual,digest)
    def walk(node):
        if isinstance(node,dict):
            if isinstance(node.get('path'),str) and isinstance(node.get('sha256'),str):identity(node['path'],node['sha256'])
            for k,v in node.items():
                if k.endswith('_sha256') and isinstance(v,dict):
                    for p,digest in v.items():
                        if isinstance(digest,str) and len(digest)==64:identity(p,digest)
                walk(v)
        elif isinstance(node,list):
            for item in node:walk(item)
    mapping={'acquisition':'acquire','additional-acquisition':'recover','host-interface':'extract','primary-source-audit':'source_audit','yse-flux-fits':'yse_fit','yse-fit-audit':'yse_audit','ozdes-acquisition':'ozdes_acquire','ozdes-band-interface':'ozdes_bands','ozdes-repeatability':'ozdes_repeat'}
    records={}
    for name,script in mapping.items():
        p=OUT/(name+'.json');r=json.loads(p.read_text());identity(code/(script+'.py'),r['code_sha256']);walk(r);records[name]=r
        identity(p,sha(p))
    for name,design in {'acquisition':'design.json','yse-flux-fits':'yse-design.json','ozdes-acquisition':'ozdes-design.json','ozdes-band-interface':'ozdes-design.json','ozdes-repeatability':'ozdes-repeat-design.json'}.items():
        identity(code/design,records[name]['design_sha256'])
    identity(WORK/'frankenblast-record.json',records['acquisition']['record_sha256'])
    for p in code.glob('*.py'):ast.parse(p.read_text());identity(p,sha(p))
    for p in code.glob('*.json'):json.loads(p.read_text());identity(p,sha(p))
    summary=pd.read_csv(WORK/'host-summary.csv').set_index('name');cross=pd.read_csv(WORK/'yse-host-crosswalk.csv');primary=cross[cross.primary_candidate]
    assert len(summary)==4385 and len(primary)==66 and len(records['host-interface']['invalid'])==114
    sfh_error=[];cov_error=[]
    for name in sorted(primary.name)[:12]:
        f=np.load(WORK/'joint-host-draws'/(name+'.npz'));chain=f['raw_chain'][:24];labels=list(f['theta_labels']);age=float(WMAP9.age(summary.loc[name,'obs_z']).value)
        edges=np.r_[1e-9,10**7.4772/1e9,np.geomspace(.1,.9*age,5),age]
        ratios=10**chain[:,[labels.index('logsfr_ratios_'+str(i))for i in range(1,7)]]
        sf=np.ones((len(chain),7))
        for j in range(1,7):sf[:,j]=sf[:,j-1]/ratios[:,j-1]
        mass=sf*np.diff(edges);mass/=mass.sum(axis=1,keepdims=True);direct=mass@((edges[:-1]+edges[1:])/2)
        sfh_error.append(float(np.max(abs(direct-f['physical_draws'][:24,0]))));cov_error.append(float(np.max(abs(np.cov(f['physical_draws'],rowvar=False)-f['covariance']))))
    assert max(sfh_error)<1e-12 and max(cov_error)<1e-12
    checks['same_row_host_draws']=dict(hosts=12,draws_per_host=24,maximum_age_difference_Gyr=max(sfh_error),maximum_covariance_difference=max(cov_error))
    rng=np.random.default_rng(927084);h=rng.normal(size=(5,100,2));r=rng.normal(size=5);v=np.linspace(.01,.1,5);b=np.array([.02,-.04])
    result=conditional_brightness_loglike(r,v,h,b,host_model_unchanged=True);brute=[]
    for i in range(5):brute.append(float(np.log(np.mean(np.exp(-.5*(r[i]-h[i]@b)**2/v[i])/np.sqrt(2*np.pi*v[i])))))
    assert np.allclose(result['per_object'],brute,rtol=1e-12,atol=1e-12)
    none=conditional_brightness_loglike(r,v,h,b,log_population_ratio=np.full((5,100),-np.inf));assert none['loglike']==-np.inf and np.all(none['population_ratio_ess']==0)
    y=np.array([-1.,2.]);mu=np.array([0.,1.]);cov=np.array([[2.,.3],[.3,1.]])
    direct=-.5*((y-mu)@np.linalg.solve(cov,y-mu)+np.linalg.slogdet(cov)[1]+2*np.log(2*np.pi));assert abs(direct-host_flux_loglike(y,mu,cov))<1e-12
    z=np.linspace(.01,1.1,1088);a=43+5*np.log10(z);trial=a+.2*z/(1+z);grey,offset=grey_compensator(trial,a,np.ones(len(z)));assert np.max(abs(trial+grey+offset-a))<1e-12
    checks['likelihood_algebra']=dict(joint_draw_integral_cases=5,negative_flux_full_covariance=True,zero_support_returns_minus_infinity_and_zero_ESS=True,grey_compensator_max_error=float(np.max(abs(trial+grey+offset-a))),grey_demo='Algebraic identity only, not empirical bound or physical cosmology fit')
    archive=np.load(WORK/'ozdes-band-likelihood.npz');targets=list(archive['target']);assert len(targets)==len(set(targets))==1088
    full=json.loads((WORK/'ozdes-band-records.json').read_text());ratio_checks=[];covchecks=[]
    selected=sorted(range(1088),key=lambda i:hashlib.sha256(targets[i].encode()).hexdigest())[:20]
    for i in selected:
        target=targets[i];f=np.load(WORK/'ozdes-band-kernels'/(target+'.npz'));w=f['operator'];pixel=f['pixel_flux'];var=f['pixel_variance']
        value=np.array([sum(float(x)*float(y)for x,y in zip(row,pixel))for row in w]);c=np.zeros((7,7))
        for j in range(7):
            for k in range(7):c[j,k]=sum(float(x)*float(y)*float(t)for x,y,t in zip(w[j],w[k],var))
        assert np.allclose(value,archive['flux'][i],rtol=1e-12,atol=1e-12) and np.allclose(c,archive['covariance'][i,0],rtol=1e-12,atol=1e-12)
        covchecks.append(float(np.max(abs(c-archive['covariance'][i,0]))))
        ix=[1,2];flux=value[ix];cc=c[np.ix_(ix,ix)];inv=np.linalg.inv(cc)
        native,delta=ratio_profile(flux,cc,archive['ratio_grid']);amplitude=native['amplitude_at_mle'];R=native['mle']
        opt=minimize(lambda x:float((flux-x[0]*np.array([1.,x[1]]))@inv@(flux-x[0]*np.array([1.,x[1]]))),[amplitude,R],method='Nelder-Mead',bounds=[(0,None),(0,4)],options={'xatol':1e-10,'fatol':1e-10,'maxiter':5000})
        # Grid minimum is an intentionally discretized approximation. Compare its
        # exact quadratic objective, and allow <=half a grid spacing in R.
        assert opt.fun<=native['minimum_chi2']+1e-7
        assert abs(opt.x[1]-R)<=.001251 or amplitude==0 or R in [0.,4.]
        ratio_checks.append(dict(target=target,grid_minus_continuous_minimum=float(native['minimum_chi2']-opt.fun),ratio_difference=float(R-opt.x[1])))
    checks['observed_spectral_vectors']=dict(distinct_hosts=1088,independent_pixel_sum_cases=20,maximum_covariance_absolute_difference=max(covchecks),ratio_profiles=ratio_checks)
    hd=json.loads((WORK/'ozdes-Hdelta-records.json').read_text());hdchecks=[]
    for i in selected:
        target=targets[i];f=np.load(WORK/'ozdes-repeat'/(target+'-Hdelta.npz'));op=f['operator'];wave=f['rest_wavelength_A'];row=next(x for x in hd if x['target']==target)
        C=op@np.diag(f['variance'])@op.T;obs=op@f['flux'];assert np.allclose(obs,[row['continuum'],row['absorption_contrast']],rtol=1e-12,atol=1e-12)
        assert abs(C[0,1]-row['covariance_continuum_contrast'])<1e-10
        assert abs(op[1]@np.ones(len(wave)))<1e-12 and abs(op[1]@wave)<1e-8
        hdchecks.append(float(C[0,1]))
    checks['local_continuum_projection']=dict(cases=20,constant_and_linear_continuum_cancel=True,maximum_absolute_continuum_contrast_covariance=max(abs(x)for x in hdchecks))
    repeated=json.loads((WORK/'ozdes-repeat-records.json').read_text());rr=[r for r in repeated if r['usable_exposures']>=2];repeat=records['ozdes-repeatability']
    assert len(rr)==repeat['hosts_with_repeated_usable_exposures']==1081
    assert sum(r['nominal_tail']<.01 for r in rr)==repeat['nominal_chi2_tail_below001']==439
    for row in repeat['bootstrap'][:8]:
        f=np.load(WORK/'ozdes-repeat'/(row['target']+'-exposures.npz'));fit=common_ratio(f['flux'],f['covariance']);inv=np.linalg.inv(f['covariance']);grid=np.linspace(0,4,4001);d=np.c_[np.ones(len(grid)),grid]
        denom=np.einsum('gi,nij,gj->gn',d,inv,d);num=np.einsum('gi,nij,nj->gn',d,inv,f['flux']);amp=np.maximum(0,num/denom)
        resid=f['flux'][None,:,:]-amp[:,:,None]*d[:,None,:];vals=np.einsum('gni,nij,gnj->g',resid,inv,resid)
        assert fit['chi2']<=vals.min()+1e-6 and abs(fit['ratio']-grid[vals.argmin()])<.00101
    checks['observed_repeatability']=dict(hosts=1081,individual_exposures=repeat['total_usable_individual_exposures'],independent_dense_grid_cases=8,nominal_tail_below001=439,parametric_hosts=len(repeat['bootstrap']),parametric_draws=sum(r['replicates']for r in repeat['bootstrap']),parametric_tail_below001=repeat['bootstrap_tail_below001'])
    yse=records['yse-flux-fits'];assert all(not r['fits']['positive_truncated']['analysis_gate']for r in yse['records']if r['status']=='completed')
    checks['brightness_gate']=dict(candidates=66,completed=62,numerical_and_sampling=51,nominal_residual_screen=9,approved_analysis=0,physical_host_coefficient=None)
    note=ROOT/'studies/unified_cosmology/notes/host-likelihood.md';links=0
    for p in [code/'README.md',note]:
        if not p.exists():continue
        identity(p,sha(p))
        for target in re.findall(r'\[[^\]]*\]\(([^)]+)\)',p.read_text()):
            if target.startswith(('https://','http://','#')):continue
            path=(p.parent/target.split('#')[0]).resolve()
            # This report is written atomically at successful completion below.
            assert path.exists() or path==OUT/'validation.json',(p,target)
            links+=1
    result=dict(completed_utc=datetime.now(timezone.utc).isoformat(),status='passed',code_sha256=sha(__file__),identity_count=len(identities),identities_sha256=identities,checks=checks,local_markdown_links=links,
       scientific_status='Observed host likelihood inputs recovered and arithmetic validated. YSE brightness/host inference fails model-adequacy gate; OzDES count-spectrum response/noise is not calibrated; no empirical residual-luminosity B(z) prior or unified correction exported.')
    temporary=OUT/'validation.json.partial'
    temporary.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    temporary.replace(OUT/'validation.json')
    print(json.dumps({'status':'passed','identities':len(identities),'checks':list(checks)}))
if __name__=='__main__':main()
