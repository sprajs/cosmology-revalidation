import ast, hashlib, io, json, math, sys, tarfile, time
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
from scipy.integrate import quad
ROOT=Path('/home/szymon/Documents/ChatGPT/supernova')
OUT=ROOT/'runs/assumption_audit/pantheon'
sys.path.insert(0,str(ROOT/'scripts/assumption_audit'))
from pantheon_age_prior_validation import violations
from pantheon_age_semantics import moments

def bad_constant(s): raise ValueError(s)
def read(p): return json.loads(p.read_text(),parse_constant=bad_constant)
s=read(OUT/'r19-global-validated-age-summary.json')
v=read(OUT/'r19-global-prior-violations.json')
f=pd.read_csv(OUT/'r19-global-validated-age-audit.csv').set_index('CID')
assert np.isfinite(f.to_numpy()).all()
assert len(f)==f.index.nunique()==len(v)==103
assert set(f.index)=={r['CID'] for r in v}
photo=pd.read_csv(ROOT/'sources/repos/benjaminrose__mc-age/data/campbell_global.tsv',sep='\t').set_index('SNID')
assert set(f.index)==set(photo.index)
assert np.array_equal(f.draws_total,f.draws_valid+f.draws_excluded)
assert np.allclose(f.invalid_fraction,f.draws_excluded/f.draws_total,rtol=0,atol=1e-15)
assert s['draws_total']==int(f.draws_total.sum())==105059947==103*1020000-53
assert s['draws_valid']==int(f.draws_valid.sum())==103771926
assert s['draws_excluded']==int(f.draws_excluded.sum())==1288021
assert s['hosts_with_excluded_draws']==int((f.draws_excluded>0).sum())==86
assert f.loc[20048,'draws_total']==1019947 and (f.drop(20048).draws_total==1020000).all()
assert sum(r['archive_nonfinite_age'] for r in v)==0
assert sum(r['rows_with_negative_tau'] for r in v)==5186
for r in v:
    x=f.loc[r['CID']]
    assert (r['draws_total'],r['draws_used'],r['draws_excluded'])==(x.draws_total,x.draws_valid,x.draws_excluded)
    assert max(r['violations'].values())<=r['draws_excluded']<=sum(r['violations'].values())
for k,n in s['all_violations'].items(): assert n==sum(r['violations'][k] for r in v)
c=f.fsps_valid_age_median-f.original_valid_age_median
assert np.allclose(c,f.median_age_change_valid,rtol=0,atol=3e-15)
assert np.allclose(f.original_valid_age_median-f.original_full_age_median,f.median_effect_of_removing_invalid,atol=3e-15,rtol=0)
g=s['age_change_conditional_on_valid_draws']
for k,z in {'mean':c.mean(),'median':c.median(),'min':c.min(),'max':c.max(),'n_abs_gt_p1':(abs(c)>.1).sum(),'n_abs_gt_p5':(abs(c)>.5).sum(),'n_abs_gt_1':(abs(c)>1).sum(),'original_median_below4':(f.original_valid_age_median<4).sum(),'fsps_median_below4':(f.fsps_valid_age_median<4).sum(),'n_crossing_4Gyr':((f.original_valid_age_median<4)!=(f.fsps_valid_age_median<4)).sum()}.items(): assert math.isclose(g[k],z,rel_tol=0,abs_tol=3e-15),(k,g[k],z)
source=ROOT/'sources/repos/benjaminrose__mc-age/calculateAge.py'
t=ast.parse(source.read_text()); fn=next(n for n in t.body if isinstance(n,ast.FunctionDef) and n.name=='lnprior')
ns={'np':np,'cosmo':FlatLambdaCDM(H0=70,Om0=.27,Tcmb0=0)}
exec(compile(ast.Module(body=[fn],type_ignores=[]),str(source),'exec'),ns)
prior=ns['lnprior']; checks=0
# Original-source lnprior checked at all boundaries, adjacent floats, and nonfinite parameters.
z=.15; ca=float(ns['cosmo'].age(z).value); base=np.array([-.5,.3,1.,1.,9.,0.,-25.])
thresholds={0:[-2.5,.5],1:[0.,.9],2:[.1,10.],3:[.5,7.],4:[2.5,ca],5:[-1.520838,1.520838],6:[-45.,-5.]}
for dim,edges in thresholds.items():
    for val in edges+[np.nan,np.inf,-np.inf]:
        for x in [val,np.nextafter(val,-np.inf),np.nextafter(val,np.inf)] if np.isfinite(val) else [val]:
            d=base.copy();d[dim]=x
            accepted=not any(a[0] for a in violations(d[None,:],ca).values())
            assert accepted==bool(np.isfinite(prior(d,z))), (dim,x)
            checks+=1
for r in v:
    z=float(photo.loc[r['CID'],'redshift']); ca=float(ns['cosmo'].age(z).value)
    for e in r['example_invalid_rows']:
        a=np.array(e['parameters'])
        assert not np.isfinite(prior(a,z))
        assert any(x[0] for x in violations(a[None,:],ca).values())
        checks+=1
print('Aggregate denominators, finite files, and',checks,'actual-source prior checks passed',flush=True)
# Independently re-open the raw CID20048 archive member to verify its anomalous denominator.
target='campbellG/SN20048_campbellG_chain.tsv'; raw=None
with tarfile.open(OUT/'sources/rose2019-campbellG.tar.gz',mode='r|gz') as tf:
    for member in tf:
        if member.isfile() and member.name==target:
            raw=tf.extractfile(member).read();break
assert raw is not None
rawhash=hashlib.sha256(raw).hexdigest()
r=next(x for x in v if x['CID']==20048)
assert rawhash==r['member_sha256']
d=pd.read_csv(io.BytesIO(raw),sep='\t',comment='#',header=None).to_numpy()
assert d.shape==(1019947,8) and np.isfinite(d).all()
z=float(photo.loc[20048,'redshift']);ca=float(ns['cosmo'].age(z).value)
# Test actual source, rather than vectorized copy, on a uniformly spaced 1021-row subset.
ids=np.unique(np.linspace(0,len(d)-1,1021,dtype=int))
assert all(np.isfinite(prior(d[j,:7],z)) for j in ids)
a,b,c0,e,h,p,k=d[:,:7].T
valid=(a>-2.5)&(a<.5)&(b>=0)&(b<=.9)&(c0>.1)&(c0<10)&(e>.5)&(e<h-2)&(h>2.5)&(h<=ca)&(p>-1.520838)&(p<1.520838)&(k>-45)&(k<-5)
assert valid.all()
new,_=moments(d[:,2],d[:,3],d[:,4],np.tan(d[:,5]),ca,True)
assert np.isfinite(new).all() and (new>=0).all() and (new<=ca).all()
assert np.allclose(np.quantile(new,[.16,.5,.84]),f.loc[20048,['fsps_valid_age_q16','fsps_valid_age_median','fsps_valid_age_q84']].to_numpy(),rtol=0,atol=3e-14)
# Independently integrate native FSPS normalization at 16 fixed posterior draws.
max_error=0.
for j in np.unique(np.linspace(0,len(d)-1,16,dtype=int)):
    tau,st,tr,phi=d[j,2:6];T=tr-st;end=ca-st;sl=np.tan(phi);rate=T/tau*np.exp(-T/tau)
    def sf(t): return t/tau*np.exp(-t/tau) if t<=T else max(0.,rate*(1+sl*(t-T)))
    points=[T]
    if sl<0 and T<T-1/sl<end: points.append(T-1/sl)
    mass=quad(sf,0,end,points=points,epsabs=1e-12,epsrel=1e-11)[0]
    age=quad(lambda t:(end-t)*sf(t),0,end,points=points,epsabs=1e-12,epsrel=1e-11)[0]/mass
    max_error=max(max_error,abs(age-new[j]))
assert max_error<1e-8
result={'verdict':'PASS for bounded independent review','scope':'All 103 summary rows and prior-violation records; raw CID20048 plus original-source boundary/example/sample prior checks; no complete second 105-million-row parsing pass.','denominators':{'hosts':103,'draws_total':105059947,'valid':103771926,'excluded':1288021,'CID20048_rows':1019947,'CID20048_missing_from_2400x425':53},'strict_finite_json_and_csv':True,'original_lnprior_boundary_and_invalid_example_checks':checks,'CID20048_original_lnprior_sample_checks':len(ids),'CID20048_native_FSPS_quadrature_checks':16,'CID20048_max_quadrature_error_Gyr':max_error,'CID20048_member_sha256':rawhash,'caveats':['Age-shift aggregates are differences of per-host valid-subset posterior medians, with equal host weighting; they are not median paired draw shifts.','Valid support does not establish convergence, independent samples, or restored correct posterior weights; exclusions cannot repair a defective chain.','All original global hosts select Tcmb0=0 from archived-column consistency; the selection is empirical version reconstruction, not a new fitted cosmological parameter.','Source-column reproduction has rare original quadrature outliers; tail reproduction is not uniformly machine precision.','Neither the revised C25 ages nor any cosmological effect is established.'],'reviewed_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [OUT/'r19-global-validated-age-summary.json',OUT/'r19-global-validated-age-audit.csv',OUT/'r19-global-prior-violations.json',ROOT/'scripts/assumption_audit/pantheon_age_prior_validation.py',ROOT/'scripts/assumption_audit/pantheon_age_semantics.py',source]}}
path=ROOT/'runs/assumption_audit/des/r19-prior-peer-review.json';path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps(result,indent=2,allow_nan=False),flush=True)
