#!/usr/bin/env python3
"""Paired gray/chromatic injections before a declared synthetic detection rule."""
from pathlib import Path
import sys
import json
import hashlib
from datetime import datetime,timezone
import numpy as np
import pandas as pd
import extinction
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT))
from lib.flux_engine import Engine, BUNDLE, LB
from lib.flux_fit import Case
HERE=Path(__file__).resolve().parent
OUT=ROOT/'studies/host_ages/results/age_recovery'
WORK=ROOT/'.work/age-recovery'

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def native_flux(case, x):
    p=x.copy();p[3]+=case.tref
    return case.engine.sncosmo_flux(p,case.b,case.t,case.z,case.ebv)
def selected(y,C,bands):
    snr=y/np.sqrt(np.diag(C))
    return sum(np.count_nonzero((bands==b)&(snr>5))>=2 for b in np.unique(bands))>=2

def fit(case, y, scale):
    L=np.linalg.cholesky(case.cov)
    lo=np.array([case.x0[0]-5,-5.,-.5,-10.]); hi=np.array([case.x0[0]+3,5.,.5,10.])
    lo[3]=max(lo[3], float(np.max(case.t-case.tref-(1+case.z)*case.engine.tmax)+1e-6))
    hi[3]=min(hi[3], float(np.min(case.t-case.tref-(1+case.z)*case.engine.tmin)-1e-6))
    initial=case.x0+np.array([np.log(scale),.15,-.01,.15])
    def fun(p):
        return np.r_[solve_triangular(L,y-native_flux(case,p),lower=True),case.prior_residuals(p)]
    sols=[]
    for shift in [np.zeros(4),np.array([.05,-.3,.02,-.3])]:
        sols.append(least_squares(fun,np.clip(initial+shift,lo+1e-6,hi-1e-6),bounds=(lo,hi),xtol=1e-9,ftol=1e-9,gtol=1e-7,max_nfev=160))
    s=min(sols,key=lambda v:v.fun@v.fun)
    bound=bool(np.any(np.minimum(s.x-lo,hi-s.x)<1e-4))
    return dict(success=bool(s.success),boundary=bound,chi2=float(s.fun@s.fun),two_start_chi2_gap=float(max(v.fun@v.fun for v in sols)-min(v.fun@v.fun for v in sols)),x=s.x)

def main():
    protocol=json.loads((HERE/'protocol.json').read_text())
    rng=np.random.default_rng(protocol['flux_experiment']['seed'])
    engine=Engine(); rows=[]; closure=[]
    for file in sorted(BUNDLE.glob('objective_*.npz')):
        cid=file.stem.split('_')[1]; case=Case(engine,cid)
        noises=rng.normal(size=(protocol['flux_experiment']['noise_realizations'],len(case.y)))@np.linalg.cholesky(case.cov).T
        base=case.flux(case.x0,interpolation='scipy')
        closure.append(float(np.max(abs(base-native_flux(case,case.x0))/np.maximum(abs(base),1))))
        pp=case.p0.copy()
        phase=engine.flux(np.r_[pp,.03],case.b,case.t,case.z,case.grid,'phase_colour','scipy')
        grids={b:{**g,'weights':g['weights']*10**(-.4*.03*g['f99'])} for b,g in case.grid.items()}
        dust=engine.flux(pp,case.b,case.t,case.z,grids,'salt','scipy')
        ab=.03*float(extinction.fitzpatrick99(np.array([LB]),3.1,3.1)[0])
        for scale in protocol['flux_experiment']['flux_scales']:
            means={'zero':base,'gray_plus':base*10**(-.4*.09),'gray_minus':base*10**(.4*.09),'f99_colour':dust,'phase_colour':phase,'f99_screen':dust*10**(-.4*ab)}
            for arm,m in means.items():
                for k,noise in enumerate(noises):
                    y=scale*m+noise
                    sel=bool(selected(y,case.cov,case.b))
                    row=dict(CID=cid,z=case.z,flux_scale=scale,arm=arm,draw=k,selected=sel,success=False,boundary=False)
                    if sel:
                        try:
                            f=fit(case,y,scale)
                            delta=f.pop('x')-case.x0-np.array([np.log(scale),0,0,0])
                            row.update(f)
                            row.update(dict(delta_mB=float(-2.5/np.log(10)*delta[0]),delta_x1=float(delta[1]),delta_c=float(delta[2]),delta_t0=float(delta[3])))
                            row['delta_Tripp']=row['delta_mB']+.148*row['delta_x1']-3.112*row['delta_c']
                        except Exception as e:
                            row['error']=type(e).__name__+': '+str(e)
                    rows.append(row)
        print(cid,'complete',flush=True)
    frame=pd.DataFrame(rows); frame.to_csv(WORK/'flux-event-ledger.csv',index=False)
    report={'scope':'Conditional injection into 12 previously selected DES objective cadences. Recomputed synthetic detection is not DES discovery, classification, clipping, selection normalization or BBC. Frozen covariance and trained SALT3 surface; no physical measured ages for these objects. Numerical integrators share the SALT3 surface; F99 and phase interventions change generated spectra relative to SALT inference, but there is no independent surface training.', 'F99_B_normalized_colour_warning':'f99_colour removes A_B by construction and is a signed chromatic basis, not a physical dust screen; f99_screen includes the gray dimming.', 'F99_screen_A_B_mag':ab, 'eligible_event_draws':len(frame),'physical_cadences':12,'noise_realizations':protocol['flux_experiment']['noise_realizations'],'generation_inference_max_relative_flux_difference':max(closure),'arms':[]}
    for (scale,arm),g in frame.groupby(['flux_scale','arm']):
        baseline=frame[(frame.arm=='zero')&(frame.flux_scale==scale)]
        joined=g.merge(baseline,on=['CID','draw'],suffixes=('','_zero'),validate='one_to_one')
        good=joined.selected&joined.selected_zero&joined.success&joined.success_zero&~joined.boundary&~joined.boundary_zero
        jj=joined[good]
        row=dict(flux_scale=scale,arm=arm,eligible=len(g),selected=int(g.selected.sum()),success=int(g.success.sum()),boundaries=int(g.boundary.sum()),paired_qualified=len(jj),selection_gained=int((joined.selected&~joined.selected_zero).sum()),selection_lost=int((~joined.selected&joined.selected_zero).sum()),max_two_start_chi2_gap=float(g.two_start_chi2_gap.max()))
        for col in ['delta_mB','delta_x1','delta_c','delta_t0','delta_Tripp']:
            delta=jj[col]-jj[col+'_zero']
            # Cluster each cadence before reporting variation across physical objects.
            per_cid=delta.groupby(jj.CID).mean()
            row[col]=dict(paired_mean=float(delta.mean()),noise_pair_sd=float(delta.std(ddof=1)),cadence_mean_range=[float(per_cid.min()),float(per_cid.max())],cadences=int(len(per_cid)))
        report['arms'].append(row)
    report['provenance']={'time_utc':datetime.now(timezone.utc).isoformat(),'protocol_sha256':digest(HERE/'protocol.json'),'code_sha256':{str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__), ROOT/'lib/flux_engine.py',ROOT/'lib/flux_fit.py']},'input_sha256':{str(p.relative_to(ROOT)):digest(p) for p in list(BUNDLE.glob('objective_*.npz'))+[BUNDLE/'exact_prior_setup.json']+list((ROOT/'data/des/model').glob('*'))},'event_ledger_sha256':digest(WORK/'flux-event-ledger.csv')}
    (OUT/'flux-recovery.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':main()
