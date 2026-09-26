"""Declared one-object A/B fixed-mask, fixed-C native profile experiment."""
from pathlib import Path
import argparse,json,time,hashlib
import numpy as np
from scipy.optimize import minimize_scalar
from astropy.io import fits
import paired_refit as p
import native_profile_gate_v2 as ng
from native_mean_oracle import Oracle,OUT
from check_results import table
G=OUT/'global-profile';CID=ng.CID

def prepare():
    G.mkdir(exist_ok=False)
    grid=fits.getdata(p.REL/'model/snoopy.B18/SNooPy_B18.fits','LUMI-GRID').field(0).astype(float)
    a=next(r for r in table(p.OUT/'fits/full/fixed/A/fit.FITRES.TEXT') if r['CID']==CID)
    apar=np.array([float(a[k]) for k in ['DLMAG','STRETCH','AV','PKMJD']],float)
    protocol=dict(status='Frozen before any global profile outcomes',CID=CID,arms=['A author positive','B author signed'],
      conditioning='Fixed released header peak, redshift, calibration, SNooPy model and host RV1.518; modern exactMWlaw99. Empirical AV is allowed negative.',
      primary_covariance='Exact B reference/native.npz from native-profile-gate-v2. For A invert the marginal C_B[A,A], never W_B[A,A]. Same covariance rule on shared raw rows.',
      alternate_covariance='Native full-B covariance exported once at pre-existing original fixed-A/start1 FITRES coordinates, explicitly float32 seed. Retain all B rows; do not refit/tune covariance.',
      alternate_anchor_parameters=apar.astype(np.float32).astype(float).tolist(),
      objective='For each fixed shape/AV, h is exact native mean at fixed Dref. Profile positive amplitude analytically using frozen C; data Gaussian quadratic only, no changing probe prior, no fitted peak penalty. This is an underlying-Gaussian metric, not the positive-selection-conditioned likelihood.',
      domain=dict(shape_knots=grid.tolist(),AV=[-1.,2.],DLMAG_verified_flat_interval=[10.,60.]),
      grid_coarse=dict(shape_subdivisions_each_cell=2,AV_points=33),grid_fine=dict(shape_subdivisions_each_cell=4,AV_points=65),
      continuous_AV='At every shape grid point, preserve endpoints and all sampled local-minimum AV brackets, polish each with bounded scalar minimization xatol1e-7; no shift-sign based branch choice.',
      continuous_shape='For every interior local minimum of each fine shape-profile curve, bounded polish in its neighboring sampled shape interval; each trial runs the full65-point AV scan plus all detected local AV brackets. Cell knots/endpoints remain explicit candidates.',
      gates=dict(coarse_fine_minimum_Q_absolute=1e-3,coarse_fine_DLMAG_absolute=1e-3,amplitude_normal_equation_scaled=1e-10,
        AV_edge='Report if either AV endpoint is within9 of its arm minimum; do not expand domain automatically.',shape_edge='Report any endpoint minimum or endpoint within1 of minimum.',
        amplitude='If b<=0 or profiled D outside10–60, flag unsupported rather than clip.',oracle='Previously passed standalone and randomized replay; repeat fixed baseline after grid, require exact equality.',
        resource='At most100000 unique native mean vectors and20min elapsed; preserve partial arrays and label numerical gate incomplete if exceeded.'),
      reporting='Retain every branch, fixed-C anchor sensitivity, coarse/fine change, boundary and finite-grid limitation. No posterior, confidence level, Bayes factor, physical dust shift or cosmology correction. Delta B−A is conditional fitted-distance sensitivity only after numerical gates.',
      hashes={str(f.relative_to(p.ROOT)):p.sha(f) for f in [Path(__file__),p.OUT/'native_mean_oracle.py',OUT/'build-manifest.json',OUT/'oracle-gate/manifest.json',p.OUT/'native-profile-gate-v2/reference/native.npz',p.OUT/'native-profile-gate-v2/fixed/native.npz',p.OUT/'fits/full/fixed/A/fit.FITRES.TEXT',p.OUT/'data/RSR_A'/f'{CID}.snana.dat',p.OUT/'data/RSR_B'/f'{CID}.snana.dat']})
    p.save(G/'protocol.json',protocol);print(p.sha(G/'protocol.json'))

def run():
    protocol=json.loads((G/'protocol.json').read_text())
    for f,h in protocol['hashes'].items():assert p.sha(p.ROOT/f)==h,f
    ref=np.load(p.OUT/'native-profile-gate-v2/reference/native.npz');nom=np.load(p.OUT/'native-profile-gate-v2/fixed/native.npz');pars=nom['parameters']
    ng.config(G/'A_acceptance','RSR_A',pars);A=ng.run(G/'A_acceptance')
    ng.config(G/'A_covariance_on_B','RSR_B',np.array(protocol['alternate_anchor_parameters']));alt=ng.run(G/'A_covariance_on_B')
    keys=lambda o:[(str(b),float(t),float(y),float(e)) for b,t,y,e in zip(o['band'],o['MJD'],o['data_flux'],o['data_fluxerr'])]
    kb=keys(ref);ka=keys(A);assert len(set(kb))==len(kb) and len(set(ka))==len(ka)
    ai=np.array([kb.index(k) for k in ka]);assert len(ai)==65 and len(kb)==74
    for k in ['band','MJD','data_flux','data_fluxerr']:assert np.array_equal(alt[k],ref[k])
    metrics=[]
    for anchor,C in [('B',ref['C']),('A',alt['C'])]:
      for arm,ix in [('A',ai),('B',np.arange(len(kb)))]:
        sub=C[np.ix_(ix,ix)];L=np.linalg.cholesky((sub+sub.T)/2);T=np.linalg.solve(L,np.eye(len(ix)));y=T@ref['data_flux'][ix]
        metrics.append(dict(label=anchor+'anchor_'+arm,indices=ix,T=T,y=y))
    oracle=Oracle(G/'stream');assert oracle.n==len(kb)
    cache={};history=[];vectors=[];started=time.monotonic();branches=[];violations=[]
    class Limit(Exception):pass
    def evaluate(s,av):
      key=(float(s),float(av))
      if key not in cache:
        if len(cache)>=100000 or time.monotonic()-started>1200:raise Limit()
        pp=pars.copy();pp[1:3]=key;h=oracle.mean(pp);vec=[]
        for m in metrics:
          hw=m['T']@h[m['indices']];b=float(hw@m['y']);q=float(hw@hw)
          if b<=0:
            vec.append((float('inf'),float('nan'),b/q));violations.append(dict(shape=s,AV=av,metric=m['label'],reason='nonpositive amplitude'));continue
          amp=b/q;D=float(pars[0]-2.5*np.log10(amp));r=m['y']-amp*hw
          value=float(r@r)
          assert abs(float(hw@r))<1e-10*max(1,q)
          if not(10<D<60):violations.append(dict(shape=s,AV=av,metric=m['label'],reason='D outside verified flat support'));value=float('inf')
          vec.append((value,D,amp))
        cache[key]=np.array(vec);history.append(key);vectors.append(h)
      return cache[key]
    def prof(s,avgrid,metric,stage):
      sampled=np.array([evaluate(s,av)[metric,0] for av in avgrid]);candidates=[(float(sampled[0]),float(avgrid[0])),(float(sampled[-1]),float(avgrid[-1]))]
      for j in range(1,len(avgrid)-1):
        if sampled[j]<=sampled[j-1] and sampled[j]<=sampled[j+1]:
          opt=minimize_scalar(lambda av:evaluate(s,av)[metric,0],bounds=(avgrid[j-1],avgrid[j+1]),method='bounded',options={'xatol':1e-7,'maxiter':60})
          assert opt.success;candidates.append((float(opt.fun),float(opt.x)))
      # All discrete candidates remain available if a bracket degenerates.
      jj=int(np.argmin(sampled));candidates.append((float(sampled[jj]),float(avgrid[jj])))
      for val,av in candidates:
        vv=evaluate(s,av)[metric];branches.append(dict(stage=stage,metric=metrics[metric]['label'],shape=float(s),AV=av,Q=val,DLMAG=float(vv[1]),amplitude=float(vv[2])))
      q,av=min(candidates);vv=evaluate(s,av)[metric]
      return dict(stage=stage,metric=metrics[metric]['label'],shape=float(s),AV=av,Q=float(vv[0]),DLMAG=float(vv[1]),amplitude=float(vv[2]),AV_low_Q=float(sampled[0]),AV_high_Q=float(sampled[-1]))
    grids=[];shapegrid=np.array(protocol['domain']['shape_knots']);status='complete';failure=None
    try:
      baseline=oracle.mean(pars)
      for stage,subdiv,nav in [('coarse',2,33),('fine',4,65)]:
        sg=np.unique(np.concatenate([np.linspace(a,b,subdiv+1) for a,b in zip(shapegrid[:-1],shapegrid[1:])]))
        avgrid=np.linspace(-1.,2.,nav)
        for i,s in enumerate(sg):
          for k in range(4):grids.append(prof(s,avgrid,k,stage))
          if i%50==0:print(json.dumps(dict(stage=stage,shape_index=i,shape_count=len(sg),native_vectors=len(cache),seconds=time.monotonic()-started)),flush=True)
      for k in range(4):
        rr=[r for r in grids if r['stage']=='fine' and r['metric']==metrics[k]['label']]
        for j in range(1,len(rr)-1):
          if rr[j]['Q']<=rr[j-1]['Q'] and rr[j]['Q']<=rr[j+1]['Q']:
            objective=lambda s:prof(s,np.linspace(-1.,2.,65),k,'shape_polish')['Q']
            opt=minimize_scalar(objective,bounds=(rr[j-1]['shape'],rr[j+1]['shape']),method='bounded',options={'xatol':1e-7,'maxiter':45})
            assert opt.success;grids.append(prof(opt.x,np.linspace(-1.,2.,65),k,'shape_polish_final'))
      assert np.array_equal(oracle.mean(pars),baseline)
    except Limit:
      status='resource limit; incomplete numerical gate'
    except Exception as exc:
      status='exception; incomplete numerical gate';failure=repr(exc)
    finally:
      final=oracle.close()
      for k in ['band','MJD','data_flux','data_fluxerr']:assert np.array_equal(final[k],nom[k])
      assert np.array_equal(final['C'],nom['C'])
      np.savez_compressed(G/'native-profiles.npz',coordinates=np.array(history),model_means=np.array(vectors),metric_results=np.array([cache[k] for k in history]),A_indices=ai,
        covariance_B=ref['C'],covariance_A_anchor=alt['C'],reference_parameters=pars,data_flux=ref['data_flux'],data_fluxerr=ref['data_fluxerr'],band=ref['band'],MJD=ref['MJD'])
      if grids:p.csvsave(G/'shape-profiles.csv',grids)
      if branches:p.csvsave(G/'all-AV-branches.csv',branches)
    summary={}
    for m in metrics:
      label=m['label'];rr=[r for r in grids if r['metric']==label];co=[r for r in rr if r['stage']=='coarse'];fi=[r for r in rr if r['stage']=='fine'];po=[r for r in rr if r['stage']=='shape_polish_final']
      if not co or not fi:continue
      bc=min(co,key=lambda r:r['Q']);bf=min(fi+po,key=lambda r:r['Q']);dq=bc['Q']-bf['Q'];dd=bc['DLMAG']-bf['DLMAG']
      summary[label]=dict(coarse=bc,finest=bf,coarse_minus_finest_Q=dq,coarse_minus_finest_DLMAG=dd,
        coarse_fine_tolerance_pass=abs(dq)<1e-3 and abs(dd)<1e-3,
        AV_edge_within9=min(r[k] for r in fi for k in ['AV_low_Q','AV_high_Q'])-bf['Q']<=9,
        shape_edge_within1=min(fi[0]['Q'],fi[-1]['Q'])-bf['Q']<=1)
    pair={anchor:summary[anchor+'anchor_B']['finest']['DLMAG']-summary[anchor+'anchor_A']['finest']['DLMAG'] for anchor in ['A','B'] if anchor+'anchor_B' in summary and anchor+'anchor_A' in summary}
    result=dict(status=status,exception=failure,native_vectors=len(cache),seconds=time.monotonic()-started,rows_A=len(ai),rows_B=len(kb),metrics=summary,
      delta_DLMAG_B_minus_A=pair,unsupported_amplitude_or_distance=violations,protocol_sha256=p.sha(G/'protocol.json'),
      interpretation='Conditional fixed-C profile geometry only; no selected-data likelihood, physical bias correction, posterior, or cosmology result.')
    p.save(G/'result.json',result);p.save(G/'manifest.json',dict(files_sha256={str(f.relative_to(G)):p.sha(f) for f in G.rglob('*') if f.is_file() and f.name!='manifest.json'}))
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['prepare','run']);a=ap.parse_args();prepare() if a.action=='prepare' else run()
