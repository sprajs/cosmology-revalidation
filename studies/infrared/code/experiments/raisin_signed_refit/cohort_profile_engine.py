"""Declared one-object A/B fixed-mask, fixed-C native profile experiment."""
from pathlib import Path
import argparse,json,time,hashlib,re,os
from collections import defaultdict,deque
import numpy as np
from scipy.optimize import minimize_scalar
from astropy.io import fits
import paired_refit as p
import native_profile_gate_v2 as ng
from native_mean_oracle import Oracle,OUT
from check_results import table
G=Path(os.environ['RAISIN_PROFILE_CASE']);CID=os.environ['RAISIN_PROFILE_CID'];ng.CID=CID
original_config=ng.config
def short_output_config(*args,**kwargs):
    nml=original_config(*args,**kwargs)
    text=nml.read_text();text=re.sub(r'(?m)^\s*TEXTFILE_PREFIX\s*=.*$'," TEXTFILE_PREFIX = 'fit'",text)
    nml.write_text(text);return nml
ng.config=short_output_config

def run():
    protocol=json.loads((G/'protocol.json').read_text())
    assert p.sha(Path(__file__))==protocol['engine_sha256']
    for f,h in protocol['hashes'].items():assert p.sha(p.ROOT/f)==h,f
    ng.config(G/'reference','RSR_B',None);ref=ng.run(G/'reference')
    original=next(r for r in table(p.OUT/'fits/full/fixed/B/fit.FITRES.TEXT') if r['CID']==CID)
    rerun=table(G/'reference/fit.FITRES.TEXT')[0]
    deltas={k:float(rerun[k])-float(original[k]) for k in ['DLMAG','STRETCH','AV','PKMJD','FITCHI2','NDOF']}
    assert np.allclose([deltas[k] for k in ['DLMAG','STRETCH','AV','PKMJD']],0,rtol=0,atol=[1e-5,1e-5,1e-5,1e-4]),deltas
    assert abs(deltas['FITCHI2'])<1e-6*max(1,float(original['FITCHI2'])) and deltas['NDOF']==0
    ng.config(G/'fixed_reference','RSR_B',ref['parameters'].astype(np.float32).astype(float));nom=ng.run(G/'fixed_reference')
    p.save(G/'reference-closure.json',dict(status='PASS',original_current_minus_audit=deltas,reference_epochs=len(ref['MJD'])))
    ref=np.load(G/'reference/native.npz');nom=np.load(G/'fixed_reference/native.npz');pars=nom['parameters']
    ng.config(G/'A_acceptance','RSR_A',pars);A=ng.run(G/'A_acceptance')
    ng.config(G/'A_covariance_on_B','RSR_B',np.array(protocol['alternate_anchor_parameters']));alt=ng.run(G/'A_covariance_on_B')
    keys=lambda o:[(str(b),float(t),float(y),float(e)) for b,t,y,e in zip(o['band'],o['MJD'],o['data_flux'],o['data_fluxerr'])]
    kb=keys(ref);ka=keys(A);queues=defaultdict(deque);groups=defaultdict(list)
    for j,key in enumerate(kb):queues[key].append(j);groups[key].append(j)
    ai=np.array([queues[key].popleft() for key in ka]);assert len(ai)==protocol['expected_A_epochs'] and len(kb)==protocol['expected_B_epochs']
    duplicate_checks=[]
    for key,indices in groups.items():
      if len(indices)<2:continue
      for index in indices[1:]:
        permutation=np.arange(len(kb));permutation[indices[0]],permutation[index]=permutation[index],permutation[indices[0]]
        for anchor,C in [('B',ref['C']),('A',alt['C'])]:
          error=float(np.max(abs(C[np.ix_(permutation,permutation)]-C))/np.max(abs(C)))
          duplicate_checks.append(dict(indices=[indices[0],index],anchor=anchor,relative_swap_error=error));assert error<1e-12
    p.save(G/'mask-mapping.json',dict(positive_A_epochs=len(ai),signed_B_epochs=len(kb),A_indices=ai.tolist(),duplicate_swap_checks=duplicate_checks,
      original_row_multiplicity_retained=True,shared_raw_float32_values_exact=True))
    for k in ['band','MJD','data_flux','data_fluxerr']:assert np.array_equal(alt[k],ref[k])
    metrics=[];asymmetry={}
    for anchor,C in [('B',ref['C']),('A',alt['C'])]:
      asymmetry[anchor]=float(np.max(abs(C-C.T))/np.max(abs(C)));assert asymmetry[anchor]<1e-12
      for arm,ix in [('A',ai),('B',np.arange(len(kb)))]:
        sub=C[np.ix_(ix,ix)];L=np.linalg.cholesky((sub+sub.T)/2);T=np.linalg.solve(L,np.eye(len(ix)));y=T@ref['data_flux'][ix]
        metrics.append(dict(label=anchor+'anchor_'+arm,indices=ix,T=T,y=y))
    oracle=Oracle(G/'stream');assert oracle.n==len(kb)
    cache={};history=[];vectors=[];started=time.monotonic();branches=[];violations=[]
    class Limit(Exception):pass
    def evaluate(s,av):
      key=(float(s),float(av))
      if key not in cache:
        if oracle.counter>=100000 or time.monotonic()-started>1200:raise Limit()
        pp=pars.copy();pp[1:3]=key;h=oracle.mean(pp);vec=[]
        badmean=bool(np.any(h<=0) or np.any(h<=1e-5) or np.any(h>=1e9))
        if badmean:violations.append(dict(shape=s,AV=av,reason='native nonpositive/clipped mean domain',rows=np.flatnonzero((h<=1e-5)|(h>=1e9)).tolist()))
        for m in metrics:
          hw=m['T']@h[m['indices']];b=float(hw@m['y']);q=float(hw@hw)
          if q<=0 or b<=0 or badmean:
            vec.append((float('inf'),float('nan'),float('nan')));violations.append(dict(shape=s,AV=av,metric=m['label'],reason='nonpositive amplitude or information or clipped mean'));continue
          amp=b/q;D=float(pars[0]-2.5*np.log10(amp));r=m['y']-amp*hw
          value=float(r@r)
          assert abs(float(hw@r))<1e-10*max(1,q)
          if not(10<D<60) or np.any(amp*h<=1e-5) or np.any(amp*h>=1e9):violations.append(dict(shape=s,AV=av,metric=m['label'],reason='D or scaled mean outside verified support'));value=float('inf')
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
    grids=[];shapegrid=np.array(protocol['domain']['shape_knots']);status='complete';failure=None;boundary_checks=[];amplitude_checks=[]
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
      native_step=float(np.float32(np.float32(shapegrid[-1]-shapegrid[0])/np.float32(len(shapegrid)-1)))
      algorithm=shapegrid[0]+np.arange(1,len(shapegrid)-1)*native_step
      for k in range(4):
        rr=[r for r in grids if r['stage']=='fine' and r['metric']==metrics[k]['label']]
        for boundary in algorithm:
          near=min(rr,key=lambda r:abs(r['shape']-boundary))
          for offset in [-1e-8,0.,1e-8]:
            s=float(boundary+offset);av=near['AV'];vv=evaluate(s,av)[k]
            candidate=dict(stage='algorithm_boundary_probe',metric=metrics[k]['label'],shape=s,AV=av,Q=float(vv[0]),DLMAG=float(vv[1]),amplitude=float(vv[2]))
            boundary_checks.append(candidate);grids.append(candidate)
        best=min([r for r in grids if r['metric']==metrics[k]['label']],key=lambda r:r['Q'])
        if best['stage']=='algorithm_boundary_probe':grids.append(prof(best['shape'],np.linspace(-1.,2.,65),k,'boundary_polish_final'))
      for k in range(4):
        rr=[r for r in grids if r['metric']==metrics[k]['label']];minimum=min(r['Q'] for r in rr)
        for r in rr:
          if r['Q']>minimum+9:continue
          if oracle.counter>=99998:raise Limit()
          pp=pars.copy();pp[1:3]=[r['shape'],r['AV']];h=oracle.mean(pp);pp[0]=r['DLMAG'];actual=oracle.mean(pp)
          expected=r['amplitude']*h;error=float(np.max(abs(actual-expected))/max(np.max(abs(expected)),1e-30))
          amplitude_checks.append(dict(metric=r['metric'],shape=r['shape'],AV=r['AV'],DLMAG=r['DLMAG'],relative_error=error))
          assert error<2e-8,('competitive native amplitude domain',r,error)
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
      if boundary_checks:p.csvsave(G/'algorithm-boundary-probes.csv',boundary_checks)
      if amplitude_checks:p.csvsave(G/'competitive-amplitude-checks.csv',amplitude_checks)
    summary={}
    for m in metrics:
      label=m['label'];rr=[r for r in grids if r['metric']==label];co=[r for r in rr if r['stage']=='coarse'];fi=[r for r in rr if r['stage']=='fine'];po=[r for r in rr if r['stage'] in ['shape_polish_final','algorithm_boundary_probe','boundary_polish_final']]
      if not co or not fi:continue
      bc=min(co,key=lambda r:r['Q']);bf=min(fi+po,key=lambda r:r['Q']);dq=bc['Q']-bf['Q'];dd=bc['DLMAG']-bf['DLMAG']
      summary[label]=dict(coarse=bc,finest=bf,coarse_minus_finest_Q=dq,coarse_minus_finest_DLMAG=dd,
        coarse_fine_tolerance_pass=abs(dq)<1e-3 and abs(dd)<1e-3,
        AV_edge_within9=min(r[k] for r in fi for k in ['AV_low_Q','AV_high_Q'])-bf['Q']<=9,
        shape_edge_within1=min(fi[0]['Q'],fi[-1]['Q'])-bf['Q']<=1)
    pair={anchor:summary[anchor+'anchor_B']['finest']['DLMAG']-summary[anchor+'anchor_A']['finest']['DLMAG'] for anchor in ['A','B'] if anchor+'anchor_B' in summary and anchor+'anchor_A' in summary}
    computational_completion=status=='complete'
    numerical_gate_pass=computational_completion and not violations and len(summary)==4 and all(r['coarse_fine_tolerance_pass'] and not r['AV_edge_within9'] and not r['shape_edge_within1'] for r in summary.values())
    if computational_completion and not numerical_gate_pass:status='computationally complete; numerical/support gate not passed'
    result=dict(status=status,computational_completion=computational_completion,numerical_gate_pass=numerical_gate_pass,exception=failure,native_vectors=len(cache),total_oracle_calls=oracle.counter,seconds=time.monotonic()-started,rows_A=len(ai),rows_B=len(kb),metrics=summary,covariance_relative_asymmetry=asymmetry,
      competitive_amplitude_checks=len(amplitude_checks),max_competitive_amplitude_relative_error=max([r['relative_error'] for r in amplitude_checks],default=None),
      delta_DLMAG_B_minus_A=pair,unsupported_amplitude_or_distance=violations,protocol_sha256=p.sha(G/'protocol.json'),
      cohort_protocol_sha256=protocol['cohort_protocol_sha256'],reference_gate=deltas,duplicate_swap_checks=duplicate_checks,
      interpretation='Conditional fixed-C profile geometry only; no selected-data likelihood, physical bias correction, posterior, or cosmology result.')
    p.save(G/'result.json',result);p.save(G/'manifest.json',dict(files_sha256={str(f.relative_to(G)):p.sha(f) for f in G.rglob('*') if f.is_file() and f.name!='manifest.json'}))
    print(json.dumps(result,indent=2))

if __name__=='__main__':run()
