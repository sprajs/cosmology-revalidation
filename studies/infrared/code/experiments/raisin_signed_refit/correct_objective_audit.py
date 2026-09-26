"""Preserve the first audit; distinguish native data chi2 from its minimand."""
import json,re
from pathlib import Path
import numpy as np
import paired_refit as p
from check_results import table,gate
from check_stability import mask
O=p.OUT
D=O/'qualified-stability'

def priors(d):
    text=(d/'fit.log').read_text()
    return {cid:float(v) for cid,v in re.findall(r'FIT MINIMUM \(([^)]+)\) returns PRIOR\s+CHI2/dof=\s*([^\s/]+)\s*/',text)}

def main():
    D.mkdir(exist_ok=True)
    old=['stability-result.json','multistart-comparison.csv','multistart-fit-ledger.csv','paired-start-sensitivity.csv','shape-profile.csv','check_stability.py','stability-protocol.json','result.json']
    preservation={str((O/f).relative_to(p.ROOT)):p.sha(O/f) for f in old}
    rows=[];cases=[];inputs={}
    for timing in ['free','fixed']:
      for arm in ['R','A','B','H']:
        ds=[O/'fits/full'/timing/arm]+[O/'fits/multistart'/str(s)/timing/arm for s in [.85,1.15]]
        bystart=[]
        for start,d in zip([1.,.85,1.15],ds):
          pp=priors(d); rr=gate(d,arm); bystart.append({r['CID']:r for r in rr})
          for f in ['fit.log','fit.nml','fit.LCPLOT.TEXT']:
            inputs[str((d/f).relative_to(p.ROOT))]=p.sha(d/f)
          for r in rr:
            assert r['CID'] in pp
            rows.append(dict(start=start,**r,printed_final_prior_chi2=pp[r['CID']],display_precision_minimand=r['FITCHI2']+pp[r['CID']],minimand_interpretation='Data chi2 plus rounded native prior; iteration-specific C and prior center, not a shared likelihood'))
        for cid in bystart[0]:
          rr=[r[cid] for r in bystart];mm=[mask(d,cid) for d in ds]
          cases.append(dict(CID=cid,timing=timing,arm=arm,same_accepted_mask=all(m==mm[0] for m in mm[1:]),same_C_and_prior_not_established=True,common_objective_comparable=False,
            data_chi2_range=float(np.ptp([r['FITCHI2'] for r in rr])),**{k+'_range':float(np.ptp([r[k] for r in rr])) for k in ['DLMAG','AV','STRETCH','PKMJD']}))
    p.csvsave(D/'fit-ledger.csv',rows);p.csvsave(D/'case-ledger.csv',cases)
    profiles=[];protocol=json.loads((O/'stability-protocol.json').read_text())
    for j in protocol['jobs']:
      if j['kind']!='profile':continue
      d=(p.ROOT/j['nml']).parent;r=table(d/'fit.FITRES.TEXT')[0];cid=r['CID'];original=O/'fits/full/free/B'
      assert abs(float(r['STRETCH'])-j['fixed_shape'])<2e-7
      pp=priors(d)[cid]
      profiles.append(dict(CID=cid,offset_cells=j['offset_cells'],fixed_shape=j['fixed_shape'],same_mask_as_original=mask(d,cid)==mask(original,cid),accepted=len(mask(d,cid)),
        **{k:float(r[k]) for k in ['STRETCH','STRETCHERR','AV','DLMAG','PKMJD','NDOF','FITCHI2','ERRFLAG_FIT']},printed_final_prior_chi2=pp,display_precision_minimand=float(r['FITCHI2'])+pp,common_objective_comparable=False))
      for f in ['fit.log','fit.nml','fit.LCPLOT.TEXT']:inputs[str((d/f).relative_to(p.ROOT))]=p.sha(d/f)
    assert len(profiles)==20
    p.csvsave(D/'shape-profile.csv',profiles)
    profile_summary={cid:dict(points=len(rr),same_mask_as_original=sum(r['same_mask_as_original'] for r in rr),accepted_counts=sorted({r['accepted'] for r in rr}),DLMAG_range=float(np.ptp([r['DLMAG'] for r in rr])),data_chi2_range=float(np.ptp([r['FITCHI2'] for r in rr]))) for cid in sorted({r['CID'] for r in profiles}) for rr in [[r for r in profiles if r['CID']==cid]]}
    result=dict(status='Numerical readiness FAIL; do not infer a preferred correction, global minimum, or likelihood interval',native_multistart_fits=len(rows),native_profile_fits=len(profiles),
      cases=len(cases),same_mask_cases=sum(c['same_accepted_mask'] for c in cases),common_likelihood_proven_cases=0,
      all_native_ERRFLAG_zero=all(r['ERRFLAG_FIT']==0 for r in rows+profiles),all_multistart_cov_positive=all(r['covariance_min_eigenvalue']>0 for r in rows),
      max_parameter_range={k:max(c[k+'_range'] for c in cases) for k in ['DLMAG','AV','STRETCH','PKMJD']},
      cases_DLMAG_spread_above_0p001=sum(c['DLMAG_range']>.001 for c in cases),cases_DLMAG_spread_above_0p01=sum(c['DLMAG_range']>.01 for c in cases),
      printed_prior_range=[min(r['printed_final_prior_chi2'] for r in rows),max(r['printed_final_prior_chi2'] for r in rows)],profiles=profile_summary,
      corrections=['FITCHI2 excludes CHI2INI; prior added here is rounded log output, not exact FCN dump.',
      'Peak prior center updates to previous iteration fit, not immutable header peak. The header sets initial seed.',
      'Final C uses previous iteration model and MW error, even with OPT_COVAR_FLUX=0. Same mask does not establish a common objective.',
      'The original fields best_start_by_objective/comparable_objective/conditionally_preferred must not be used for scientific ranking.',
      '20 local profile points were executed; the original resource prose said 21 but one out-of-grid endpoint was omitted before execution.'])
    p.save(D/'result.json',result)
    for path,h in preservation.items():assert p.sha(p.ROOT/path)==h
    inputs.update(preservation)
    for f in [Path(__file__),p.BUILD/'src/snlc_fit.F90']:
      inputs[str(f.relative_to(p.ROOT))]=p.sha(f)
    p.save(D/'manifest.json',dict(inputs_sha256=inputs,preserved_legacy_artifacts_sha256=preservation,outputs_sha256={f.name:p.sha(f) for f in D.glob('*') if f.is_file() and f.name!='manifest.json'}))
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
