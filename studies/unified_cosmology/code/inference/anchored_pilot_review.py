"""Independent synthetic exercise of the actual pilot driver; no physical calls."""
from contextlib import ExitStack
import hashlib
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from scipy.special import logsumexp
from scipy.stats import multivariate_normal, multivariate_t

ROOT=Path(__file__).resolve().parents[4]
HERE=ROOT/'studies/unified_cosmology/code/inference'
sys.path.insert(0,str(HERE))
import anchored_pilot as p
import camb
import cobaya.model
import independence_runtime
from independence_proposal import FrozenMixture

def forbidden(*args,**kwargs):raise AssertionError('No physical theory/background/model call permitted.')

mixture=FrozenMixture(np.array([.3,-.4]),np.array([[.8,.17],[.17,1.2]]),['x','y'])
design=json.loads(p.DESIGN.read_text())
request=p.requests(mixture,'lcdm',design)
independent=logsumexp(np.array([
    np.log(.9)+multivariate_normal.logpdf(request['points'],mixture.mean,1.05**2*mixture.covariance),
    np.log(.1)+multivariate_t.logpdf(request['points'],mixture.mean,4*3/5*mixture.covariance,df=5)]),axis=0)
error=float(np.max(abs(independent-request['logq'])));assert error<1e-12
cases=[]
with patch.object(camb,'get_results',forbidden),patch.object(camb,'get_background',forbidden),patch.object(camb,'get_transfer_functions',forbidden):
  for scenario in ['all400','four_fallbacks','wall_after_completed_call','exception_is_not_rejection','nan_is_not_rejection']:
    with tempfile.TemporaryDirectory(dir=ROOT/'.work') as td:
      parent=Path(td);out=parent/'pilot'
      valid=parent/'validation.json'
      valid.write_text(json.dumps({'status':'passed_synthetic_pilot_density_and_retention_checks','dependency_sha256':p.dependencies()}))
      clock=SimpleNamespace(now=0.)
      clock.monotonic=lambda:clock.now
      theory=SimpleNamespace(exact_calls=0,surrogate_calls=0)
      truth=[]
      class FakeModel:
        def __init__(self):
          self.theory={'spectral_surrogate':theory}
          self.likelihood=['synthetic_gaussian']
          self.parameterization=SimpleNamespace(sampled_params=lambda:['x','y'])
        def __enter__(self):return self
        def __exit__(self,*args):return False
        def logposterior(self,point):
          i=len(truth);truth.append(dict(point));clock.now+=201 if scenario=='wall_after_completed_call' else .125
          theory.surrogate_calls+=1
          if scenario=='four_fallbacks':theory.exact_calls+=1
          if scenario=='exception_is_not_rejection' and i==4:raise ValueError('Synthetic unexpected solver failure')
          if scenario=='nan_is_not_rejection' and i==4:return SimpleNamespace(logpost=np.nan,logprior=0.,loglikes=[np.nan])
          if i<2:lp=-np.inf
          else:lp=float(multivariate_normal.logpdf([point['x'],point['y']],[-.2,.1],[[1.,-.1],[-.1,.7]]))
          return SimpleNamespace(logpost=lp,logprior=0.,loglikes=[lp])
      info={'params':{'x':{'prior':{}},'y':{'prior':{}}}}
      synthetic_identity={'identity':'explicit-synthetic-consumer-fixture'}
      provenance={'proposal_kind':'anchored_bridge_weighted_v1','scope':'Synthetic driver fixture, not actual training.'}
      with ExitStack() as stack:
        for obj,name,value in [(p,'VALIDATION',valid),(p,'time',clock)]:stack.enter_context(patch.object(obj,name,value))
        stack.enter_context(patch.object(p.s,'runtime',return_value={'synthetic':True}))
        stack.enter_context(patch.object(p.s,'validation_guard',return_value={'synthetic':True}))
        stack.enter_context(patch.object(p.s,'configuration',return_value=info))
        stack.enter_context(patch.object(p.s,'identify',return_value=synthetic_identity))
        stack.enter_context(patch.object(p.s,'load_proposal',return_value=(mixture,provenance)))
        stack.enter_context(patch.object(p.os,'sched_getaffinity',return_value={0}))
        stack.enter_context(patch.object(p.os,'sched_setaffinity',return_value=None))
        warm=stack.enter_context(patch.object(independence_runtime,'native_initialize',return_value={'synthetic':True}))
        factory=stack.enter_context(patch.object(cobaya.model,'get_model',side_effect=lambda info:FakeModel()))
        result=p.pilot({'model':'lcdm'},parent/'proposal',out,0)
        assert warm.call_count==factory.call_count==1
      rows=[json.loads(v)for v in (out/'records.jsonl').read_text().splitlines()]
      assert len(truth)==len(rows)==result['target_calls_recorded']
      current=None;current_lp=None;current_lq=None;accepted=0;occupied=[]
      for row in rows:
        i=row['index'];assert row['point']==dict(zip(['x','y'],request['points'][i]))
        if row.get('status')=='failed_evaluation':
          assert row['occupied_candidate'] is None
          assert scenario in ['exception_is_not_rejection','nan_is_not_rejection'];continue
        lp=row['target_logpost']
        if lp is None:expected=False
        elif current is None:expected=True
        else:expected=request['logu'][i]<min(0.,lp-current_lp+current_lq-independent[i])
        assert expected==row['accepted']
        if expected:current=i;current_lp=lp;current_lq=independent[i];accepted+=1
        assert row['occupied_candidate']==current
        if current is not None:occupied.append(current)
      if scenario=='all400':assert len(rows)==400 and result['status']=='completed_declared_request_budget'
      if scenario=='four_fallbacks':assert len(rows)==4 and result['successful_native_fallbacks']==4 and result['status']=='stopped_native_fallback_budget'
      if scenario=='wall_after_completed_call':assert len(rows)==3 and result['seconds']==603 and result['status']=='stopped_wall_budget_after_completed_call'
      if scenario in ['exception_is_not_rejection','nan_is_not_rejection']:assert len(rows)==5 and result['status']=='failed_evaluation'
      assert result['accepted_transitions_after_initialization']==accepted-1
      assert result['occupied_states_recorded_including_rejections']==len(occupied)
      assert result['qualified_scientific_measurement'] is False
      assert not {'posterior','posterior_mean','posterior_covariance','qualified_under_declared_numerical_gates'}&set(result)
      assert result['explicit_native_initialization_calls']==1
      cases.append({'scenario':scenario,'target_calls':len(rows),'status':result['status'],
                    'retained_occupied_states':len(occupied),'preinitialization_rejections':2,
                    'acceptance_and_repetition_independently_rebuilt':True,'no_measurement_fields':True})
paths=[Path(__file__)]+[HERE/n for n in ['anchored_pilot.py','anchored-pilot-design.json','anchored_pilot_validate.py',
    'anchored_sampling.py','anchored_correction.py','anchored_measurement.py','anchored_proposal.py','anchored-proposal-design.json','anchored-sampling-design.json']]
output={'status':'passed_independent_synthetic_pilot_driver_controls','cases':cases,
        'normalized_mixture_density_error':error,'real_model_or_background_or_native_calls':0,
        'synthetic_monkeypatch_scope':'Only isolated fake identities, prerequisites, model, initialization and affinity for driver-control tests. No scientific manifest was modified or qualified.',
        'source_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()for path in paths}}
(ROOT/'studies/unified_cosmology/results/inference/anchored-pilot-review-controls.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps(output,indent=2))
