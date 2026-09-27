"""Independent normalized-density/MH/retained-state controls; no physical calls."""
import argparse
import json
from pathlib import Path
from unittest.mock import patch
import numpy as np
from scipy.stats import multivariate_normal,multivariate_t

import anchored_pilot as p
from independence_proposal import FrozenMixture


def validate():
    rng=np.random.default_rng(8274421)
    q=FrozenMixture(np.array([.2,-.1]),np.array([[1.,.2],[.2,.8]]),['x','y'])
    design=json.loads(p.DESIGN.read_text())
    first=p.requests(q,'lcdm',design);second=p.requests(q,'lcdm',design)
    assert all(np.array_equal(first[k],second[k]) for k in first)
    assert first['points'].shape==(400,2)
    independent=np.logaddexp(np.log(.9)+multivariate_normal.logpdf(first['points'],q.mean,q.covariance*1.05**2),
                             np.log(.1)+multivariate_t.logpdf(first['points'],q.mean,q.covariance*4*3/5,df=5))
    density_error=float(np.max(abs(independent-first['logq'])));assert density_error<1e-12
    x=rng.normal(size=(10000,2));y=rng.normal(size=(10000,2))
    lp_x=multivariate_normal.logpdf(x,np.zeros(2),np.eye(2));lp_y=multivariate_normal.logpdf(y,np.zeros(2),np.eye(2))
    lq_x=q.logpdf(x);lq_y=q.logpdf(y)
    balance=[]
    for a,b,qa,qb in zip(lp_x,lp_y,lq_x,lq_y):
        _,forward=p.transition(b,a,qa,qb,-1.)
        _,reverse=p.transition(a,b,qb,qa,-1.)
        balance.append(abs(a+qb+forward-b-qa-reverse))
    assert max(balance)<1e-12
    # A bounded independent known-Gaussian target verifies normalized draws,
    # the MH direction and the retention of all rejected states together.
    n=60000;mu=np.array([-.15,.25]);tcov=np.array([[.75,-.12],[-.12,1.1]])
    current=None;current_lp=None;current_lq=None;occupied=[];accepted=0
    for i in range(n):
        candidate,_=q.draw(rng);lp=float(multivariate_normal.logpdf(candidate,mu,tcov));lq=q.logpdf(candidate)
        yes,_=p.transition(lp,current_lp,current_lq,lq,float(np.log(rng.random())))
        if yes:current=candidate;current_lp=lp;current_lq=lq;accepted+=1
        occupied.append(current.copy())
    sample=np.array(occupied)[2000:]
    mean_error=float(np.max(abs(sample.mean(0)-mu)))
    covariance_error=float(np.max(abs(np.cov(sample,rowvar=False)-tcov)))
    assert mean_error<.025 and covariance_error<.035
    # Explicit support rejection before and after initialization never redraws.
    sequence=[-np.inf,-1.,-np.inf,-1.2,-np.inf]
    cur=None;clp=None;clq=None;ids=[];rows=[]
    for i,lp in enumerate(sequence):
        yes,alpha=p.transition(lp,clp,clq,-2.,-10.)
        if yes:cur=i;clp=lp;clq=-2.
        ids.append(cur)
        rows.append({'index':i,'occupied_candidate':cur,'finite_target':bool(np.isfinite(lp)),
                     'target_logpost':lp if np.isfinite(lp) else None,'logq':-2.,'accepted':yes})
    assert ids==[None,1,1,3,3]
    metric=p.metrics(rows,{'points':np.arange(10).reshape(5,2),'names':['x','y']})
    assert metric['occupied_states_recorded_including_rejections']==4
    assert metric['longest_hold']==2 and metric['accepted_transitions_after_initialization']==1
    assert metric['qualified_scientific_measurement'] is False
    json.dumps(metric,allow_nan=False)
    for lp in [np.nan,np.inf]:
        try:p.transition(lp,-1.,-1.,-1.,-1.)
        except AssertionError:pass
        else:raise AssertionError('Undefined density accepted.')
    # A missing prerequisite must fail before constructing a model or touching
    # an output directory. This is deliberately a synthetic refusal fixture.
    with patch.object(p.s,'validation_guard',side_effect=ValueError('unqualified-prerequisite')):
        try:p.pilot({},Path('absent'),Path('absent'),0)
        except ValueError as error:assert str(error)=='unqualified-prerequisite'
        else:raise AssertionError('Prerequisite bypassed.')
    return {'status':'passed_synthetic_pilot_density_and_retention_checks',
        'normalized_mixture_logpdf_error':density_error,'detailed_balance_max_error':max(balance),
        'known_Gaussian_MH_draws':n,'known_Gaussian_mean_error':mean_error,
        'known_Gaussian_covariance_error':covariance_error,'retained_support_rejection_ids':ids,
        'native_or_background_calls':0,'dependency_sha256':p.dependencies(),
        'source_sha256':{p.s.relative(Path(__file__)):p.s.digest(__file__)},
        'limitations':'Only synthetic target arithmetic and prerequisite refusal; no real anchored pilot, sampler speed claim or scientific measurement.'}


def main():
    a=argparse.ArgumentParser();a.add_argument('--output',type=Path,required=True);args=a.parse_args()
    result=validate();p.s.dump_new(args.output,result);print(json.dumps({k:v for k,v in result.items() if k not in ['dependency_sha256','source_sha256']}))

if __name__=='__main__':main()
