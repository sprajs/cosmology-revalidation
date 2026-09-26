"""Retrospective object-cross-fit test of existing held-out epoch residuals.

See residual-protocol.md; all interpretation is fixed-design Gaussian/local.
"""
from pathlib import Path
import sys, json, hashlib, platform
import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular
from scipy.stats import beta

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts/phase2/independent_flux'))
from fit import Case, Engine, jacobian

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def spectral_templates(case, rec):
    p=np.array(rec['x']);p[3]+=case.tref
    phase=(case.t-p[3])/(1+case.z)
    result=np.empty((len(phase),3));engine=case.engine
    for b in np.unique(case.b):
        ii=np.flatnonzero(case.b==b);g=case.grid[b]
        shape=(engine.snmodel.source._model['M0'](phase[ii],g['rest'])+
               p[1]*engine.snmodel.source._model['M1'](phase[ii],g['rest']))/1e-12
        flux_elements=np.exp(p[0])*shape*10**(-.4*p[2]*g['cl'])*g['weights']
        features=[np.exp(-.5*((g['rest']-4300)/600)**2),
                  np.exp(-.5*((g['rest']-6000)/600)**2),
                  np.tanh(phase[ii,None]/20)*g['cl']]
        for k,feature in enumerate(features):
            result[ii,k]=-.4*np.log(10)*np.sum(flux_elements*feature,axis=1)
    return result

def prepare():
    file=ROOT/'phase2/independent_flux/heldout_predictions.json'
    records=[r for r in json.loads(file.read_text()) if r['split']=='epoch' and r['family']=='salt']
    assert len(records)==12
    engine=Engine();rs=[];families={'observer':[],'rest_phase':[]};info=[];start=0;inputs=[file]
    for rec in records:
        assert rec['success'] and not rec['boundary'] and rec['prediction']['valid']
        c=Case(engine,rec['CID']);x=np.array(rec['x'])
        mu=c.flux(x);j=jacobian(lambda p:c.flux(p),x)
        tr=np.array(rec['train_rows']);te=np.array(rec['test_rows'])
        Ctr=c.cov[np.ix_(tr,tr)]
        A=np.linalg.solve(Ctr,c.cov[np.ix_(tr,te)]).T
        jp=j[te]-A@j[tr]
        condmean=mu[te]+A@(c.y[tr]-mu[tr])
        h=np.array(rec['hessian_covariance'])
        pcov=c.cov[np.ix_(te,te)]-A@c.cov[np.ix_(tr,te)]+jp@h@jp.T
        assert np.allclose(condmean,rec['prediction']['test_model'],rtol=1e-10,atol=1e-9)
        assert np.allclose(np.diag(pcov),rec['prediction']['test_variance'],rtol=1e-9,atol=1e-8)
        chol=np.linalg.cholesky(pcov)
        r=solve_triangular(chol,c.y[te]-condmean,lower=True)
        rs.append(r)
        prior=np.diag([0.,0.,0.,.01])
        normal=j[tr].T@np.linalg.solve(Ctr,j[tr])+prior
        def adjusted(G):
            rtrain=np.linalg.solve(normal,j[tr].T@np.linalg.solve(Ctr,G[tr]))
            return solve_triangular(chol,G[te]-A@G[tr]-jp@rtrain,lower=True)
        gray=adjusted((-.4*np.log(10)*mu)[:,None])
        assert np.max(abs(gray))<1e-5
        bandcontrast=np.column_stack([(c.b=='g').astype(float)-(c.b=='r'),
                                     (c.b=='i').astype(float)-(c.b=='r'),
                                     (c.b=='z').astype(float)-(c.b=='r')])
        obs=-.4*np.log(10)*mu[:,None]*bandcontrast
        rest=spectral_templates(c,rec)
        for name,G in [('observer',obs),('rest_phase',rest)]:
            families[name].append(adjusted(G))
        info.append({'CID':c.cid,'zHEL':float(c.z),'ntrain':len(tr),'ntest':len(te),
                     'row_start':start,'row_end':start+len(te),
                     'baseline_conditional_chi2':float(r@r),
                     'max_gray_response_after_refit':float(abs(gray).max()),
                     'bands':''.join(sorted(set(c.b))),
                     'condition_number_normal':float(np.linalg.cond(normal))})
        start+=len(te)
        inputs.append(ROOT/f'phase2/official/portable_pilot/objective_{c.cid}.npz')
    return np.concatenate(rs),{k:np.vstack(v) for k,v in families.items()},info,inputs

def score_operator(T,info,sigma):
    n=len(T);Q=np.zeros((n,n));constant=0.;folds=[]
    for row in info:
        te=np.arange(row['row_start'],row['row_end']);tr=np.setdiff1d(np.arange(n),te)
        V=np.linalg.inv(np.eye(T.shape[1])/sigma**2+T[tr].T@T[tr])
        pred=np.eye(len(te))+T[te]@V@T[te].T
        E=np.zeros((len(te),n));E[np.arange(len(te)),te]=1
        error=E.copy();error[:,tr]-=T[te]@V@T[tr].T
        qi=E.T@E-error.T@np.linalg.solve(pred,error)
        ci=-.5*np.linalg.slogdet(pred)[1]
        Q+=qi;constant+=ci
        folds.append((qi,ci,V,te,tr))
    return .5*(Q+Q.T),constant,folds

def tail(vals,threshold):
    n=len(vals);k=int(np.sum(vals>=threshold));return {'exceedances':k,'draws':n,
      'plus_one_tail':(k+1)/(n+1),
      'binomial_95_interval':[0. if k==0 else float(beta.ppf(.025,k,n-k+1)),1. if k==n else float(beta.ppf(.975,k+1,n-k))]}

def main():
    if (OUT/'residual-transfer.json').exists():raise RuntimeError('Preserve completed run')
    r,templates,info,inputs=prepare();n=len(r);seed=26092677;rng=np.random.default_rng(seed)
    noise=rng.normal(size=(10000,n));results={};operators={};mc=[];per=[]
    for family,T in templates.items():
        for sigma in [.01,.02,.05]:
            label=f'{family}_{sigma:g}'
            Q,k,folds=score_operator(T,info,sigma)
            val=float(.5*r@Q@r+k)
            sims=.5*np.einsum('bi,ij,bj->b',noise,Q,noise)+k
            mc.append(sims);operators[label]=(Q,k)
            foldvals=[];maxcoef=0.
            for row,(qf,kf,V,te,tr) in zip(info,folds):
                coef=V@T[tr].T@r[tr]
                maxcoef=max(maxcoef,float(abs(coef).max()))
                score=float(.5*r@qf@r+kf);foldvals.append(score)
                per.append({'candidate':label,'CID':row['CID'],'delta_log_score':score,
                            'coefficient_0':coef[0],'coefficient_1':coef[1],'coefficient_2':coef[2]})
            assert np.isclose(sum(foldvals),val,atol=1e-10)
            results[label]={'sum_delta_log_score':val,'positive_objects':int(np.sum(np.array(foldvals)>0)),
              'local_null_tail':tail(sims,val),'largest_absolute_training_posterior_mean':maxcoef,
              'template_singular_values':np.linalg.svd(T,compute_uv=False).tolist()}
    mcs=np.array(mc);observed_max=max(v['sum_delta_log_score'] for v in results.values())
    maxnoise=mcs.max(axis=0);threshold=float(np.quantile(maxnoise,.95,method='higher'))
    power={}
    for fam,T in templates.items():
        for j in range(3):
            for amp in [.02,.05]:
                sims=noise+amp*T[:,j]
                scores=np.array([.5*np.einsum('bi,ij,bj->b',sims,Q,sims)+k for Q,k in operators.values()])
                power[f'{fam}_mode{j}_amp{amp:g}']={'known_direction_snr':float(abs(amp)*np.linalg.norm(T[:,j])),
                 'global_5percent_detection_fraction':float(np.mean(scores.max(axis=0)>=threshold))}
    pd.DataFrame(per).to_csv(OUT/'residual-transfer-object-scores.csv',index=False)
    pd.DataFrame(info).to_csv(OUT/'residual-transfer-objects.csv',index=False)
    np.savez_compressed(OUT/'residual-transfer-arrays.npz',residual=r,**templates,**{k.replace('.','p')+'_Q':q for k,(q,c) in operators.items()})
    inputs += [OUT/'residual-protocol.md',ROOT/'scripts/phase2/independent_flux/engine.py',ROOT/'scripts/phase2/independent_flux/fit.py']
    out={'scope':'Retrospective cross-object prediction of pre-existing epoch-held-out residuals; fixed local Gaussian design; no physical correction or cosmology inference.',
      'objects':len(info),'heldout_epochs':n,'baseline_conditional_chi2':float(r@r),
      'baseline_chi2_qualification':'Not an exact chi-square null: mask/covariance/fitted design and selection were conditioned on.',
      'results':results,'maximum_score_gain':observed_max,
      'maximum_search_null_tail':tail(maxnoise,observed_max),'global_95percent_threshold':threshold,
      'injection_power':power,'seed':seed,'null_draws':len(noise),
      'interpretation_limits':['Shared SALT/calibration errors across objects omitted from available covariance.',
       'Published accepted mask and frozen covariance use full-data information.',
       'Templates/Jacobians are evaluated at archived fits; Gaussian simulations do not reestimate those quantities.',
       'No prospective selection or nonlinear flux perturbation pipeline is run.',
       'Ridge scales are diagnostic regularization, not empirical systematic priors.'],
      'inputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},
      'source_sha256':sha(__file__),'python':platform.python_version(),'numpy':np.__version__}
    (OUT/'residual-transfer.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['objects','heldout_epochs','baseline_conditional_chi2','results','maximum_score_gain','maximum_search_null_tail','global_95percent_threshold','injection_power']},indent=2))

if __name__=='__main__':main()
