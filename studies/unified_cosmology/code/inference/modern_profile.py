"""Bounded penalized-profile candidates with native exact-CAMB verification.

Surrogate optimization supplies candidate coordinates only. It does not certify
native stationarity, posterior mass, Bayesian evidence or discovery significance.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('MKL_NUM_THREADS','1')
os.environ.setdefault('CLIPY_NOJAX','1')
from pathlib import Path
import datetime, hashlib, json, time
import numpy as np
from scipy.optimize import minimize
from cobaya.model import get_model
from modern_fast import configuration, identify
from late_geometry import sample_path

ROOT=Path(__file__).resolve().parents[4]
WORK=ROOT/'.work/unified-cosmology/inference/modern-profile'
RESULT=ROOT/'studies/unified_cosmology/results/inference/modern-profile.json'
SURROGATE=ROOT/'.work/unified-cosmology/inference/surrogate/cubic-0509.npz'
AUTHOR=ROOT/'studies/unified_cosmology/results/external_probes/author-configuration.json'
MAX_CALLS_PER_START=700
MAX_SECONDS_PER_START=480
MAX_FALLBACKS_PER_MODEL=2
BOX=4.


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


class BudgetReached(RuntimeError):pass


def evaluate(model,point,gaussians):
    started=time.monotonic();r=model.logposterior(point)
    if not np.isfinite(r.logpost):return None
    penalty=sum(((point[k]-mean)/sigma)**2 for k,(mean,sigma) in gaussians.items())
    likes=dict(zip(model.likelihood,map(float,r.loglikes)))
    score=-2*sum(likes.values())+penalty
    return {'point':{k:float(v)for k,v in point.items()},'loglikes':likes,'sum_loglikelihood':sum(likes.values()),'gaussian_calibration_penalty':float(penalty),'penalized_objective':float(score),'logposterior_with_prior_constants_audit_only':float(r.logpost),'logprior_with_constants_audit_only':float(sum(r.logpriors)),'objective_plus_twice_logposterior_constant':float(score+2*r.logpost),'derived':dict(zip(model.parameterization.derived_params(),map(float,r.derived))),'seconds':time.monotonic()-started}


def fit_model(label,author):
    info=configuration(model=label,evolution='none',sample='dovekie',calibration='official_planck',surrogate=SURROGATE)
    info['debug']=False
    identity=identify(info,sample_path('dovekie'),SURROGATE)
    gaussians={k:(v['prior']['loc'],v['prior']['scale'])for k,v in info['params'].items()if isinstance(v,dict)and isinstance(v.get('prior'),dict)and v['prior'].get('dist')=='norm'}
    assert set(gaussians)=={'A_planck','P_act','Tcal','Ecal'}
    source=ROOT/author['proposal_only']['path'];raw_names=source.open().readline().lstrip('#').split();cov=np.loadtxt(source)
    mean=author['proposal_only']['transformed_mean'];results=[]
    with get_model(info) as model:
        names=list(model.parameterization.sampled_params());ix=[raw_names.index(n)for n in names]
        centre=np.array([mean[n]for n in names]);chol=np.linalg.cholesky(cov[np.ix_(ix,ix)])
        theory=model.theory['spectral_surrogate'];cache={};all_valid=[]

        def coordinates(x):return dict(zip(names,map(float,centre+chol@x)))

        reference=evaluate(model,coordinates(np.zeros(len(names))),gaussians)
        assert reference is not None
        print(json.dumps({'model':label,'first_proposal_seconds':reference['seconds'],'first_penalized_objective':reference['penalized_objective'],'surrogate_calls':theory.surrogate_calls,'exact_fallback_calls':theory.exact_calls}),flush=True)
        rng=np.random.default_rng(927731 if label=='lcdm' else 927732)
        starts=[np.zeros(len(names)),rng.normal(0,.35,len(names))]
        for number,initial in enumerate(starts):
            begin=time.monotonic();calls=0;invalid=0;best=reference;budget_reason=None
            exact_before=theory.exact_calls;surrogate_before=theory.surrogate_calls

            def objective(x):
                nonlocal calls,invalid,best
                if calls>=MAX_CALLS_PER_START:raise BudgetReached('proposal_call_budget')
                if time.monotonic()-begin>MAX_SECONDS_PER_START:raise BudgetReached('wall_time_budget')
                if theory.exact_calls>=MAX_FALLBACKS_PER_MODEL:raise BudgetReached('exact_fallback_budget')
                calls+=1;key=np.asarray(x,dtype='<f8').tobytes()
                if key not in cache:
                    previous=theory.exact_calls;record=evaluate(model,coordinates(x),gaussians)
                    if record is not None:
                        record['used_exact_fallback']=theory.exact_calls>previous
                        record['optimizer_coordinates']=np.asarray(x).tolist()
                        all_valid.append(record)
                        if record['penalized_objective']<best['penalized_objective']:best=record
                    cache[key]=record
                record=cache[key]
                if record is None:invalid+=1;return 1e20
                if record['penalized_objective']<best['penalized_objective']:best=record
                if calls%100==0:print(json.dumps({'model':label,'start':number,'calls':calls,'best':best['penalized_objective'],'seconds':time.monotonic()-begin,'fallbacks':theory.exact_calls}),flush=True)
                return record['penalized_objective']

            try:
                fitted=minimize(objective,initial,method='L-BFGS-B',bounds=[(-BOX,BOX)]*len(names),options={'maxfun':MAX_CALLS_PER_START,'maxiter':90,'ftol':1e-10,'gtol':2e-4,'eps':2e-4,'maxls':16})
                status={'success':bool(fitted.success),'message':str(fitted.message),'iterations':int(fitted.nit),'nfev':int(fitted.nfev),'reported_objective':float(fitted.fun),'reported_coordinates':fitted.x.tolist()}
            except BudgetReached as e:
                budget_reason=str(e);status={'success':False,'message':budget_reason}
            result={'start':number,'initial_coordinates':initial.tolist(),'optimizer':status,'calls':calls,'invalid_prior_or_theory_evaluations':invalid,'seconds':time.monotonic()-begin,'exact_fallbacks':theory.exact_calls-exact_before,'surrogate_calls':theory.surrogate_calls-surrogate_before,'best_candidate':best,'budget_reason':budget_reason}
            results.append(result)
            (WORK/f'{label}-proposal-start-{number}.json').write_text(json.dumps(result,indent=2)+'\n')
            print(json.dumps({'model':label,'start_complete':number,'calls':calls,'objective':best['penalized_objective'],'optimizer':status}),flush=True)
        constants=[x['objective_plus_twice_logposterior_constant']for x in all_valid]+[reference['objective_plus_twice_logposterior_constant']]
        assert np.ptp(constants)<1e-8, 'Unaccounted parameter-dependent prior factor.'
        if all_valid:
            with (WORK/f'{label}-evaluations.jsonl').open('w')as f:
                for r in all_valid:f.write(json.dumps(r)+'\n')
        proposal_counts={'exact_fallbacks':theory.exact_calls,'surrogate_calls':theory.surrogate_calls,'valid_unique_likelihood_evaluations':len(all_valid)+1,'constant_prior_normalization_range':float(np.ptp(constants))}
    return {'configuration_identity':identity,'gaussian_calibration_priors':gaussians,'parameter_names':names,'author_mean':mean,'author_covariance_file':str(source.relative_to(ROOT)),'author_covariance_sha256':sha(source),'reference':reference,'starts':results,'proposal_counts':proposal_counts,'computational_trust_box_standardized_cholesky_coordinates':[-BOX,BOX],'scope':'Author moments set numerical starting/scaling coordinates only; box is a declared optimization restriction, not a change to the scientific target.'}


def exact_check(label,fit,extra=None):
    info=configuration(model=label,evolution='none',sample='dovekie',calibration='official_planck')
    info['debug']=False;identity=identify(info,sample_path('dovekie'))
    candidates=[('author_mean',fit['reference']['point'])]+[(f'surrogate_start_{r["start"]}',r['best_candidate']['point'])for r in fit['starts']]
    if extra is not None:candidates.append(('embedded_LCDM_candidate',dict(extra,w=-1.,wa=0.)))
    output=[]
    with get_model(info)as model:
        for name,point in candidates:
            record=evaluate(model,point,fit['gaussian_calibration_priors'])
            assert record is not None,(label,name)
            record['candidate_label']=name;output.append(record)
            print(json.dumps({'native_model':label,'candidate':name,'penalized_objective':record['penalized_objective'],'seconds':record['seconds']}),flush=True)
    return {'configuration_identity':identity,'candidates':output,'native_exact_evaluations':len(output),'best_evaluated_candidate':min(output,key=lambda x:x['penalized_objective']),'stationarity_certified':False}


def main():
    WORK.mkdir(parents=True,exist_ok=True)
    parent_files=[Path(__file__).with_name(p)for p in ['modern_run.py','modern_fast.py','spectral_surrogate.py','likelihood.py','late_geometry.py','target_identity.py']]+[Path(__file__).parents[1]/'external_probes'/p for p in ['modern_adapter.py','adapter.py','fast_lensing.py']]
    before={str(p.relative_to(ROOT)):sha(p)for p in parent_files}
    design={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'models':['lcdm','cpl'],'objective':'S = -2 sum(native component loglikelihoods) + sum[((Gaussian calibration parameter - prior mean)/sigma)^2]. Uniform-prior support is enforced, but prior-density normalization constants are not part of S.','fixed_sample':'Dovekie','calibration':'official_planck','luminosity_evolution':'none','surrogate':str(SURROGATE.relative_to(ROOT)),'surrogate_sha256':sha(SURROGATE),'starts_per_model':2,'proposal_calls_per_start_cap':MAX_CALLS_PER_START,'wall_seconds_per_start_cap':MAX_SECONDS_PER_START,'exact_fallbacks_per_model_cap':MAX_FALLBACKS_PER_MODEL,'native_verification_points_planned':7,'maximum_total_native_calls_including_fallbacks':11,'optimization_method':'L-BFGS-B, finite-difference step2e-4 in author Cholesky coordinates, trustbox±4; no active source changes','native_optimization':'None: evaluate candidate positions only; no exact stationarity certification.'}
    (WORK/'execution-design.json').write_text(json.dumps(design,indent=2)+'\n')
    authors=json.loads(AUTHOR.read_text())['models'];fits={};exact={}
    for label in ['lcdm','cpl']:
        fits[label]=fit_model(label,authors[label])
        extra=exact['lcdm']['best_evaluated_candidate']['point']if label=='cpl'else None
        exact[label]=exact_check(label,fits[label],extra)
        (WORK/f'{label}-complete.json').write_text(json.dumps({'proposal':fits[label],'native':exact[label]},indent=2)+'\n')
    for p,h in before.items():assert sha(ROOT/p)==h,f'Active target source changed: {p}'
    lcdm=exact['lcdm']['best_evaluated_candidate'];cpl=exact['cpl']['best_evaluated_candidate']
    embedded=next(r for r in exact['cpl']['candidates']if r['candidate_label']=='embedded_LCDM_candidate')
    nesting=abs(lcdm['penalized_objective']-embedded['penalized_objective'])
    assert nesting<1e-6,nesting
    out={'schema':'modern-penalized-profile-v1','status':'native_verified_candidates_not_certified_native_optima','design':design,'proposal':fits,'native':exact,'native_comparison':{'LCDM_minus_CPL_penalized_objective':lcdm['penalized_objective']-cpl['penalized_objective'],'LCDM_minus_CPL_negative_twice_loglikelihood':-2*(lcdm['sum_loglikelihood']-cpl['sum_loglikelihood']),'LCDM_minus_CPL_calibration_penalty':lcdm['gaussian_calibration_penalty']-cpl['gaussian_calibration_penalty'],'LCDM_embedding_exact_objective_discrepancy':nesting,'meaning':'Difference between best native-evaluated candidates. Not certified profile minima, posterior probability, Bayes factor, or sigma. Native component normalization is retained, so objective is not a global goodness-of-fit chi-square.'},'source_sha256':dict(before,**{str(Path(__file__).relative_to(ROOT)):sha(__file__)}),'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [AUTHOR,SURROGATE,sample_path('dovekie')]},'work_records_sha256':{str(p.relative_to(ROOT)):sha(p)for p in sorted(WORK.glob('*'))if p.is_file()}}
    RESULT.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out['native_comparison'],indent=2),flush=True)


if __name__=='__main__':main()
