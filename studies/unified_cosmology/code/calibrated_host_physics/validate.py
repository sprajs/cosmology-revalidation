#!/usr/bin/env python3
"""Independent objective, endpoint, grid and synthetic checks of conditional fits."""
import ast
import hashlib
import json
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2
from model import ROOT,CODE,WORK,INPUT,Library,problem,solve,whiten,nonnegative,velocity_smooth

RESULT=ROOT/'studies/unified_cosmology/results/calibrated_host_physics'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def independent_extrema(case, maximum_columns=35):
    """Age-constrained quadratic profiles, without Charnes-Cooper variables.

    A deterministic reduced physical sublibrary keeps the independent solver
    inexpensive. This checks the same equations on actual data, not a global
    certificate for the larger production grid (whose dual gates are separate).
    """
    original=np.flatnonzero(case['mass']>0)
    selected=np.unique(original[np.linspace(0,len(original)-1,maximum_columns).astype(int)])
    selected=np.r_[selected,np.flatnonzero(case['mass']==0)]
    a=case['matrix'][:,selected];y=case['y'];cov=case['cov']
    ages=case['ages'][selected];mass=case['mass'][selected]
    conic=solve(a,y,cov,ages,mass)
    if conic['status']!='compatible':
        return dict(status='conditional_sublibrary_incompatible',chi2=conic['best_chi2'])
    aw,yw=whiten(a,y,cov)
    norm=np.linalg.norm(aw,axis=0);aa=aw/norm
    threshold=conic['threshold'];rows=[]
    def profile(target,sign):
        # The ratio inequality is homogeneous and linear in original mass:
        # age<=target iff (ages-target*mass).coeff<=0, for positive mass.
        g=sign*(ages-target*mass)/norm
        g/=np.linalg.norm(g)
        eligible=g<=0
        start=np.zeros(len(norm))
        start[eligible],value,_=nonnegative(aa[:,eligible],yw)
        if np.all(g>=0):
            # At a library age boundary the homogeneous constraint forces
            # all forbidden coefficients exactly to zero. Reduce that face
            # explicitly instead of asking a generic optimizer to discover
            # a constraint with no strictly feasible point.
            return value
        constraint={'type':'ineq','fun':lambda x:-g@x,'jac':lambda x:-g}
        ans=minimize(lambda x:.5*np.sum((aa@x-yw)**2),start,
            jac=lambda x:aa.T@(aa@x-yw),method='SLSQP',bounds=[(0,None)]*len(norm),
            constraints=constraint,options={'ftol':1e-10,'maxiter':4000})
        value=float(np.sum((aa@ans.x-yw)**2))
        assert ans.success and g@ans.x<1e-7, (ans.message,g@ans.x,target,sign)
        return value
    for endpoint,sign in [('lower',1.),('upper',-1.)]:
        lo=float(ages[mass>0].min());hi=float(ages[mass>0].max())
        boundary=lo if sign==1 else hi
        if profile(boundary,sign)<=threshold:
            value=boundary;iterations=0
        else:
            for iterations in range(24):
                middle=(lo+hi)/2;feasible=profile(middle,sign)<=threshold
                if (feasible and sign==1) or (not feasible and sign==-1):hi=middle
                else:lo=middle
            value=hi if sign==1 else lo
        delta=abs(value-conic[endpoint]['age_Gyr'])
        rows.append(dict(endpoint=endpoint,profile_status='all_quadratic_subproblems_solved',
            independent_age_Gyr=value,conic_age_Gyr=conic[endpoint]['age_Gyr'],
            absolute_difference_Gyr=delta,bisection_iterations=iterations))
        assert delta<5e-4,rows[-1]
    return dict(status='passed',columns=len(selected),endpoints=rows)


def main():
    records=json.loads((WORK/'fit-records.json').read_text())
    refined=json.loads((WORK/'fit-records-refined.json').read_text())
    checked={}
    for name in ['acquisition','stellar-library','stellar-library-refined','fits','fits-refined','independent-review','resolution-reference','summary']:
        record=json.loads((RESULT/f'{name}.json').read_text())
        for key in ['dependencies_sha256','output_sha256','fixture_sha256','data_files_sha256','source_sha256','input_sha256','code_sha256']:
            if not isinstance(record.get(key),dict):continue
            for path,digest in record.get(key,{}).items():
                assert sha(ROOT/path)==digest,(name,path)
                checked[path]=digest
    acquisition=json.loads((RESULT/'acquisition.json').read_text())
    for path,digest in [(CODE/'acquire.py',acquisition['code_sha256']),
                        (CODE/'design.json',acquisition['design_sha256']),
                        (CODE/'requirements-lock.txt',acquisition['requirements_sha256'])]:
        assert sha(path)==digest;checked[str(path.relative_to(ROOT))]=digest
    extension=list((WORK/'.venv').glob('lib/python*/site-packages/fsps/_fsps*.so'))
    assert len(extension)==1 and sha(extension[0])==acquisition['compiled_extension_sha256']
    checked[str(extension[0].relative_to(ROOT))]=sha(extension[0])
    endpoints=0;objective_error=[];endpoint_error=[];max_gap=0.;retry_count=0
    for row in records+refined:
        case=np.load(ROOT/row['fixture']);assert sha(ROOT/row['fixture'])==row['fixture_sha256']
        a,y,cov=case['matrix'],case['y'],case['cov']
        assert np.linalg.eigvalsh(cov).min()>0
        r=a@case['best_coefficients']-y
        direct=float(r@np.linalg.solve(cov,r))
        error=abs(direct-row['fit']['best_chi2']);objective_error.append(error)
        assert error<1e-7
        assert (direct<=row['fit']['threshold'])==row['fit']['compatible']
        assert row['fit']['status']!='numerical_failure',row
        for key in ['lower','upper']:
            if row['fit'][key] is None:continue
            ep=row['fit'][key];coeff=case[key+'_coefficients']
            assert len(coeff)>0 and min(coeff)>=0
            residual=a@coeff-y;direct=float(residual@np.linalg.solve(cov,residual))
            value=float(case['ages']@coeff/(case['mass']@coeff))
            assert direct<=row['fit']['threshold']+2e-5
            assert abs(value-ep['primal_age_Gyr'])<1e-6
            assert ep['duality_gap_Gyr']<=1e-5
            assert (ep['age_Gyr']<=value+1e-7 if key=='lower' else ep['age_Gyr']>=value-1e-7)
            endpoint_error.append(abs(direct-ep['chi2']));max_gap=max(max_gap,ep['duality_gap_Gyr'])
            endpoints+=1;retry_count+=len(ep['retries'])
    plan=json.loads((CODE/'implementation-design.json').read_text())
    idlist=sorted(p.stem for p in (INPUT/'spectra').glob('*.npz'))
    assert len(idlist)==55
    assert len(records)==55*len(plan['scenarios'])
    for scenario in plan['scenarios']:
        assert sorted(r['targetid'] for r in records if r['scenario']==scenario)==idlist
    index={(r['targetid'],r['scenario']):r for r in records}
    refinement=[]
    for r in refined:
        base=index[(r['targetid'],r['scenario'])]
        assert r['fit']['best_chi2']<=base['fit']['best_chi2']+1e-6
        comparison=dict(targetid=r['targetid'],scenario=r['scenario'],
            primary_status=base['fit']['status'],refined_status=r['fit']['status'],
            best_chi2_change=r['fit']['best_chi2']-base['fit']['best_chi2'])
        if r['fit']['compatible'] and base['fit']['compatible']:
            assert r['fit']['lower']['age_Gyr']<=base['fit']['lower']['age_Gyr']+2e-5
            assert r['fit']['upper']['age_Gyr']>=base['fit']['upper']['age_Gyr']-2e-5
            comparison.update(lower_age_change_Gyr=r['fit']['lower']['age_Gyr']-base['fit']['lower']['age_Gyr'],
                              upper_age_change_Gyr=r['fit']['upper']['age_Gyr']-base['fit']['upper']['age_Gyr'])
        refinement.append(comparison)
    independent=[]
    for tid in plan['refined_targetids'][:4]:
        for scenario in ['primary','original_free_emission']:
            row=index[(tid,scenario)]
            independent.append(dict(targetid=tid,scenario=scenario,
                **independent_extrema(np.load(ROOT/row['fixture']))))
    # Synthetic data use actual response matrices and the observed covariance.
    # On every draw whose true mean lies in the declared ellipsoid, the true
    # mixture age must lie in the returned set. This is a numerical recovery
    # check, not an empirical validation of the missing calibration covariance.
    rng=np.random.default_rng(927550)
    injection=[]
    row=index[(plan['refined_targetids'][0],'primary')]
    case=np.load(ROOT/row['fixture']);a=case['matrix'];cov=case['cov'];ages=case['ages'];mass=case['mass']
    for j in range(40):
        active=rng.choice(a.shape[1],3,replace=False);weights=rng.dirichlet(np.ones(3))
        coeff=np.zeros(a.shape[1]);coeff[active]=weights
        target=a@coeff
        amplitude=np.linalg.norm(case['y'])/np.linalg.norm(target)
        target*=amplitude;true_age=float(ages@coeff)
        noise=rng.normal(size=5);y=target+np.linalg.cholesky(cov)@noise
        fit=solve(a,y,cov,ages,mass);inside=float(noise@noise)<=chi2.ppf(.95,5)
        captured=fit['status']=='compatible' and fit['lower']['age_Gyr']-2e-5<=true_age<=fit['upper']['age_Gyr']+2e-5
        assert not inside or captured
        injection.append(dict(true_age_Gyr=true_age,true_mean_inside_noise_ellipsoid=bool(inside),recovered=bool(captured),status=fit['status']))
    zero=solve(a,np.zeros(5),cov,ages,mass)
    assert zero['lower']['age_Gyr']==ages.min() and zero['upper']['age_Gyr']==ages.max()
    gasrow=index[(plan['refined_targetids'][0],'original_free_emission')]
    gascase=np.load(ROOT/gasrow['fixture']);aa=gascase['matrix'];mm=gascase['mass']
    gasdata=aa[:,mm==0]@np.ones((mm==0).sum())
    gasfit=solve(aa,gasdata,gascase['cov'],gascase['ages'],mm)
    assert gasfit['lower']['age_Gyr']==gascase['ages'][mm>0].min()
    assert gasfit['upper']['age_Gyr']==gascase['ages'].max()
    for p in CODE.glob('*.py'):ast.parse(p.read_text())
    history=[('fit-initial-offset-failure.log','The initial unary-minus/floor-division offset-count error aborted before fitting; repaired to the pinned author DIA convention.'),
             ('validation-initial-slsqp-failure.log','A direct fractional-age SLSQP cross-check did not converge. Replaced by independently constrained convex quadratic age profiles; production conic fits unchanged.'),
             ('validation-boundary-slsqp-failure.log','Generic quadratic SLSQP failed on an age-boundary face without strict interior. That known face is now reduced exactly before NNLS; no scientific tolerance relaxed.')]
    failures=[dict(path=str((WORK/name).relative_to(ROOT)),sha256=sha(WORK/name),interpretation=meaning)
              for name,meaning in history if (WORK/name).exists()]
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),status='passed',
        identities_checked=len(checked),cases_checked=len(records)+len(refined),
        direct_covariance_objective_max_error=max(objective_error),
        original_coordinate_endpoints_checked=endpoints,
        endpoint_chi2_max_error=max(endpoint_error),maximum_duality_gap_Gyr=max_gap,
        rejected_solver_attempts_preserved=retry_count,
        independent_direct_ratio_optimizer=independent,refinement=refinement,
        synthetic_mixtures=injection,
        zero_flux_all_ages_recovered=True,gas_only_arbitrary_stellar_age_recovered=True,
        initial_implementation_failures_preserved=failures,
        scientific_scope='Only numerical and identity validation. These checks do not establish stellar-library accuracy, noise calibration, host/progenitor equivalence or cosmological luminosity evolution.',
        dependencies_sha256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),CODE/'model.py',CODE/'implementation-design.json',
            RESULT/'fits.json',RESULT/'fits-refined.json',RESULT/'stellar-library.json',RESULT/'stellar-library-refined.json',
            RESULT/'independent-review.json',RESULT/'summary.json',RESULT/'resolution-reference.json']},
        checked_sha256=checked)
    (RESULT/'validation.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['checked_sha256','synthetic_mixtures','refinement']},indent=2))


if __name__=='__main__':main()
