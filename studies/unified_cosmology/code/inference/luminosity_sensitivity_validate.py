"""Independent quadrature, dense GLS and known-weight checks of sensitivity math."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.linalg import cho_factor, cho_solve

from luminosity_sensitivity import (IntegratedLuminosity, log_normal_interval,
    exact_background_record, weight_diagnostics, DESIGN_PATH, RESULTS, ROOT, digest)
from late_geometry import sample_path


def dense_loglike(residual,covariance):
    factor = cho_factor(covariance,lower=True)
    precision = cho_solve(factor,np.eye(len(residual)))
    unit = np.ones(len(residual)); a = unit@precision@unit
    r = residual-residual.mean()
    q = r@precision@r-(r@precision@unit)**2/a
    return float(-.5*(q+np.linalg.slogdet(covariance)[1]+np.log(a)+(len(r)-1)*np.log(2*np.pi)))


def main():
    rng = np.random.default_rng(2727868)
    cases = []; offset_errors = []
    for index in range(24):
        n = 31; z = np.geomspace(.01,1.25,n)
        matrix = rng.normal(size=(n,n))
        covariance = .01*(matrix@matrix.T/n+.2*np.eye(n))
        observed = rng.normal(size=n); prediction = observed+rng.normal(0,.2,n)
        sn = IntegratedLuminosity(z,z,observed,covariance)
        result = sn.from_prediction(prediction); residual = prediction-observed
        baseline = dense_loglike(residual,covariance)
        errors = {'baseline':abs(result['baseline_SN_loglike']-baseline)}
        f = np.log1p(z)/np.log(2)
        log_integrand = lambda eps:dense_loglike(residual+eps*f,covariance)-baseline
        optimum = minimize_scalar(lambda eps:-log_integrand(eps),bounds=(-.5,.5),method='bounded')
        peak = max(log_integrand(optimum.x),log_integrand(-.5),log_integrand(.5))
        integral,error = quad(lambda eps:np.exp(log_integrand(eps)-peak),-.5,.5,epsabs=1e-12,epsrel=1e-12)
        numerical = np.log(integral)+peak
        errors['linear_quadrature'] = abs(numerical-result['log_likelihood_ratios']['linear'])
        for name,sigma in [('smooth01',.1),('smooth03',.3)]:
            augmented = covariance+sigma*sigma*sn.basis@sn.basis.T
            ratio = dense_loglike(residual,augmented)-baseline
            errors[name+'_dense_covariance'] = abs(ratio-result['log_likelihood_ratios'][name])
        assert max(errors.values())<1e-9,(index,errors)
        shifted = sn.from_prediction(prediction+1e6)
        offset_errors.append(max(abs(shifted['log_likelihood_ratios'][k]-result['log_likelihood_ratios'][k]) for k in result['log_likelihood_ratios']))
        assert offset_errors[-1]<2e-6
        cases.append(errors)
    tail_errors = []
    for low,high in [(-41.,-40.),(-10.,-9.99),(-.01,.01),(8.,8.01),(40.,41.)]:
        mode = min(max(0.,low),high)
        value,_ = quad(lambda x:np.exp(-.5*(x*x-mode*mode)),low,high,epsabs=1e-13,epsrel=1e-13)
        independent = np.log(value)-.5*mode*mode-.5*np.log(2*np.pi)
        tail_errors.append(abs(log_normal_interval(low,high)-independent))
    assert max(tail_errors)<2e-10
    z = np.linspace(.01,1.2,31); observed = .08*np.log1p(z)/np.log(2)+19.
    signed = IntegratedLuminosity(z,z,observed,.01*np.eye(31)).from_prediction(np.zeros(31))
    assert abs(signed['linear_conditional_untruncated_mean']-.08)<1e-12
    degenerate = IntegratedLuminosity(np.full(31,.2),np.full(31,.2),np.ones(31),np.eye(31)).from_prediction(np.zeros(31))
    assert max(abs(v) for v in degenerate['log_likelihood_ratios'].values())<1e-12

    native = json.loads((RESULTS/'modern-profile.json').read_text())
    data_path = sample_path('dovekie')
    with np.load(data_path) as data:
        sn = IntegratedLuminosity(*(data[k] for k in ['zHD','zHEL','MU','covariance']))
        # Independent full augmented-covariance validation on the actual SN design.
        prediction = data['MU']+np.linspace(-.07,.1,len(data['MU']))
        result = sn.from_prediction(prediction)
        residual = prediction-data['MU']; baseline = dense_loglike(residual,data['covariance'])
        real_covariance_errors = {name:abs(dense_loglike(residual,data['covariance']+sigma*sigma*sn.basis@sn.basis.T)-baseline-result['log_likelihood_ratios'][name])
                                 for name,sigma in [('smooth01',.1),('smooth03',.3)]}
    assert max(real_covariance_errors.values())<1e-7
    background_checks = []
    for model in ['lcdm','cpl']:
        for candidate in native['native'][model]['candidates']:
            result = exact_background_record(candidate['point'],{'model':model},sn)
            difference = result['baseline_SN_loglike']-candidate['loglikes']['released_sn']
            background_checks.append({'model':model,'candidate':candidate['candidate_label'],
                                      'SN_loglike_difference':difference})
    assert max(abs(r['SN_loglike_difference']) for r in background_checks)<1e-5
    design = json.loads(DESIGN_PATH.read_text()); gates = design['overlap_gates']
    n = 20000; x = rng.normal(size=n); shift = .4
    check = weight_diagnostics(x*shift-.5*shift**2,{'x':x},np.repeat(np.arange(4),n//4),gates)
    mean_error = check['weighted_summaries_for_diagnostics']['x']['mean']-shift
    assert abs(mean_error)<.035
    assert abs(check['raw_weight_ESS']/n-np.exp(-shift**2))<.02
    assert check['qualified_overlap'],check['failed_gates']
    constant = weight_diagnostics(np.zeros(n),{'x':x},np.repeat(np.arange(4),n//4),gates)
    assert constant['Pareto_k'] is None and constant['qualified_overlap']
    concentrated = weight_diagnostics(np.r_[100.,np.zeros(n-1)],{'x':x},np.repeat(np.arange(4),n//4),gates)
    assert not concentrated['qualified_overlap'] and 'raw_weight_ESS' in concentrated['failed_gates']
    collapsed = weight_diagnostics(np.r_[1000.,np.zeros(n-1)],{'x':x},np.repeat(np.arange(4),n//4),gates)
    assert not collapsed['qualified_overlap'] and collapsed['per_chain']['1']['raw_weight_ESS']==0
    assert collapsed['weighted_chain_stability']['x']['independent_chain_weighted_means']['1'] is None
    json.dumps(collapsed,allow_nan=False)
    sources = [Path(__file__),Path(__file__).with_name('luminosity_sensitivity.py'),DESIGN_PATH]
    output = {'status':'passed','seed':2727868,'random_GLS_cases':len(cases),
        'maximum_errors':{key:max(r[key] for r in cases) for key in cases[0]},
        'maximum_large_offset_ratio_difference':max(offset_errors),
        'remote_normal_interval_max_log_error':max(tail_errors),
        'positive_injection_recovered_coefficient':signed['linear_conditional_untruncated_mean'],
        'constant_modes_exactly_cancel_after_M':True,'actual_covariance_logratio_errors':real_covariance_errors,
        'native_background_SN_checks':background_checks,'known_weight_target_mean_error':mean_error,
        'known_weight_expected_ESS_fraction':float(np.exp(-shift**2)),
        'known_weight_observed_ESS_fraction':check['raw_weight_ESS']/n,
        'constant_and_concentrated_weight_gates_checked':True,
        'zero_weight_chains_fail_without_clipping_or_nonfinite_JSON':True,
        'source_sha256':{str(p.relative_to(ROOT)):digest(p) for p in sources},
        'SN_data_sha256':digest(data_path),'native_profile_result_sha256':digest(RESULTS/'modern-profile.json'),
        'scope':'Independent numerical identities and interface checks; does not establish cosmological overlap or posterior validity.'}
    (RESULTS/'luminosity-sensitivity-validation.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))


if __name__=='__main__':main()
