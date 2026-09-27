"""Independent block-conditioning audit; no cosmological model or fit.

Synthetic Gaussian systems establish normalization and fixed-contrast semantics.
Released-covariance checks, when enabled below, use arbitrary positive distance
curves only to verify likelihood algebra, never to infer cosmology.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
from contextlib import ExitStack
from unittest.mock import patch
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.linalg import solve
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT/'studies/unified_cosmology/results/inference/calibration-holdout-review.json'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rel(path):
    return str(Path(path).resolve().relative_to(ROOT))


def flat_logl(r, cov):
    one = np.ones(len(r))
    x = solve(cov, np.column_stack((r, one)), assume_a='gen')
    a = one @ x[:, 1]
    chi = r @ x[:, 0] - (one @ x[:, 0])**2/a
    sign, ld = np.linalg.slogdet(cov)
    assert sign == 1
    return -.5*(chi + ld + np.log(a) + (len(r)-1)*np.log(2*np.pi))


def proper_logl(r, cov):
    sign, ld = np.linalg.slogdet(cov)
    assert sign == 1
    return -.5*(r @ solve(cov, r, assume_a='gen') + ld + len(r)*np.log(2*np.pi))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUT)
    args = parser.parse_args()
    import camb
    import cobaya.model
    def forbidden(*unused, **kwargs):
        raise AssertionError('No cosmological model, background or spectra permitted in this review.')
    guard = ExitStack()
    for name in ('get_results', 'get_background', 'get_transfer_functions'):
        guard.enter_context(patch.object(camb, name, forbidden))
    guard.enter_context(patch.object(cobaya.model, 'get_model', forbidden))
    import calibration_holdout as producer
    dependencies = [Path(__file__), Path(producer.__file__), producer.DESIGN,
                    ROOT/'studies/unified_cosmology/code/inference/anchored_adapter.py',
                    ROOT/'studies/unified_cosmology/code/distance_ladder/calibration_interface.py']
    source_hashes = {rel(p):sha(p) for p in dependencies}
    rng = np.random.default_rng(273761)
    errors = {k: 0. for k in ('log_density', 'determinant_identity', 'contrast_mean',
                              'contrast_variance', 'offset_invariance', 'direct_M_quadrature',
                              'producer_A', 'producer_b', 'producer_V0', 'producer_V',
                              'producer_v', 'producer_sigma', 'producer_log_density',
                              'producer_contrast', 'producer_logcdf', 'producer_logsf',
                              'producer_offset_invariance', 'noncalibrator_model_shift',
                              'producer_noncalibrator_model_shift')}
    rows = []
    for i in range(64):
        nn, nc = int(rng.integers(3, 18)), int(rng.integers(1, 9))
        n = nn+nc
        mat = rng.normal(size=(n,n))
        cov = mat@mat.T/n + np.diag(rng.uniform(.2, 1., n))
        r = rng.normal(size=n)
        cnn, ccn, ccc = cov[:nn,:nn], cov[nn:,:nn], cov[nn:,nn:]
        A = solve(cnn, ccn.T, assume_a='gen').T
        oneN, oneC = np.ones(nn), np.ones(nc)
        b = oneC-A@oneN
        V0 = ccc-A@ccn.T
        p1 = solve(cnn, oneN, assume_a='gen')
        aN = oneN@p1
        mN = p1@r[:nn]/aN
        t = r[nn:]-A@r[:nn]-b*mN
        V = V0+np.outer(b,b)/aN
        full = flat_logl(r,cov)
        noncal = flat_logl(r[:nn],cnn)
        conditional = proper_logl(t,V)
        errors['log_density'] = max(errors['log_density'], abs(full-noncal-conditional))
        afull = np.ones(n)@solve(cov,np.ones(n),assume_a='gen')
        det_error = abs(np.linalg.slogdet(cov)[1]+np.log(afull)
                        -np.linalg.slogdet(cnn)[1]-np.log(aN)-np.linalg.slogdet(V)[1])
        errors['determinant_identity'] = max(errors['determinant_identity'],det_error)
        # The differential calibration offset is along one_C. The vector b
        # instead describes the latent M response at fixed noncalibrator data.
        q = solve(V,oneC,assume_a='gen')
        v = q/(oneC@q)
        sigma = 1/np.sqrt(oneC@q)
        errors['contrast_mean'] = max(errors['contrast_mean'],abs(v@oneC-1))
        errors['contrast_variance'] = max(errors['contrast_variance'],abs(v@V@v-sigma**2))
        shift = 53.27
        shifted_mN = p1@(r[:nn]+shift)/aN
        shifted_t = r[nn:]+shift-A@(r[:nn]+shift)-b*shifted_mN
        errors['offset_invariance'] = max(errors['offset_invariance'],float(np.max(abs(shifted_t-t))))
        # mu_N -> mu_N+k means r_N -> r_N-k, with fixed calibrated values.
        k=.31
        rn_k=r[:nn]-k
        t_k=r[nn:]-A@rn_k-b*(p1@rn_k/aN)
        errors['noncalibrator_model_shift']=max(errors['noncalibrator_model_shift'],
                                                float(np.max(abs(t_k-t-k))),abs(v@t_k-v@t-k))
        # Independently permute rows so implementation cannot rely on contiguous blocks.
        permutation = rng.permutation(n)
        mask = np.arange(n)[permutation] >= nn
        order_n = permutation[~mask]
        order_c = permutation[mask]-nn
        kernel = producer.GaussianHoldout(cov[np.ix_(permutation,permutation)],mask)
        for key, expected in [('A',A[np.ix_(order_c,order_n)]), ('b',b[order_c]),
                              ('V0',V0[np.ix_(order_c,order_c)]), ('V',V[np.ix_(order_c,order_c)]),
                              ('v',v[order_c]), ('sigma',sigma)]:
            errors['producer_'+key] = max(errors['producer_'+key],float(np.max(abs(getattr(kernel,key)-expected))))
        fixed_v = kernel.v.copy()
        actual = kernel.evaluate(r[permutation])
        for key, expected in [('full_loglike',full), ('noncalibrator_loglike',noncal),
                              ('conditional_calibrator_loglike',conditional)]:
            errors['producer_log_density'] = max(errors['producer_log_density'],abs(actual[key]-expected))
        errors['producer_contrast'] = max(errors['producer_contrast'],abs(actual['calibration_contrast_mag']-v@t))
        for key, expected in [('logcdf',norm.logcdf(v@t/sigma)), ('logsf',norm.logsf(v@t/sigma))]:
            errors['producer_'+key] = max(errors['producer_'+key],abs(actual['conditional_'+key]-expected))
        shifted = kernel.evaluate(r[permutation]+shift)
        for key in ('full_loglike','noncalibrator_loglike','conditional_calibrator_loglike','calibration_contrast_mag'):
            errors['producer_offset_invariance'] = max(errors['producer_offset_invariance'],abs(shifted[key]-actual[key]))
        changed_c = r[permutation].copy(); changed_c[mask] += rng.normal(size=nc)*100
        changed = kernel.evaluate(changed_c)
        assert changed['noncalibrator_loglike'] == actual['noncalibrator_loglike']
        assert np.array_equal(kernel.v,fixed_v)
        shifted_n=r[permutation].copy();shifted_n[~mask]-=k
        changed_n=kernel.evaluate(shifted_n)
        errors['producer_noncalibrator_model_shift']=max(errors['producer_noncalibrator_model_shift'],
            abs(changed_n['calibration_contrast_mag']-actual['calibration_contrast_mag']-k),
            abs(changed_n['noncalibrator_loglike']-actual['noncalibrator_loglike']))
        # Independent scalar posterior-M integral before completing its square.
        # Standard-normal x parametrizes M=mN+x/sqrt(aN).
        if i < 8:
            center = conditional
            def f(x):
                res = r[nn:]-A@r[:nn]-b*(mN+x/np.sqrt(aN))
                return np.exp(norm.logpdf(x)+proper_logl(res,V0)-center)
            integral, err = quad(f,-12,12,epsabs=1e-12,epsrel=1e-12)
            errors['direct_M_quadrature'] = max(errors['direct_M_quadrature'],abs(np.log(integral)))
        rows.append({'N':nn,'C':nc,'log_density_error':float(full-noncal-conditional),
                     'amplitude':float(v@t),'sigma':float(sigma)})
    # Known mixture demonstrates why tail probabilities must be averaged first.
    # Same observed scalar 0; predicted component means -3,+3, common sigma=1.
    weights = np.array([.5,.5]); residual = np.array([3.,-3.])
    cdf = float(weights@norm.cdf(residual))
    tail = 2*min(cdf,1-cdf)
    wrong = float(weights@(2*np.minimum(norm.cdf(residual),norm.sf(residual))))
    independent_cdf = quad(lambda y: .5*norm.pdf(y,-3,1)+.5*norm.pdf(y,3,1),
                           -np.inf,0,epsabs=1e-12,epsrel=1e-12)[0]
    assert max(errors.values()) < 1e-10
    assert abs(cdf-independent_cdf) < 1e-12 and tail == 1 and wrong < .003
    # Invalid inputs must be refused. b=0 is valid for the one_C contrast.
    refusals = []
    def refused(name, operation):
        try:
            operation()
        except (AssertionError, ValueError, np.linalg.LinAlgError):
            refusals.append(name)
        else:
            raise AssertionError('Did not refuse '+name)
    zero_b = np.array([[1.,0.,.5],[0.,1.,.5],[.5,.5,2.]])
    mask3 = np.array([False,False,True])
    zero_b_kernel=producer.GaussianHoldout(zero_b,mask3)
    assert np.array_equal(zero_b_kernel.b,np.zeros(1))
    assert abs(zero_b_kernel.sigma-np.sqrt(1.5))<1e-14
    zero_baseline=zero_b_kernel.evaluate([.2,.4,.8])
    zero_shifted=zero_b_kernel.evaluate([.2-.31,.4-.31,.8])
    assert abs(zero_shifted['calibration_contrast_mag']-zero_baseline['calibration_contrast_mag']-.31)<1e-14
    refused('nonsquare_covariance',lambda:producer.GaussianHoldout(np.eye(3)[:2],mask3))
    refused('nonfinite_covariance',lambda:producer.GaussianHoldout(np.eye(3)*np.nan,mask3))
    refused('nonboolean_mask',lambda:producer.GaussianHoldout(np.eye(3),mask3.astype(int)))
    refused('indefinite_covariance',lambda:producer.GaussianHoldout(-np.eye(3),mask3))
    kernel = producer.GaussianHoldout(np.eye(3),mask3)
    refused('wrong_residual_shape',lambda:kernel.evaluate(np.zeros(4)))
    refused('nonfinite_residual',lambda:kernel.evaluate([0.,np.inf,0.]))
    copied_C=np.eye(3); copied_mask=mask3.copy()
    copied = producer.GaussianHoldout(copied_C,copied_mask)
    previous=copied.evaluate([0.,1.,2.])
    copied_C[:]=np.nan; copied_mask[:]=False
    assert copied.evaluate([0.,1.,2.])==previous
    # Four chronological groups, >=10 points in each of40 blocks.
    groups = np.repeat(np.arange(4),500)
    z = np.tile(np.array([-3.,3.]),1000)
    predictive_rows = [{'conditional_logcdf':float(norm.logcdf(x)),
                        'conditional_logsf':float(norm.logsf(x)),
                        'conditional_calibrator_loglike':-12.} for x in z]
    prediction = producer.predictive_summary(predictive_rows,np.zeros(2000),groups)
    assert prediction['status']=='qualified_predictive_precision'
    assert abs(prediction['calibration_contrast_predictive_tail']['two_sided_equal_tail_probability']-1)<1e-14
    constant = producer.contribution_precision(np.zeros(2000),np.full(2000,-700.),groups,
                                               json.loads(producer.DESIGN.read_text())['predictive_precision'])
    assert constant['qualified_predictive_precision']
    assert max(constant['relative_batch_MCSE'].values())<1e-12
    assert constant['relative_independent_chain_MCSE']<1e-12
    chain_values=np.repeat([.2,.2,.2,.5],500)
    chain_shift=producer.contribution_precision(np.zeros(2000),np.log(chain_values),groups,
                                                json.loads(producer.DESIGN.read_text())['predictive_precision'])
    # Independently from chain means: SE(mean)/mean = std([.2,.2,.2,.5])/sqrt(4)/.275.
    expected_chain_mcse=np.std([.2,.2,.2,.5],ddof=1)/np.sqrt(4)/np.mean([.2,.2,.2,.5])
    assert abs(expected_chain_mcse-3/11)<1e-15
    assert abs(chain_shift['relative_independent_chain_MCSE']-expected_chain_mcse)<1e-13
    assert not chain_shift['qualified_predictive_precision']
    assert max(chain_shift['relative_batch_MCSE'].values())<.2
    assert 'relative_independent_chain_MCSE' in chain_shift['failed_gates']
    concentrated_rows = [{'conditional_logcdf':float(norm.logcdf(x)),
                          'conditional_logsf':float(norm.logsf(x)),
                          'conditional_calibrator_loglike':-12.} for x in np.r_[0.,np.full(1999,-40.)]]
    concentrated=producer.predictive_summary(concentrated_rows,np.zeros(2000),groups)
    assert concentrated['calibration_contrast_predictive_tail'] is None
    assert concentrated['status']=='predictive_precision_requires_followup'
    # Independent influence/block-sum calculation with unequal weights/integrands.
    lw=rng.normal(0,.1,2000); integrand=rng.normal(0,.07,2000)-3
    from scipy.special import logsumexp
    limits=json.loads(producer.DESIGN.read_text())['predictive_precision']
    cp=producer.contribution_precision(lw,integrand,groups,limits)
    w=np.exp(lw-logsumexp(lw)); values=np.exp(integrand)
    rho=w*values/(w@values)
    batch_error=0.
    for count in (10,20,40):
        blocks=[]
        for group in range(4):
            positions=np.flatnonzero(groups==group)
            sizes=np.full(count,len(positions)//count)
            sizes[:len(positions)%count]+=1
            edges=np.r_[0,np.cumsum(sizes)]
            blocks.extend(positions[edges[j]:edges[j+1]] for j in range(count))
        sums=np.array([np.sum(w[idx]*(values[idx]/(w@values)-1)) for idx in blocks])
        expected=np.sqrt(len(sums)/(len(sums)-1)*np.sum((sums-sums.mean())**2))
        batch_error=max(batch_error,abs(expected-cp['relative_batch_MCSE'][str(count)]))
    assert batch_error<1e-13
    independent_chain_sums=np.array([np.sum(w[groups==g]*(values[groups==g]/(w@values)-1)) for g in range(4)])
    expected_chain_error=np.sqrt(4*np.var(independent_chain_sums,ddof=1))
    chain_mcse_error=abs(expected_chain_error-cp['relative_independent_chain_MCSE'])
    assert chain_mcse_error<1e-13
    # Actual covariance/interface, but twelve fixed arbitrary curves, no cosmology.
    released=producer.ReleasedHoldout()
    release_inputs={key:item['sha256'] for key,item in released.release.input_records.items()}
    release_errors={'full_to_released':0.,'noncal_plus_conditional':0.,'H0_scaling_contrast_sign':0.}
    zz=released.release.z_hd_noncalibrator
    for i in range(12):
        da=4500*zz/(1+zz)*(1+(i-5.5)*.006*zz/(1+zz)+.004*np.sin((i+1)*zz))
        result=released.evaluate(da)
        release_errors['full_to_released']=max(release_errors['full_to_released'],abs(result['released_full_loglike_closure']))
        release_errors['noncal_plus_conditional']=max(release_errors['noncal_plus_conditional'],abs(result['full_loglike_closure']))
        # At unchanged shape, H0 -> 1.05 H0 gives D_A -> D_A/1.05.
        # These are algebraic copies of the synthetic curves, not model calls.
        scaled=released.evaluate(da/1.05)
        difference=scaled['calibration_contrast_mag']-result['calibration_contrast_mag']
        assert difference<0
        release_errors['H0_scaling_contrast_sign']=max(release_errors['H0_scaling_contrast_sign'],
                                                       abs(difference+5*np.log10(1.05)))
    assert max(release_errors.values())<1e-7
    assert all(sha(ROOT/path)==expected for path,expected in {**source_hashes,**release_inputs}.items())
    guard.close()
    output = {
        'status':'passed_independent_calibration_holdout_review',
        'seed':273761,'synthetic_systems':64,'released_fixed_distance_curves':12,
        'released_algebraic_H0_scaling_controls':12,
        'observational_fits':0,'background_calls':0,'spectra_calls':0,'model_constructions':0,
        'max_errors':errors,'cases':rows,
        'refused_cases':refusals,'fixed_contrast_unchanged':True,
        'zero_b_SPD_case_accepts_identifiable_common_calibration_contrast':True,
        'withheld_value_changes_leave_noncalibrator_loglike_bitwise_unchanged':True,
        'input_covariance_and_mask_are_copied':True,
        'released_normalization_max_errors':release_errors,
        'predictive_controls':{'balanced_mixture':prediction,'constant_integrand':constant,
                               'concentrated_tail':concentrated,'independent_batch_MCSE_max_error':batch_error,
                               'chain_constant_mismatch':chain_shift,
                               'chain_constant_mismatch_independent_expected_MCSE':float(expected_chain_mcse),
                               'independent_chain_MCSE_max_error':float(chain_mcse_error)},
        'mixture_control':{'predictive_CDF':cdf,'direct_integral_CDF':independent_cdf,
                           'correct_equal_tail_fraction':tail,
                           'incorrect_average_component_equal_tails':wrong},
        'interpretation':[
            'Proper C|N,theta factor exactly multiplies the same improper-flat-M noncalibrator likelihood.',
            'v is fixed when covariance and row selection are fixed; v^T one_C=1 and Var(v^T t)=sigma^2.',
            'Shifting only mu_N by k shifts t by k one_C and the calibration contrast by k; H0 rise at fixed shape gives k=-5 log10(H0_new/H0_old).',
            'Tail uses the CMB+BAO+N posterior, never a posterior conditioned on C.',
            'Equal-tail predictive diagnostic is not a frequentist p value or a density-ranked tail.',
            'b=0 is allowed: b describes latent common-M response, not the direction of a differential calibration offset.',
            'This tests the released conditional Gaussian construction, not physical independence of its covariance inputs.'
        ],
        'source_sha256':source_hashes,'input_sha256':release_inputs,
        'pre_observation_estimand_correction':{
            'old_direction':'b=one_C-A one_N',
            'old_review_sha256':'d24d966d4569cbea7b66f645b59fbd92bb60e796855a516aa27f61b252c2aea0',
            'old_validator_sha256':'aaad9bc08afe8238bb37acda4c0dd9690af1409f5a6b3a05e2ed90c195b7dad6',
            'preserved_files':'.work/unified-cosmology/calibrator-conditional-review/latent-M-contrast-preserved.json',
            'reason':'The old direction was a valid fixed linear contrast, but not the intended common calibrator/noncalibrator offset. It was replaced before observational bridge evaluation.'},
        'pre_observation_precision_correction':{
            'old_review_sha256':'1c5609fa27a1d4d62cea113b4f3ad7857aa8bf4286f145e0004d2b0e22b691f6',
            'preserved_files':'.work/unified-cosmology/calibrator-conditional-review/pre-chain-precision-preserved.json',
            'reason':'Within-chain block subdivision can hide a chain-constant predictive-integrand discrepancy. Independent-chain MCSE now uses the same unchanged 0.2 precision threshold.'},
        'dependencies':{'numpy':np.__version__,'scipy':scipy.__version__},
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':output['status'],'max_errors':errors,'mixture_control':output['mixture_control']}))


if __name__=='__main__':main()
