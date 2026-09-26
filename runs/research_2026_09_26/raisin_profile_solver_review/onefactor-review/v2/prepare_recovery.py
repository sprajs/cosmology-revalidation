"""Freeze a synthetic recovery design; generate no noise and launch no fits."""
from pathlib import Path
import hashlib,importlib.util,json,warnings
import numpy as np
from scipy.special import log_ndtr

ROOT=Path(__file__).resolve().parents[5];OUT=Path(__file__).resolve().parent
p=ROOT/'scripts/research_2026_09_26/onefactor_sign_likelihood.py'
s=importlib.util.spec_from_file_location('active',p);k=importlib.util.module_from_spec(s);s.loader.exec_module(k)
t=np.array([-10.,-3.,-1.,0.,1.,3.,10.]);expected=np.exp(-.5*t*t-.5*np.log(2*np.pi)-log_ndtr(t));err=float(np.max(abs(k.inverse_mills(t)-expected)))
assert err<1e-12
with warnings.catch_warnings():
    warnings.simplefilter('error')
    for a in [-1e12,-1e150]:
        try:k.log_orthant([a,a],[1.,1.],[.1,.1],[1.,1.])
        except ValueError:pass
        else:raise AssertionError('remote mode must be refused')
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
(OUT/'recovery-kernel.py').write_bytes(p.read_bytes())
source=ROOT/'runs/research_2026_09_26/astra_design/raisin_signed_refit/fixed-c-profile/global-profile/native-profiles.npz'
with np.load(source) as z:C=z['covariance_B'].copy()
d,v=k.decompose_covariance(C);h=.5*np.sqrt(np.diag(C))
np.savez_compressed(OUT/'recovery-design.npz',covariance=C,diagonal_variance=d,loading=v,template=h,true_amplitude=np.array(1.))
checks=dict(active_kernel_sha256=sha(p),Mills_moderate_max_abs_error=err,remote_modes_refused_with_warnings_as_errors=True,note='Original v2 review remains against its unchanged snapshot. This limited check covers active v2_1 hardening only.')
(OUT/'active-kernel-check.json').write_text(json.dumps(checks,indent=2)+'\n')
protocol=dict(status='Frozen design only; no recovery draws or estimator fits launched',truth=dict(amplitude=1.,template='h_i=0.5*sqrt(C_ii), all74coordinates; not a supernova SED or observed fitted mean'),
    covariance='Fixed B reference C only, diag(d)+vvT; conditional engineering generator, not established extraction covariance',realizations=128,
    randomness=dict(bit_generator='PCG64',new_noise_seed=202609260731,independent_archive_seed=202609260732,per_draw='One scalar Z then74epsilon; y=h+vZ+sqrt(d)*epsilon; same y across all estimators'),
    four_estimators=['full signed Gaussian','naive retained-positive Gaussian','joint censored known-cadence likelihood','full-sign-pattern-conditioned likelihood: joint censor minus full orthant log probability'],
    extra_control='One independently generated archived-positive mask; new signed noise on it, ordinary marginal Gaussian',amplitude_domain=[0.,4.],
    optimization=dict(analytic=['signed Gaussian','naive positive Gaussian','frozen-mask Gaussian'],proper_likelihoods='129 then257 point amplitude grids0..4; retain endpoints/every sampled local maximum; continuously polish all brackets; preserve ties/nonidentification',xtol=1e-8,fine_coarse_loglikelihood_tolerance=1e-7,fine_coarse_amplitude_tolerance=1e-5,quadrature_orders=[64,128,256],max_loglikelihood_quadrature_error=1e-8,independent_adaptive_checks='At truth and all candidate optima for both proper likelihoods; refuse unsupported cases rather than score as zero'),
    boundaries='Retain zero/upper boundary counts and unidentified all-negative conditional cases; no regenerated draws; distance undefined at zero, report finite subset membership',
    outputs=['per-draw noiseless mean/noise/y/signmask and archive mask','amplitude bias,RMSE,paired differences with MCSE','score at true amplitude and MCSE','deltaD=-2.5log10(a) separately with zero/failed/tied counts','minima,quadrature/support/resource failures and domain boundaries'],
    score_step=dict(central=1e-4,half=5e-5,absolute_derivative_agreement=1e-6),
    resource=dict(workers=1,initial_benchmark_realizations=8,maximum_elapsed_seconds=600,numerical_scalar_fits=256,analytic_amplitude_evaluations=384,preserve_partial=True),
    limitations='Finite-sample MLE bias/Jensen effects remain with a correct likelihood. No SN model fitting, population transport, event selection or physical/cosmology correction.',
    inputs_sha256={str(q.relative_to(ROOT)):sha(q) for q in [OUT/'recovery-kernel.py',OUT/'recovery-design.npz',source,OUT/'review.py',OUT/'result.json',Path(__file__)]})
(OUT/'recovery-protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
print(json.dumps(dict(active_check=checks,recovery_protocol_sha256=sha(OUT/'recovery-protocol.json')),indent=2))
