"""Independent likelihood algebra, synthetic coverage, and q/j equation checks."""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.integrate import quad
from likelihood import ROOT, ReleasedDistances, expansion_diagnostics


def main():
    rng = np.random.default_rng(270926)
    like = ReleasedDistances()
    distance = 10**((like.observed-25)/5)
    class Provider:
        def get_angular_diameter_distance(self, z):
            # Generated data need not have identical distances at repeated zHD;
            # use a smooth positive synthetic background instead.
            return 3000*z/(1+z)/(1+.3*z)
    like.provider = Provider()
    da = like.provider.get_angular_diameter_distance(like.background_z)[like.inverse_z]
    mu = 5*np.log10(da*(1+like.z)*(1+like.zhel))+25
    residual = mu-like.observed
    c = np.load(like.data_file)['covariance']
    chol = cholesky(c,lower=True)
    one = solve_triangular(chol,np.ones(len(c)),lower=True)
    white = solve_triangular(chol,residual,lower=True)
    best = np.dot(one,white)/np.dot(one,one)
    independent_chi2 = np.sum((white-best*one)**2)
    derived = {}
    lp = like.logp(_derived=derived)
    assert abs(derived['sn_chi2']-independent_chi2)<1e-7
    like.observed += 100
    lp_offset = like.logp()
    like.observed -= 100
    assert abs(lp-lp_offset)<1e-7
    # Numerical 1D integration of the constant offset at the optimum, stripping
    # the known Gaussian residual term to avoid underflow.
    amplitude = np.dot(one,one)
    numerical = quad(lambda t:np.exp(-.5*amplitude*t*t),-10/np.sqrt(amplitude),10/np.sqrt(amplitude),epsabs=1e-11)[0]
    integral_relative = abs(numerical/np.sqrt(2*np.pi/amplitude)-1)
    assert integral_relative<1e-10
    # Exact Gaussian injections for joint intercept+luminosity-slope GLS. This
    # validates statistical algebra under this model, not survey-model adequacy.
    design = np.column_stack([np.ones(len(c)),like.evolution])
    x = solve_triangular(chol,design,lower=True)
    fisher = x.T@x
    inverse = np.linalg.inv(fisher)
    projector = np.linalg.solve(fisher,x.T)
    draws = np.array([.13,.08])[:,None] + projector@rng.normal(size=(len(c),4096))
    sigma = np.sqrt(np.diag(inverse))
    zscore = (draws-np.array([.13,.08])[:,None])/sigma[:,None]
    coverage = np.mean(abs(zscore)<1.959963984540054,axis=1)
    assert np.all(abs(coverage-.95)<.015)
    # Independent analytic flatLCDM including radiation: q=.5*m+r-Lambda
    # divided byE², j=1+2r/E². q,j do not depend on H0 units.
    redshift = np.concatenate([np.arange(5)*.001+a for a in [0.,.5,1.]])
    om, radiation = .315,.00009
    de = 1-om-radiation
    e2 = om*(1+redshift)**3+radiation*(1+redshift)**4+de
    diagnostics = expansion_diagnostics(redshift,67.4*np.sqrt(e2))
    errors = {}
    for zz,label in [(0,'0'),(.5,'05'),(1,'1')]:
        m,r = om*(1+zz)**3,radiation*(1+zz)**4
        expectedq = (.5*m+r-de)/(m+r+de)
        expectedj = 1+2*r/(m+r+de)
        errors['q'+label] = abs(diagnostics['q'+label]-expectedq)
        errors['j'+label] = abs(diagnostics['j'+label]-expectedj)
    assert max(errors.values())<1e-6
    report = {'status':'passed','independent_sn_chi2_absolute_error':abs(derived['sn_chi2']-independent_chi2),
              'constant_offset_invariance':abs(lp-lp_offset),'offset_integral_relative_error':integral_relative,
              'conditional_gaussian_injections':4096,'95_percent_coverage_intercept_and_slope':coverage.tolist(),
              'standardized_estimate_means':zscore.mean(axis=1).tolist(),
              'q_jerk_absolute_errors':errors,
              'scope':'Independent algebra and coverage conditional on released Gaussian covariance and assumed mean; does not validate physical survey selection or host model.',
              'code_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in Path(__file__).parent.glob('*.py')}}
    out = ROOT/'studies/unified_cosmology/results/inference'
    out.mkdir(parents=True,exist_ok=True)
    (out/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
