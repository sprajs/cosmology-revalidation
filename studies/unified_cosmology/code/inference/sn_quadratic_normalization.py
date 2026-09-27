"""Distinguish the residual quadratic from the author's marginalized score.

An independently integrated matter-plus-Lambda distance calculation and two
linear solvers check this accounting. No covariance rescaling or CMB calls.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve, lu_factor, lu_solve
from scipy.optimize import minimize_scalar
from scipy.integrate import quad

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DATA = ROOT/'.work/unified-cosmology/survey-selection/normalized/dovekie-total.npz'
AUTHOR = ROOT/'.work/unified-cosmology/external-probes/author-chains/Dovekie_cosmosis_likelihood.py'
ACQUISITION = ROOT/'studies/unified_cosmology/results/external_probes/author-acquisition.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def distances(z, zhel, matter, nodes=256):
    gx, gw = np.polynomial.legendre.leggauss(nodes)
    t = z[:, None]*(gx+1)/2
    radial = z/2*((1/np.sqrt(matter*(1+t)**3+1-matter))@gw)
    return 5*np.log10((1+zhel)*299792.458/70*radial)+25


def run():
    acquisition = json.loads(ACQUISITION.read_text())
    def hashes(value):
        if isinstance(value, dict):
            if value.get('path') == relative(AUTHOR):
                yield value
            for item in value.values():
                yield from hashes(item)
        elif isinstance(value, list):
            for item in value:
                yield from hashes(item)
    sources = list(hashes(acquisition))
    assert len(sources) == 1 and sources[0]['sha256'] == digest(AUTHOR)
    # Import only the reviewed numerical function, not the CosmoSIS module.
    parsed = ast.parse(AUTHOR.read_text())
    fn = [n for n in parsed.body if isinstance(n, ast.FunctionDef) and n.name == 'cov_log_likelihood']
    assert len(fn) == 1
    namespace = {'np': np}
    exec(compile(ast.Module(body=fn, type_ignores=[]), str(AUTHOR), 'exec'), namespace)
    official = namespace['cov_log_likelihood']
    with np.load(DATA) as data:
        z, zh, observed, covariance = [data[k] for k in ['zHD', 'zHEL', 'MU', 'covariance']]
    factor = cho_factor(covariance, lower=True)
    precision = cho_solve(factor, np.eye(len(z)))
    ones = np.ones(len(z)); u = precision@ones; information = float(ones@u)
    offset = float(np.log(information/(2*np.pi)))
    lu = lu_factor(covariance)
    def score(prediction):
        residual = prediction-observed
        residual -= residual.mean()
        return float(residual@precision@residual-(residual@u)**2/information)
    optimum = minimize_scalar(lambda om: score(distances(z, zh, om)),
                              bounds=(.01, .99), method='bounded', options={'xatol': 1e-11})
    assert optimum.success
    tests = []
    for matter in [.2, .3, .33, .4, float(optimum.x)]:
        prediction = distances(z, zh, matter)
        radial = np.array([quad(lambda t: 1/np.sqrt(matter*(1+t)**3+1-matter),
                                0, value, epsabs=1e-12, epsrel=1e-12)[0] for value in z])
        direct = 5*np.log10((1+zh)*299792.458/70*radial)+25
        residual = direct-observed; residual -= residual.mean()
        solution = lu_solve(lu, residual); unit = lu_solve(lu, ones)
        independent = float(residual@solution-(ones@solution)**2/(ones@unit))
        pure = score(prediction)
        author = float(-2*official(prediction, observed, precision))
        shifted = float(-2*official(prediction+1., observed+1., precision))
        checks = {'quadrature_magnitude_max_abs': float(np.max(abs(prediction-direct))),
                  'quadratic_LU_Cholesky_abs': abs(independent-pure),
                  'author_score_minus_quadratic_minus_normalization': author-pure-offset,
                  'common_observed_and_predicted_shift_error': shifted-author}
        assert checks['quadrature_magnitude_max_abs'] < 1e-10
        assert checks['quadratic_LU_Cholesky_abs'] < 1e-7
        assert abs(checks['author_score_minus_quadratic_minus_normalization']) < 1e-7
        assert abs(checks['common_observed_and_predicted_shift_error']) < 1e-7
        tests.append({'Omega_m': matter, 'residual_quadratic': pure,
                      'author_marginalized_score': author, 'checks': checks})
    bound = {relative(p): digest(p) for p in [DATA, AUTHOR, ACQUISITION]}
    return {'status': 'passed_residual_quadratic_normalization_accounting',
            'observations': len(z), 'fitted_geometry': 'Flat matter plus Lambda, radiation omitted; H0=70 is an arbitrary intercept.',
            'optimum_Omega_m': float(optimum.x), 'minimum_residual_quadratic': float(optimum.fun),
            'minimum_author_marginalized_score': float(optimum.fun+offset),
            'flat_magnitude_information': information, 'author_additive_log_A_over_2pi': offset,
            'paper_reference': {'url': 'https://arxiv.org/html/2511.07517v3#S10.T10',
                                'table': 10, 'printed_SN_only_LCDM_score': 1640.3,
                                'printed_Omega_m': '.330 +/- .015',
                                'BEAMS_effective_data_points': 1684,
                                'BEAMS_reference': 'Section 10.5 and Appendix C: sum of BBC BEAMS probabilities.'},
            'minimum_author_score_minus_printed_value': float(optimum.fun+offset-1640.3),
            'checks': tests, 'input_sha256': bound, 'source_sha256': {relative(__file__): digest(__file__)},
            'CMB_or_CAMB_background_calls': 0, 'covariance_modified': False,
            'interpretation': 'The author code names Q+log(A/2pi) chi2 after integrating the flat magnitude. Only Q has the fixed-geometry N-1 Gaussian quadratic reference. Printed best-fit scores and posterior-averaged Q also refer to different geometries. The paper separately reports 1684 effective data points from summed BBC BEAMS probabilities. This differs from the 1820-dimensional fixed Gaussian replication model; neither count may be silently substituted for the other. Accounting agreement does not validate covariance construction, and a small Gaussian-replication tail alone does not diagnose a problem with the original BEAMS mixture generator.'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    assert not args.output.exists()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: result[k] for k in ['status', 'optimum_Omega_m', 'minimum_residual_quadratic',
                    'author_additive_log_A_over_2pi', 'minimum_author_marginalized_score', 'minimum_author_score_minus_printed_value']}))


if __name__ == '__main__':
    main()
