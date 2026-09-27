"""Independently check all 47 constrained BayeSN synthetic coordinates.

Uses a synthetic linear forward kernel, no SED, photometry, or posterior fit.
The reference sums SciPy densities and directly integrates the broad D prior.
Run in the isolated BayeSN environment described in bayesn_heldout/README.md.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                  JAX_PLATFORM_NAME='cpu', JAX_ENABLE_X64='true',
                  XLA_FLAGS='--xla_cpu_multi_thread_eigen=false --xla_force_host_platform_device_count=1')
os.sched_setaffinity(0, [min(os.sched_getaffinity(0))])
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent/'bayesn_heldout'))

import jax.numpy as jnp
import numpy as np
from numpyro.infer.util import log_density
from scipy.integrate import quad
from scipy.stats import expon, norm, uniform
import model


def validate():
    rng = np.random.default_rng(273811)
    matrix = rng.normal(size=(3, 47))*.02

    class SyntheticForward:
        def forward(self, parameters):
            return jnp.asarray(matrix)@parameters

    data = dict(metadata=dict(mu_LCDM=35., sigma_external=.1),
                errors=[.4, .6, .8], flux=[1., -.2, .3])
    records = []
    for arm in ('LCDM', 'broad'):
        distances = ([34.6, 34.9, 35., 35.1, 35.4, 35.7] if arm=='LCDM'
                     else [19.8, 20.05, 34.8, 35., 49.95, 50.2])
        for index, distance in enumerate(distances):
            parameters = np.r_[distance, .05+.12*index, 1.3+.7*index, -.5+.2*index,
                               rng.normal(size=42), -8+5*index]
            sites = dict(D=parameters[0], AV=parameters[1], RV=parameters[2], theta=parameters[3],
                         epsilon_white=jnp.array(parameters[4:46]), tau=parameters[46])
            actual = float(log_density(model.numpyro_model(SyntheticForward(), data, arm), (), {}, sites)[0])
            if arm=='LCDM':
                prior_distance = norm.logpdf(distance, 35., np.hypot(.1, .088))
            else:
                integrated = quad(lambda mu: norm.pdf(distance, mu, .088)/30, 20, 50,
                                  points=[np.clip(distance, 20, 50)], epsabs=1e-14, epsrel=1e-12)[0]
                prior_distance = np.log(integrated)
            expected = (prior_distance + expon.logpdf(parameters[1], scale=.329)
                        + uniform.logpdf(parameters[2], loc=1.2, scale=4.8)
                        + norm.logpdf(parameters[3]) + norm.logpdf(parameters[4:46]).sum()
                        + uniform.logpdf(parameters[46], loc=-10, scale=30)
                        + norm.logpdf(data['flux'], matrix@parameters, data['errors']).sum())
            error = abs(actual-expected)
            assert error < 1e-10
            records.append(dict(arm=arm, D=distance, absolute_log_density_error=error))
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    return dict(status='passed_independent_47_coordinate_density_check', cases=records,
                max_error=max(row['absolute_log_density_error'] for row in records), native_SED_calls=0,
                source_sha256={str(Path(__file__).resolve().relative_to(ROOT)): sha(Path(__file__)),
                               str(Path(model.__file__).relative_to(ROOT)): sha(Path(model.__file__))},
                scope='Independent synthetic linear forward kernel; proper priors, coordinate ordering and signed Gaussian likelihood. Does not test mixing, SED physics, or observed data.')


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists(), 'Preserve previous validation records.'
    result = validate()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(result, indent=2))
