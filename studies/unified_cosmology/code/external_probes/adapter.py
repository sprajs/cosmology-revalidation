"""Evaluated public CMB/BAO likelihood configuration; no distance-prior shortcut.

The caller may replace cosmological parameter priors and add a SN likelihood.
All observables then share the same CAMB provider and derived sound horizon.
"""
from pathlib import Path
import os

ROOT = Path(__file__).resolve().parents[4]
WORK = ROOT / '.work/unified-cosmology/external-probes'
PACKAGES = WORK / 'packages'


def external_info(mode='full', model='lcdm', bao=True, lensing=True):
    if mode not in ('full', 'lite') or model not in ('lcdm', 'cpl'):
        raise ValueError((mode, model))
    os.environ.setdefault('CLIPY_NOJAX', '1')
    high = ('planck_2018_highl_plik.TTTEEE' if mode == 'full' else
            'planck_2018_highl_plik.TTTEEE_lite_native')
    likes = {high: {}, 'planck_2018_lowl.TT': {}, 'planck_2018_lowl.EE': {}}
    if lensing:
        likes['planck_2018_lensing.native'] = {}
    if bao:
        likes['bao.desi_dr2'] = {}
    # Broad declared computational defaults, not likelihood-implied priors.
    params = {
        'H0': {'prior': {'min': 40, 'max': 100}, 'ref': 67.36, 'proposal': .3},
        'ombh2': {'prior': {'min': .005, 'max': .1}, 'ref': .02237, 'proposal': .00015},
        'omch2': {'prior': {'min': .001, 'max': .99}, 'ref': .1200, 'proposal': .001},
        'logA': {'prior': {'min': 1.61, 'max': 3.91}, 'ref': 3.044, 'proposal': .01, 'drop': True},
        'As': {'value': 'lambda logA: 1.e-10*np.exp(logA)'},
        'ns': {'prior': {'min': .8, 'max': 1.2}, 'ref': .9649, 'proposal': .004},
        'tau': {'prior': {'min': .01, 'max': .8}, 'ref': .0544, 'proposal': .006},
        'omegam': {'derived': True}, 'rdrag': {'derived': True},
    }
    if model == 'cpl':
        params['w'] = {'prior': {'min': -3, 'max': 1}, 'ref': -1, 'proposal': .04}
        params['wa'] = {'prior': {'min': -5, 'max': 3}, 'ref': 0, 'proposal': .1}
    else:
        params.update(w=-1., wa=0.)
    return {
        'packages_path': str(PACKAGES),
        'theory': {'camb': {'path': 'global', 'extra_args': {
            'dark_energy_model': 'ppf', 'mnu': .06,
            'num_massive_neutrinos': 1, 'nnu': 3.044, 'omk': 0.,
            'lens_potential_accuracy': 1, 'halofit_version': 'mead',
        }}},
        'likelihood': likes, 'params': params,
    }


def reference_point(model):
    """Deterministic prior-reference coordinates, not a claimed best fit."""
    import numpy as np
    values = model.prior.reference(random_state=np.random.default_rng(272601))
    return dict(zip(model.parameterization.sampled_params(), values))


if __name__ == '__main__':
    import argparse, yaml
    p = argparse.ArgumentParser()
    p.add_argument('--mode', choices=['full', 'lite'], default='full')
    p.add_argument('--model', choices=['lcdm', 'cpl'], default='lcdm')
    a = p.parse_args()
    print(yaml.safe_dump(external_info(a.mode, a.model), sort_keys=False))
