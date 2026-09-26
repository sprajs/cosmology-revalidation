"""Independent frozen-amplitude mass-threshold counterfactual verification."""
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RELEASE = ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/distances/w'
RUN = ROOT/'runs/research_2026_09_26/raisin_differential/mass-threshold'
OUT = ROOT/'runs/research_2026_09_26/raisin_mass_independent_verification'
INPUTS = {}


def record(path):
    INPUTS[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return path


def read(branch, variant):
    path = record(RELEASE/f'{branch}_dist/RAISIN_combined_FITOPT{variant:03d}.FITRES')
    lines = path.read_text().splitlines()
    columns = next(x.split()[1:] for x in lines if x.startswith('VARNAMES:'))
    rows = [dict(zip(columns, x.split()[1:])) for x in lines if x.startswith('SN:')]
    return np.array([r['CID'] for r in rows]), {key: np.array([float(r[key]) for r in rows]) for key in
             ['DLMAG', 'DLMAGERR', 'MASS_CORR', 'HOST_LOGMASS', 'zHD']}


def main():
    OUT.mkdir(exist_ok=False)
    record(Path(__file__))
    saved = dict(np.load(record(RUN/'counterfactual-arrays.npz')))
    claimed = json.loads(record(RUN/'counterfactual-result.json').read_text())
    record(RUN/'counterfactual-protocol.json')
    record(RUN/'code/raisin_cosmo/cosmo_sys.py')
    records = {}; maximum_error = 0.
    for branch in ('nir', 'optical', 'opticalnir'):
        ids, base = read(branch, 0); ids4, four = read(branch, 4); ids5, five = read(branch, 5)
        assert np.array_equal(ids, ids4) and np.array_equal(ids, ids5) and np.array_equal(ids, saved['CID'])
        assert np.array_equal(base['HOST_LOGMASS'], four['HOST_LOGMASS'])
        assert np.array_equal(base['zHD'], four['zHD']) and np.array_equal(base['zHD'], five['zHD'])
        gamma = 2*abs(four['MASS_CORR'][0])
        assert np.max(abs(abs(four['MASS_CORR'])-gamma/2)) < 1e-13
        expected_old_mass = np.where(four['HOST_LOGMASS'] > 10, gamma/2, -gamma/2)
        assert np.max(abs(expected_old_mass-four['MASS_CORR'])) < 1e-13
        new_mass = np.where(four['HOST_LOGMASS'] > 10.44, gamma/2, -gamma/2)
        repaired = four['DLMAG']+new_mass-four['MASS_CORR']
        changed = (four['HOST_LOGMASS'] > 10) & (four['HOST_LOGMASS'] <= 10.44)
        assert changed.sum() == 16 and np.sum(changed & (base['zHD'] < .1)) == 11
        old_delta = four['DLMAG']-base['DLMAG']; new_delta = repaired-base['DLMAG']
        center = lambda v, s: v - np.dot(v, 1/s**2)/np.sum(1/s**2)
        old4 = center(old_delta, four['DLMAGERR'])
        new4 = center(new_delta, four['DLMAGERR'])
        v5 = center(five['DLMAG']-base['DLMAG'], five['DLMAGERR'])
        old_cov = np.outer(old4, old4)+np.outer(v5, v5)
        new_cov = np.outer(new4, new4)+np.outer(v5, v5)
        for key, value in {'fitopt004_counterfactual': repaired, 'centered4_published': old4,
                           'centered4_counterfactual': new4, 'massstep_source_reconstruction': old_cov,
                           'covariance_contribution_change': new_cov-old_cov}.items():
            error = float(np.max(abs(value-saved[branch+'_'+key])))
            maximum_error = max(maximum_error, error)
            assert error < 1e-12, (branch, key, error)
        a = np.where(base['zHD'] > .2, 1/37, -1/42)
        assert np.array_equal(a, saved['a_high_minus_low'])
        mass_lc = np.loadtxt(record(RELEASE/f'{branch}_syst/RAISIN_massstep_lcparams_cosmosis.txt'))
        stat_lc = np.loadtxt(record(RELEASE/f'{branch}_syst/RAISIN_stat_lcparams_cosmosis.txt'))
        packed = np.loadtxt(record(RELEASE/f'{branch}_syst/RAISIN_massstep.covmat'))
        assert packed[0] == len(ids)
        published = packed[1:].reshape(len(ids), len(ids))
        published += np.diag(mass_lc[:, 5]**2-stat_lc[:, 5]**2)
        assert np.max(abs(published-saved[branch+'_massstep_published'])) < 1e-13
        # The protocol replaces the component in the supplied rounded export,
        # retaining its other contributions and printing residual exactly.
        repaired_export = published + new_cov-old_cov
        replacement_error = float(np.max(abs(repaired_export-saved[branch+'_massstep_counterfactual'])))
        assert replacement_error < 1e-12
        maximum_error = max(maximum_error, replacement_error)
        p = np.eye(79)-np.ones((79, 79))/79
        values = {'old_contrast_sd_mag': float(np.sqrt(a @ published @ a)),
                  'counterfactual_contrast_sd_mag': float(np.sqrt(a @ repaired_export @ a)),
                  'source_formula_old_contrast_sd_mag': float(np.sqrt(a @ old_cov @ a)),
                  'source_formula_counterfactual_contrast_sd_mag': float(np.sqrt(a @ new_cov @ a)),
                  'old_rank_after_intercept': int(np.linalg.matrix_rank(p @ np.column_stack([old4, v5]), tol=1e-9)),
                  'new_rank_after_intercept': int(np.linalg.matrix_rank(p @ np.column_stack([new4, v5]), tol=1e-9)),
                  'high_low_response_to_application_repair_mag': float(a @ (repaired-four['DLMAG'])),
                  'fitted_step_retained_mag': float(gamma),
                  'maximum_published_covariance_rounding_difference': float(np.max(abs(old_cov-published)))}
        assert abs(values['old_contrast_sd_mag']-claimed['branches'][branch]['old_mass_group_contrast_sd_mag']) < 1e-12
        assert abs(values['counterfactual_contrast_sd_mag']-claimed['branches'][branch]['counterfactual_mass_group_contrast_sd_mag']) < 1e-12
        records[branch] = values
    result = {'all_gates_pass': True, 'max_saved_array_error': maximum_error, 'branches': records,
              'input_sha256': INPUTS,
              'scope': 'Source-backed fixed-amplitude FITOPT004 application repair and mass-group covariance only. Published nominal distances unchanged. No full covariance update or cosmology.'}
    (OUT/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ('all_gates_pass', 'max_saved_array_error', 'branches')}, indent=2))


if __name__ == '__main__':
    main()
