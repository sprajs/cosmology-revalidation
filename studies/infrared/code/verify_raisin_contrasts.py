"""Independent raw-table reconstruction of frozen RAISIN paired statistics."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RELEASE = ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/distances/w'
RUN = ROOT/'runs/research_2026_09_26/raisin_differential'
OUT = ROOT/'runs/research_2026_09_26/raisin_independent_verification'
FILES = {}


def record(p):
    FILES[str(p.relative_to(ROOT))] = hashlib.sha256(p.read_bytes()).hexdigest()
    return p


def main():
    OUT.mkdir(exist_ok=False)
    record(Path(__file__))
    members = list(csv.DictReader(record(RUN/'frozen-membership.csv').open()))
    ids = [r['CID'] for r in members]
    tabs = {}; cov = {}; errors = {}
    export_checks = {}
    for branch in ('nir', 'optical', 'opticalnir'):
        path = record(RELEASE/f'{branch}_dist/RAISIN_combined_FITOPT000.FITRES')
        lines = path.read_text().splitlines()
        fields = next(line.split()[1:] for line in lines if line.startswith('VARNAMES:'))
        rows = [dict(zip(fields, line.split()[1:])) for line in lines if line.startswith('SN:')]
        assert [r['CID'] for r in rows] == ids and len(set(ids)) == 79
        tab = {key: np.array([float(r[key]) for r in rows]) for key in
               ('DLMAG', 'DLMAG_biascor', 'MASS_CORR', 'DLMAGERR', 'zHD', 'HOST_LOGMASS')}
        tabs[branch] = tab; errors[branch] = tab['DLMAGERR']
        stat = np.loadtxt(record(RELEASE/f'{branch}_syst/RAISIN_stat_lcparams_cosmosis.txt'))
        total = np.loadtxt(record(RELEASE/f'{branch}_syst/RAISIN_all_lcparams_cosmosis.txt'))
        packed = np.loadtxt(record(RELEASE/f'{branch}_syst/RAISIN_all.covmat'))
        assert packed[0] == 79 and len(packed)-1 == 79**2
        off = packed[1:].reshape(79, 79)
        assert np.max(abs(off-off.T)) == 0 and np.max(abs(off.diagonal())) == 0
        cov[branch] = off + np.diag(total[:, 5]**2-stat[:, 5]**2)
        mb_error = max(abs(stat[:, 4]-(tab['DLMAG']-19.36)))
        sigma_error = max(abs(stat[:, 5]-tab['DLMAGERR']))
        assert mb_error <= 5.1e-7 and sigma_error <= 5.1e-7
        assert max(abs(stat[:, 1]-tab['zHD'])) <= 5.1e-7
        export_checks[branch] = {'mb_rounding_max': float(mb_error), 'sigma_rounding_max': float(sigma_error),
                                 'systematic_min_eigenvalue_preserved': float(np.linalg.eigvalsh(cov[branch])[0])}
    expected = json.loads(record(RUN/'results.json').read_text())['contrasts']
    comparisons = []; max_error = 0.
    for row in expected:
        branch = row['branch_minus_nir']; variant = row['correction_variant']; split = row['split']
        t = tabs[branch]; n = tabs['nir']
        valid = np.ones(79, dtype=bool); high = t['zHD'] > .2; low = t['zHD'] < .1
        if split.endswith('_high_CSP_low'):
            survey = {'PS1': 'PS1MD', 'DES': 'DES'}[split.split('_')[0]]
            valid = np.array([r['survey'] == survey or r['stratum'] == 'low' for r in members])
        elif split == 'host_ge10': valid = t['HOST_LOGMASS'] >= 10
        elif split == 'host_lt10': valid = t['HOST_LOGMASS'] < 10
        high &= valid; low &= valid
        assert high.sum() == row['high_n'] and low.sum() == row['low_n']
        weights = high/high.sum() - low/low.sum()
        delta = t['DLMAG']-n['DLMAG']
        if variant in ('without_bias', 'without_bias_or_mass'):
            delta = delta + t['DLMAG_biascor']-n['DLMAG_biascor']
        if variant in ('without_mass', 'without_bias_or_mass'):
            delta = delta - t['MASS_CORR']+n['MASS_CORR']
        calculated = {'contrast_mag': float(weights @ delta),
                      'empirical_paired_se_mag': float(np.sqrt(np.var(delta[high], ddof=1)/high.sum()+np.var(delta[low], ddof=1)/low.sum())),
                      'statistical_rho0_se_mag': float(np.sqrt(weights**2 @ (errors[branch]**2+errors['nir']**2))),
                      'statistical_per_object_rho_bounds_se_mag':
                          [float(np.linalg.norm(weights*(errors[branch]-errors['nir']))),
                           float(np.linalg.norm(weights*(errors[branch]+errors['nir'])))]}
        if variant == 'released':
            vb = float(weights @ cov[branch] @ weights); vn = float(weights @ cov['nir'] @ weights)
            assert vb > 0 and vn > 0
            sb, sn = np.sqrt(vb), np.sqrt(vn)
            tb = np.sqrt(vb + weights**2 @ errors[branch]**2)
            tn = np.sqrt(vn + weights**2 @ errors['nir']**2)
            calculated.update(systematic_marginal_se_mag=[float(sb), float(sn)],
                              systematic_unknown_cross_covariance_se_bounds_mag=[float(abs(sb-sn)), float(sb+sn)],
                              total_branchwise_stat_plus_sys_unknown_cross_se_bounds_mag=[float(abs(tb-tn)), float(tb+tn)])
        error = max(float(np.max(abs(np.asarray(value)-np.asarray(row[key])))) for key, value in calculated.items())
        assert error < 1e-12, (branch, variant, split, error)
        max_error = max(max_error, error)
        comparisons.append({'branch': branch, 'variant': variant, 'split': split, 'max_error': error})
    assert len(comparisons) == 40
    result = {'checks': 40, 'maximum_arithmetic_error': max_error, 'export_checks': export_checks,
              'comparisons': comparisons, 'input_sha256': FILES,
              'scope': 'Independent raw-table arithmetic and conditional marginal-variance bounds. Tiny export-negative eigenvalues preserved; no joint covariance, physical correction, or coverage established.'}
    (OUT/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ('checks', 'maximum_arithmetic_error', 'export_checks')}, indent=2))


if __name__ == '__main__':
    main()
