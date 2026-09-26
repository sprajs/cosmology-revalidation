"""Registered descriptive signed-ancestor extraction checks; no cosmology fit."""
from collections import Counter, defaultdict
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'runs/research_2026_09_26/raisin_signed_baseline'
COHORT = ROOT / 'runs/research_2026_09_26/astra_design/raisin_signed_refit/cohort.csv'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def write_csv(path, rows):
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def freeze():
    OUT.mkdir(exist_ok=True)
    assert not (OUT / 'protocol.json').exists()
    cases = list(csv.DictReader(COHORT.open()))
    paths = [COHORT, Path(__file__)] + [ROOT / c['raw_path'] for c in cases]
    protocol = {
        'question': 'Adequacy of signed author-DIFFIMG errors and constant object-band baseline model, conditional on ten selected DES16 SNe.',
        'known_before_freeze': 'Positive-only RAISIN lineage, 33 negative header-fit-window rows, precise positive-error roundtrip, and native fitter instability are already known. This is a registered follow-up, not pristine blinding.',
        'cohort': [c['CID'] for c in cases],
        'phase': {'primary': 'MJD < released PEAKMJD - 180 observer days', 'sensitivity': 'MJD < released PEAKMJD - 365 observer days'},
        'filters': ['g', 'r', 'i', 'z'],
        'mask': 'No clipping or PHOTFLAG selection. Finite flux, positive finite quoted errors required; gate failures stop rather than silently delete rows.',
        'estimands': ['Per object-band raw inverse-variance baseline and conditional error; signed standardized median and percentiles; weighted-intercept projected squared residual per residual degree of freedom.',
                      'Within object-band all-pair flux differences divided by sqrt(sigma_i^2+sigma_j^2), by fixed lag bin and same integer MJD.',
                      'Projected residual pair products minus their expected diagonal-Gaussian intercept-projection covariance.'],
        'lag_edges_observer_days': [0, 0.5, 7, 30, 180, None],
        'interpretation': 'Descriptive only. Pairs and epochs are dependent. No population p-value, automatic error renormalization, or baseline subtraction from the SN fits. Template contamination/host subtraction/selection and flux-dependent error estimates remain alternatives.',
        'source_hashes': {str(p.relative_to(ROOT)): sha(p) for p in paths},
    }
    dump(OUT / 'protocol.json', protocol)
    (OUT / 'executed-source.py').write_bytes(Path(__file__).read_bytes())
    print(json.dumps({'protocol_sha256': sha(OUT / 'protocol.json')}))


def parse_raw(path):
    rows = []
    for lineno, line in enumerate(path.read_text().splitlines(), 1):
        words = line.split()
        if not words:
            continue
        if words[0] == 'VARLIST:':
            names = words[1:]
        elif words[0] == 'OBS:':
            assert len(words[1:]) == len(names)
            values = dict(zip(names, words[1:]))
            band = values.get('FLT', values.get('BAND'))
            if band not in 'griz':
                continue
            row = {'MJD': float(values['MJD']), 'flux': float(values['FLUXCAL']),
                   'error': float(values['FLUXCALERR']), 'band': band,
                   'PHOTFLAG': values.get('PHOTFLAG', ''), 'source_line': lineno}
            assert all(math.isfinite(row[k]) for k in ('MJD', 'flux', 'error')) and row['error'] > 0
            rows.append(row)
    return rows


def score():
    assert not (OUT / 'result.json').exists(), 'Preserve completed result'
    protocol = json.loads((OUT / 'protocol.json').read_text())
    for name, expected in protocol['source_hashes'].items():
        assert sha(ROOT / name) == expected, name
    cases = list(csv.DictReader(COHORT.open()))
    groups, pairs, ledger, flags = [], [], [], []
    lags = [(0, .5), (.5, 7), (7, 30), (30, 180), (180, np.inf)]
    for cut in (180, 365):
        for case in cases:
            cid = case['CID']; peak = float(case['peak_header'])
            data = [r for r in parse_raw(ROOT / case['raw_path']) if r['MJD'] < peak - cut]
            for r in data:
                ledger.append({'cut_days': cut, 'CID': cid, **r, 'source_path': case['raw_path']})
            for band in protocol['filters']:
                rows = [r for r in data if r['band'] == band]
                assert len(rows) > 2, (cut, cid, band)
                f = np.array([r['flux'] for r in rows]); s = np.array([r['error'] for r in rows])
                t = np.array([r['MJD'] for r in rows]); v = 1 / s; x = f / s
                norm = v @ v; mean = (v @ x) / norm; se = 1 / np.sqrt(norm)
                residual = x - v * mean
                Q = residual @ residual
                groups.append({'cut_days': cut, 'CID': cid, 'band': band, 'N': len(rows),
                               'negative': int(sum(f < 0)), 'baseline_FLUXCAL': float(mean),
                               'baseline_conditional_diag_SE': float(se),
                               'baseline_conditional_diag_z': float(mean / se),
                               'Q_after_intercept': float(Q), 'dof': len(rows) - 1,
                               'Q_per_dof': float(Q / (len(rows) - 1)),
                               'median_flux_over_error': float(np.median(x)),
                               'p05_flux_over_error': float(np.quantile(x, .05)),
                               'p95_flux_over_error': float(np.quantile(x, .95)),
                               'median_error': float(np.median(s))})
                for flag, n in sorted(Counter(r['PHOTFLAG'] for r in rows).items()):
                    chosen = np.array([r['PHOTFLAG'] == flag for r in rows])
                    flags.append({'cut_days': cut, 'CID': cid, 'band': band, 'PHOTFLAG': flag,
                                  'N': n, 'negative': int(np.sum(f[chosen] < 0)),
                                  'median_flux_over_error': float(np.median(x[chosen])),
                                  'mean_projected_residual_squared': float(np.mean(residual[chosen] ** 2))})
                i, j = np.triu_indices(len(rows), 1)
                lag = np.abs(t[i] - t[j]); same_night = np.floor(t[i]) == np.floor(t[j])
                diff = (f[i] - f[j]) / np.sqrt(s[i] ** 2 + s[j] ** 2)
                # P_ij=-v_i*v_j/(v'v) for distinct rows under diagonal Gaussian errors.
                product_excess = residual[i] * residual[j] + v[i] * v[j] / norm
                for lower, upper in lags:
                    for nightly in (False, True):
                        take = (lag >= lower) & (lag < upper) & (same_night == nightly)
                        n = int(np.sum(take))
                        if n:
                            pairs.append({'cut_days': cut, 'CID': cid, 'band': band,
                                          'lag_lower': lower, 'lag_upper': None if np.isinf(upper) else upper,
                                          'same_integer_MJD': nightly, 'N_pairs': n,
                                          'mean_normalized_difference_squared': float(np.mean(diff[take] ** 2)),
                                          'mean_projected_product_excess': float(np.mean(product_excess[take]))})
    write_csv(OUT / 'baseline-groups.csv', groups)
    write_csv(OUT / 'pair-summaries.csv', pairs)
    write_csv(OUT / 'baseline-source-rows.csv', ledger)
    write_csv(OUT / 'flag-strata.csv', flags)
    summary = []
    for cut in (180, 365):
        gs = [g for g in groups if g['cut_days'] == cut]
        q = np.array([g['Q_per_dof'] for g in gs]); means = np.array([g['baseline_FLUXCAL'] for g in gs])
        summary.append({'cut_days': cut, 'N_epochs': sum(g['N'] for g in gs), 'N_groups': len(gs),
                        'N_negative': sum(g['negative'] for g in gs),
                        'pooled_Q_per_dof': sum(g['Q_after_intercept'] for g in gs) / sum(g['dof'] for g in gs),
                        'group_Q_per_dof_median': float(np.median(q)), 'group_Q_per_dof_range': [float(q.min()), float(q.max())],
                        'baseline_FLUXCAL_median': float(np.median(means)),
                        'baseline_FLUXCAL_range': [float(means.min()), float(means.max())]})
    result = {'protocol_sha256': sha(OUT / 'protocol.json'), 'summary': summary,
              'no_global_p_value': True, 'no_fit_data_or_errors_modified': True,
              'output_sha256': {p.name: sha(p) for p in OUT.glob('*.csv')}}
    dump(OUT / 'result.json', result)
    print(json.dumps(result))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('action', choices=['freeze', 'score'])
    args = parser.parse_args()
    {'freeze': freeze, 'score': score}[args.action]()
