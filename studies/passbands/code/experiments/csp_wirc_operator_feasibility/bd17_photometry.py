#!/usr/bin/env python3
"""Frozen source-formula integration for CSP WIRC/RC2 BD17 feasibility only."""
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
P = json.loads((OUT / 'bd17-photometry-protocol.json').read_text())
for rec in P['source_files'].values():
    path = ROOT / rec['path']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == rec['sha256'], path
sed = np.loadtxt(ROOT / P['source_files']['bd17_sed']['path'])
assert np.all(np.diff(sed[:, 0]) > 0)
hc = float(P['units']['hc'].split()[0])
result = {'protocol_sha256': hashlib.sha256((OUT / 'bd17-photometry-protocol.json').read_bytes()).hexdigest(), 'bands': {}}
for name, key in [('wirc_j', 'wirc_j'), ('retrocam_rc2_j', 'retrocam_rc2_j')]:
    band = np.loadtxt(ROOT / P['source_files'][key]['path'])
    x, s = band.T
    assert np.all(np.diff(x) > 0) and x[0] >= sed[0, 0] and x[-1] <= sed[-1, 0]
    f_on_filter = np.interp(x, sed[:, 0], sed[:, 1])
    count_filter = float(np.trapezoid(s * f_on_filter * x, x) / hc)
    union = np.unique(np.r_[x, sed[(sed[:, 0] > x[0]) & (sed[:, 0] < x[-1]), 0]])
    count_union = float(np.trapezoid(np.interp(union, x, s) * np.interp(union, sed[:, 0], sed[:, 1]) * union, union) / hc)
    zp = P['source_zero_points'][key]
    m_filter = float(-2.5 * np.log10(count_filter) + zp)
    m_union = float(-2.5 * np.log10(count_union) + zp)
    result['bands'][name] = {'n_filter_rows': len(x), 'count_filter_grid': count_filter, 'count_union_grid': count_union, 'mag_filter_grid': m_filter, 'mag_union_grid': m_union, 'quadrature_diff_mag': m_union-m_filter, 'wavelength_range_angstrom': [float(x[0]), float(x[-1])]}
result['released_rc2_bd17_magref'] = 8.4192
result['rc2_difference_from_released_mag'] = result['bands']['retrocam_rc2_j']['mag_filter_grid'] - 8.4192
result['wirc_minus_rc2_mag'] = result['bands']['wirc_j']['mag_filter_grid'] - result['bands']['retrocam_rc2_j']['mag_filter_grid']
(OUT / 'bd17-photometry-result.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
print(json.dumps(result, indent=2))
