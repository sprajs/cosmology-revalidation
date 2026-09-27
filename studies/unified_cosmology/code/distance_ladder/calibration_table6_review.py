"""Compare explicit published SN-free Cepheid errors with sibling covariances.

No covariance component is inferred from subtraction; the comparison cannot
establish which unreleased construction step produced the published matrix.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
WORK = ROOT / '.work/unified-cosmology/calibration-input-search'
PDF = ROOT / '.work/unified-cosmology/distance-ladder/Riess2022v3.pdf'
FACTOR = ROOT / 'studies/unified_cosmology/results/distance_ladder/cepheid-factor.json'
DAT = ROOT / '.work/unified-cosmology/calibration-interface/Pantheon+SH0ES.dat'
COV = ROOT / '.work/unified-cosmology/calibration-interface/Pantheon+SH0ES_STATONLY.cov'
PINS = {PDF: '572e5d2c0f32bf6e80a118eec13ed0107e1e2830f9ebe90b5bc4e9d1b7fa5df1',
        DAT: '1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8',
        COV: '9f177129a332735d3637affd20054080d5260815f3ca0809120c05b2c902297f'}
HOSTS = ['N1448', 'N3147', 'N5468', 'N5643']
ALIASES = {'2008fv': '2008fv_comb'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT / 'studies/unified_cosmology/results/distance_ladder/calibration-table6-review.json')
    args = parser.parse_args()
    for path, digest in PINS.items():
        assert sha(path) == digest
    before = {str(p.relative_to(ROOT)): sha(p) for p in [*PINS, FACTOR]}
    WORK.mkdir(parents=True, exist_ok=True)
    text = subprocess.check_output(['pdftotext', '-f', '30', '-l', '30', '-layout', str(PDF), '-'], text=True)
    assert 'Table 6.' in text and 'without inclusion of any SN in any host' in text
    text_path = WORK / 'Riess2022v3-table6.txt'
    text_path.write_text(text)
    rows = []
    for line in text.splitlines():
        fields = line.split()
        if len(fields) == 12 and fields[0].isdigit() and fields[1] in HOSTS:
            rows.append({'table_row': int(fields[0]), 'host': fields[1], 'SN': fields[2],
                         'mu_b': float(fields[9]), 'sigma_b': float(fields[10])})
    assert len(rows) == 9
    factor = json.loads(FACTOR.read_text())
    H = np.asarray(factor['covariance_distance_modulus'])
    data = pd.read_csv(DAT, sep=r'\s+')
    with COV.open() as stream:
        n = int(stream.readline())
        cov = np.loadtxt(stream).reshape(n, n)
    assert n == len(data) == 1701
    results = []
    for host in HOSTS:
        printed = [r for r in rows if r['host'] == host]
        sigma = printed[0]['sigma_b']
        assert all(r['sigma_b'] == sigma and r['mu_b'] == printed[0]['mu_b'] for r in printed)
        ci = factor['host_order'].index(host)
        variance = float(H[ci, ci])
        events = [ALIASES.get(r['SN'], r['SN']) for r in printed]
        pairs = []
        for i, sn1 in enumerate(events):
            index1 = np.flatnonzero(data['CID'].eq(sn1).to_numpy() & data['IS_CALIBRATOR'].eq(1).to_numpy())
            assert len(index1)
            for sn2 in events[i+1:]:
                index2 = np.flatnonzero(data['CID'].eq(sn2).to_numpy() & data['IS_CALIBRATOR'].eq(1).to_numpy())
                assert len(index2)
                block = cov[np.ix_(index1, index2)]
                pairs.append({'SN1': sn1, 'SN2': sn2, 'original_rows1': index1.tolist(),
                              'original_rows2': index2.tolist(), 'covariance_min': float(block.min()),
                              'covariance_max': float(block.max()), 'values_mag2': np.unique(block).tolist()})
        distinct = sorted(set(v for p in pairs for v in p['values_mag2']))
        assert len(distinct) == 1
        sibling_cov = distinct[0]
        lo, hi = sigma - .0005, sigma + .0005
        results.append({'host': host, 'table6b_rows': printed,
                        'printed_sigma_rounding_interval_mag': [lo, hi],
                        'printed_sigma_squared_mag2': sigma**2,
                        'printed_variance_rounding_interval_mag2': [lo**2, hi**2],
                        'SN_free_factor_sigma_mag': variance**.5,
                        'SN_free_factor_variance_mag2': variance,
                        'factor_sigma_within_printed_rounding_interval': bool(lo <= variance**.5 <= hi),
                        'distinct_sibling_STATONLY_covariance_mag2': sibling_cov,
                        'sibling_covariance_over_printed_variance': sibling_cov / sigma**2,
                        'ratio_including_printed_rounding': [sibling_cov/hi**2, sibling_cov/lo**2],
                        'sibling_covariance_over_factor_variance': sibling_cov/variance,
                        'sqrt_sibling_covariance_mag_not_host_error_measurement': sibling_cov**.5,
                        'sibling_covariance_in_printed_variance_interval': bool(lo**2 <= sibling_cov <= hi**2),
                        'physical_SN_pair_ledger': pairs})
    assert all(r['factor_sigma_within_printed_rounding_interval'] and not r['sibling_covariance_in_printed_variance_interval'] for r in results)
    for name, digest in before.items():
        assert sha(ROOT/name) == digest
    result = {'status': 'published_table6b_errors_match_SN_free_factor_not_sibling_covariance_scale',
              'paper': {'url': 'https://arxiv.org/pdf/2112.04510v3', 'PDF_page_one_based': 30,
                        'printed_page': 30, 'table': '6', 'column': 'mu_host superscript b and adjacent sigma',
                        'definition': 'Cepheid-based distance excluding every SN in every host',
                        'visual_review': 'Page30 rendered and inspected; column/footnote verified against pdftotext before this audit.',
                        'text_path': str(text_path.relative_to(ROOT)), 'text_sha256': sha(text_path)},
              'host_comparisons': results,
              'released_table_CEPH_columns': [c for c in data.columns if 'CEPH' in c.upper()],
              'released_table_has_CEPH_DISTERR': 'CEPH_DISTERR' in data.columns,
              'interpretation': ['The factor-of-about-two variance scale is absent from the explicit Table6(b) host errors, even allowing their printed rounding.',
                                 'The released STATONLY entries are complete SN-pair covariances, not certified pure Cepheid components. This comparison localizes a difference between products but does not identify its cause.',
                                 'No factor2 correction, host-component substitution, covariance rescaling or new cosmology follows from this test. The table calls its distance parameters approximations; exact production covariance ancestry is still missing.'],
              'input_sha256': before,
              'source_sha256': {str(Path(__file__).relative_to(ROOT)): sha(__file__)},
              'cosmology_fits': 0, 'covariance_changes': 0}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': result['status'], 'comparisons': results}, indent=2))


if __name__ == '__main__':
    main()
