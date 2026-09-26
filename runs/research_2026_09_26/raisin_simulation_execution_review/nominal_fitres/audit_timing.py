#!/usr/bin/env python3
"""Descriptive timing audit of the two pinned RAISIN nominal simulation FITRES."""
import collections
import decimal
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
TRUTH_FIELDS = ('SIM_PKMJD', 'SIM_LIBID', 'SIM_TYPE_INDEX', 'SIM_STRETCH', 'SIM_AV')
PROTOCOL_SHA = '3697ece51ac4ce3c7ff0cccbb79f92c81b3d1cb25dbeec52363647b6e30d0838'


def parse(path):
    headers, names, rows = [], None, []
    with gzip.open(path, 'rt', errors='strict') as stream:
        for line_no, line in enumerate(stream, 1):
            if line.startswith('#'):
                headers.append(line.rstrip())
            elif line.startswith('VARNAMES:'):
                names = line.split()[1:]
            elif line.startswith('SN:'):
                assert names is not None
                values = line.split()[1:]
                assert len(values) == len(names), (path, line_no)
                row = dict(zip(names, values))
                row['_line'] = line_no
                rows.append(row)
    assert names and len(set(names)) == len(names)
    required = (*TRUTH_FIELDS, 'CID', 'PKMJDINI', 'PKMJD', 'PKMJDERR',
                'STRETCHERR', 'AVERR', 'ERRFLAG_FIT', 'CUTFLAG_SNANA')
    assert all(x in names for x in required)
    counts = collections.Counter(r['CID'] for r in rows)
    assert max(counts.values(), default=0) == 1, (path, [k for k, v in counts.items() if v > 1][:10])
    return {'headers': headers, 'names': names, 'rows': rows,
            'by_cid': {r['CID']: r for r in rows}}


def number(value):
    try:
        n = decimal.Decimal(value)
    except decimal.InvalidOperation:
        return None
    return n if n.is_finite() else None


def valid_peak(value):
    n = number(value)
    return n if n is not None and -8 < n < 999999 else None


def precision_counts(rows, field):
    c = collections.Counter()
    for r in rows:
        value = r[field]
        c['scientific' if 'e' in value.lower() else 'fixed'] += 1
        c['fraction_digits:' + str(len(value.split('.', 1)[1]) if '.' in value else 0)] += 1
    return dict(sorted(c.items()))


def summarize(values):
    arr = np.asarray([float(x) for x in values], dtype=float)
    if not len(arr):
        return {'n': 0}
    med = float(np.median(arr))
    return {'n': len(arr), 'mean': float(arr.mean()), 'sd_population': float(arr.std()),
            'median': med, 'mad': float(np.median(np.abs(arr - med))),
            'min': float(arr.min()), 'max': float(arr.max()),
            **{f'q{int(q * 100):02d}': float(np.quantile(arr, q)) for q in (.01, .05, .95, .99)},
            **{f'abs_le_{x:.2f}': int(np.count_nonzero(np.abs(arr) <= x + 1e-12))
               for x in (.01, .03, .05)}}


def strata(rows):
    return {'all_valid': rows,
            'errflag_zero': [r for r in rows if r['ERRFLAG_FIT'] == '0'],
            'errflag_nonzero': [r for r in rows if r['ERRFLAG_FIT'] != '0'],
            'cutflag_zero': [r for r in rows if r['CUTFLAG_SNANA'] == '0'],
            'cutflag_nonzero': [r for r in rows if r['CUTFLAG_SNANA'] != '0']}


def error_audit(rows, field):
    counts = collections.Counter()
    by_flag = collections.defaultdict(collections.Counter)
    printed = collections.Counter()
    for r in rows:
        raw = r[field]
        n = number(raw)
        key = ('nonnumeric_or_nonfinite' if n is None else
               'sentinel' if n <= -8 or n >= 999999 else
               'negative' if n < 0 else 'zero' if n == 0 else 'positive')
        counts[key] += 1
        by_flag[r['ERRFLAG_FIT']][key] += 1
        printed[raw] += 1
    return {'categories': dict(counts), 'by_errflag': {k: dict(v) for k, v in by_flag.items()},
            'unique_printed_values': dict(printed) if len(printed) <= 12 else None,
            'n_unique_printed_values': len(printed)}


def main():
    manifest = json.loads((ROOT / 'acquisition.json').read_text())
    protocol_hash = hashlib.sha256((ROOT / 'timing-audit-protocol.json').read_bytes()).hexdigest()
    assert protocol_hash == PROTOCOL_SHA
    data = {}
    for arm in ('nir', 'optnir'):
        path = ROOT / f'{arm}.FITRES.gz'
        entry = next(f for f in manifest['files'] if f['local_path'].endswith(f'/{arm}.FITRES.gz'))
        assert hashlib.sha256(path.read_bytes()).hexdigest() == entry['sha256']
        data[arm] = parse(path)
    out = {'protocol_sha256': protocol_hash, 'acquisition_sha256':
           hashlib.sha256((ROOT / 'acquisition.json').read_bytes()).hexdigest(),
           'arms': {}, 'common_cid': {}}
    for arm, d in data.items():
        rows = d['rows']
        result = {'rows': len(rows), 'headers': d['headers'], 'columns': d['names'],
                  'errflag_counts': dict(collections.Counter(r['ERRFLAG_FIT'] for r in rows)),
                  'cutflag_counts': dict(collections.Counter(r['CUTFLAG_SNANA'] for r in rows)),
                  'peak_validity': {}, 'precision': {}, 'errors': {}, 'differences_days': {}}
        for field in ('PKMJDINI', 'PKMJD', 'SIM_PKMJD'):
            result['peak_validity'][field] = dict(collections.Counter(
                'valid' if valid_peak(r[field]) is not None else 'invalid_or_sentinel' for r in rows))
            result['precision'][field] = precision_counts(rows, field)
        for field in ('PKMJDERR', 'STRETCHERR', 'AVERR'):
            result['errors'][field] = error_audit(rows, field)
        for label, first, second in [('initializer_minus_truth', 'PKMJDINI', 'SIM_PKMJD'),
                                      ('fitted_minus_truth', 'PKMJD', 'SIM_PKMJD'),
                                      ('fitted_minus_initializer', 'PKMJD', 'PKMJDINI')]:
            candidates = [r for r in rows if valid_peak(r[first]) is not None and valid_peak(r[second]) is not None]
            result['differences_days'][label] = {
                s: summarize([valid_peak(r[first]) - valid_peak(r[second]) for r in rr])
                for s, rr in strata(candidates).items()}
        out['arms'][arm] = result
    a, b = data['nir']['by_cid'], data['optnir']['by_cid']
    common = sorted(set(a) & set(b))
    out['common_cid']['n'] = len(common)
    out['common_cid']['nir_only'] = len(set(a) - set(b))
    out['common_cid']['optnir_only'] = len(set(b) - set(a))
    out['common_cid']['nir_only_cids'] = sorted(set(a) - set(b), key=int)
    out['common_cid']['optnir_only_cids'] = sorted(set(b) - set(a), key=int)
    out['common_cid']['identity_discrepancies'] = {}
    for field in TRUTH_FIELDS:
        bad = []
        for cid in common:
            aa, bb = number(a[cid][field]), number(b[cid][field])
            if aa is None or bb is None or aa != bb:
                bad.append({'CID': cid, 'nir': a[cid][field], 'optnir': b[cid][field]})
        out['common_cid']['identity_discrepancies'][field] = bad
    for label, field in [('fitted_optnir_minus_nir', 'PKMJD'),
                         ('initializer_optnir_minus_nir', 'PKMJDINI')]:
        vals = [valid_peak(b[cid][field]) - valid_peak(a[cid][field]) for cid in common
                if valid_peak(a[cid][field]) is not None and valid_peak(b[cid][field]) is not None]
        out['common_cid'][label] = summarize(vals)
    (ROOT / 'timing-result.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'rows': {k: v['rows'] for k, v in out['arms'].items()},
                      'common': out['common_cid']['n'],
                      'truth_discrepancies': {k: len(v) for k, v in out['common_cid']['identity_discrepancies'].items()},
                      'nir_initial': out['arms']['nir']['differences_days']['initializer_minus_truth']['all_valid'],
                      'nir_fit': out['arms']['nir']['differences_days']['fitted_minus_truth']['all_valid']}, indent=2))


if __name__ == '__main__':
    main()
