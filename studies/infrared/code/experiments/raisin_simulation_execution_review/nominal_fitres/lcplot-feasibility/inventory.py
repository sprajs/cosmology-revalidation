#!/usr/bin/env python3
"""Schema and coverage inventory only; no flux or timing differences."""
import collections
import gzip
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
acq = json.loads((HERE / 'acquisition.json').read_text())
assert hashlib.sha256((HERE / 'FITOPT000.LCPLOT.gz').read_bytes()).hexdigest() == acq['sha256']

with gzip.open(PARENT / 'nir.FITRES.gz', 'rt') as f:
    for line in f:
        if line.startswith('VARNAMES:'):
            fit_names = line.split()[1:]
            break
    fit_cids = [dict(zip(fit_names, line.split()[1:]))['CID']
                for line in f if line.startswith('SN:')]

row_markers = collections.Counter()
flags, bands, ifits = collections.Counter(), collections.Counter(), collections.Counter()
flag_band = collections.Counter()
per_cid_all, per_cid_data = collections.Counter(), collections.Counter()
keys = collections.Counter()
precisions = {name: collections.Counter() for name in ('MJD', 'Tobs', 'SIM_FLUXCAL', 'FLUXCAL', 'FLUXCAL_ERR')}
plot_cids, seen = [], set()
with gzip.open(HERE / 'FITOPT000.LCPLOT.gz', 'rt') as f:
    first = f.readline()
    assert first.startswith('VARNAMES:')
    names = first.split()[1:]
    assert len(set(names)) == len(names)
    for lineno, line in enumerate(f, 2):
        parts = line.split()
        marker = parts[0] if parts else 'blank'
        row_markers[marker] += 1
        if marker != 'OBS:':
            continue
        assert len(parts) - 1 == len(names), lineno
        row = dict(zip(names, parts[1:]))
        cid = row['CID']
        if cid not in seen:
            plot_cids.append(cid)
            seen.add(cid)
        per_cid_all[cid] += 1
        if row['DATAFLAG'] == '1':
            per_cid_data[cid] += 1
        flags[row['DATAFLAG']] += 1
        bands[row['BAND']] += 1
        ifits[row['IFIT']] += 1
        flag_band[(row['DATAFLAG'], row['BAND'])] += 1
        keys[(cid, row['MJD'], row['BAND'])] += 1
        for field in precisions:
            value = row[field]
            mantissa = value.lower().split('e')[0]
            decimals = len(mantissa.split('.', 1)[1]) if '.' in mantissa else 0
            precisions[field][f'mantissa_decimals_{decimals}'] += 1
            precisions[field]['scientific' if 'e' in value.lower() else 'fixed'] += 1

out = {
    'source_sha256': acq['sha256'],
    'varnames': names,
    'markers': dict(row_markers),
    'n_plot_cids': len(plot_cids),
    'plot_cids_first_last': [plot_cids[0], plot_cids[-1]],
    'plot_cids_equal_fitres_first_500_in_order': plot_cids == fit_cids[:500],
    'plot_cids_equal_numeric_1_to_500': set(map(int, plot_cids)) == set(range(1, 501)),
    'n_fitres_cids': len(fit_cids),
    'fitres_cids_without_plot': len(set(fit_cids) - seen),
    'dataflag_counts': dict(flags),
    'band_counts': dict(bands),
    'ifit_counts': dict(ifits),
    'dataflag_band_counts': {f'{a}:{b}': n for (a, b), n in sorted(flag_band.items())},
    'rows_per_cid_all_min_max': [min(per_cid_all.values()), max(per_cid_all.values())],
    'rows_per_cid_data_min_max': [min(per_cid_data.values()), max(per_cid_data.values())],
    'cid_mjd_band_key_duplicate_rows': sum(n - 1 for n in keys.values() if n > 1),
    'precision': {k: dict(v) for k, v in precisions.items()},
}
(HERE / 'inventory.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps({k: out[k] for k in ('n_plot_cids','n_fitres_cids','dataflag_counts','band_counts','dataflag_band_counts','plot_cids_equal_fitres_first_500_in_order')},indent=2))
