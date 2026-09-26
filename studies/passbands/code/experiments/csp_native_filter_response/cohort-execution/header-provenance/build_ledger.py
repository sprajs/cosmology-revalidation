#!/usr/bin/env python3
"""Read-only, source-field CSP peak provenance ledger for the frozen 42 IDs."""
import csv
import gzip
import hashlib
import json
import pathlib
import struct
from decimal import Decimal

ROOT = pathlib.Path(__file__).resolve().parents[5]
HERE = pathlib.Path(__file__).resolve().parent
PREP = HERE.parents[1] / 'cohort-preparation'
EXEC = HERE.parent
ARCHIVE = ROOT / 'runs/research_2026_09_26/raisin_profile_solver_review/csp-filter-interpretation/native-design/CSP-nominal.FITRES.gz'


def read_fit(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        lines = f.readlines()
    names = next(line.split()[1:] for line in lines if line.startswith('VARNAMES:'))
    return {v[0]: dict(zip(names, v)) for line in lines if line.startswith('SN:') for v in [line.split()[1:]]}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def header(path):
    out = {}
    for line in path.read_text().splitlines():
        if line.startswith('OBS:'):
            break
        if ':' in line and not line.startswith('#'):
            k, v = line.split(':', 1)
            out[k.strip()] = v.split('#')[0].strip()
    return out


def f32(v):
    return struct.unpack('<f', struct.pack('<f', float(v)))[0]


def main():
    nominal = PREP / 'full/data/nominal/CSPDR3_RAISIN'
    ids = [x.strip().replace('CSPDR3_', '').replace('.PKMJD.DAT', '') for x in (nominal / 'CSPDR3_RAISIN.LIST').read_text().splitlines() if x.strip() and not x.startswith('#')]
    files = {
        'archived_nir': ARCHIVE,
        'replay_nir': EXEC / 'full/fits/nominal/fit.FITRES.TEXT',
        'oldpk_nir': HERE / 'author/output/fit_nir_sys_oldpkmjds/CSP_RAISIN/CSPDR3_RAISIN/FITOPT000.FITRES.gz',
        'archived_optical': HERE / 'author/output/fit_optical_sys/CSP_RAISIN/CSPDR3_RAISIN/FITOPT000.FITRES.gz',
        'archived_opticalnir': HERE / 'author/output/fit_opticalnir_sys/CSP_RAISIN/CSPDR3_RAISIN/FITOPT000.FITRES.gz',
        'older_opticalnir': HERE / 'author/CSPDR3_RAISIN_optnir.FITRES.TEXT',
    }
    fit = {k: read_fit(v) for k, v in files.items()}
    fields = ('zHEL', 'zCMB', 'zHD', 'VPEC', 'MWEBV', 'PKMJDINI', 'PKMJD', 'STRETCH', 'AV', 'RV')
    rows = []
    for cid in ids:
        hpath = nominal / f'CSPDR3_{cid}.PKMJD.DAT'
        h = header(hpath)
        row = {'CID': cid, 'release_file': str(hpath.relative_to(ROOT)), 'release_sha256': digest(hpath),
               'release_PEAKMJD': h.get('PEAKMJD'), 'release_REDSHIFT_HELIO': h.get('REDSHIFT_HELIO'),
               'release_REDSHIFT_FINAL': h.get('REDSHIFT_FINAL'), 'release_MWEBV': h.get('MWEBV')}
        for label in files:
            rec = fit[label].get(cid, {})
            for field in fields:
                row[f'{label}_{field}'] = rec.get(field, '')
        row['archive_vs_release_peak_difference_day'] = str(Decimal(row['archived_nir_PKMJDINI']) - Decimal(row['release_PEAKMJD']))
        row['archive_vs_release_peak_float32_ulps'] = ((f32(row['archived_nir_PKMJDINI']) - f32(row['release_PEAKMJD'])) / (2**-8))
        rows.append(row)
    out = HERE / 'full42-header-fit-ledger.csv'
    with out.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    summary = {'n_ids': len(ids), 'fit_files': {k: {'path': str(v.relative_to(ROOT)), 'sha256': digest(v), 'n_rows': len(fit[k]), 'cohort_present': sum(cid in fit[k] for cid in ids)} for k,v in files.items()},
               'release_peak_vs_archived_nir_printed_different': [r['CID'] for r in rows if r['release_PEAKMJD'] != r['archived_nir_PKMJDINI']],
               'release_peak_vs_archived_nir_float32_different': [r['CID'] for r in rows if f32(r['release_PEAKMJD']) != f32(r['archived_nir_PKMJDINI'])],
               'replay_peak_vs_release_printed_different': [r['CID'] for r in rows if r['release_PEAKMJD'] != r['replay_nir_PKMJDINI']],
               'replay_peak_vs_release_float32_different': [r['CID'] for r in rows if f32(r['release_PEAKMJD']) != f32(r['replay_nir_PKMJDINI'])],
               'other_original_fields_different_archive_vs_replay': {field: [r['CID'] for r in rows if r[f'archived_nir_{field}'] != r[f'replay_nir_{field}']] for field in fields},
               'ledger_sha256': digest(out)}
    (HERE / 'full42-header-fit-summary.json').write_text(json.dumps(summary, indent=2) + '\n')


if __name__ == '__main__':
    main()
