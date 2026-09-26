"""Summarize every historical-source native result without choosing fit minima."""
from collections import Counter, defaultdict
from pathlib import Path
import csv
import gzip
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'runs/research_2026_09_26/raisin_historical_native'
MODERN = ROOT / 'runs/research_2026_09_26/astra_design/raisin_signed_refit'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table(path, prefix):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as handle:
        lines = handle.read().splitlines()
    columns = next(line.split()[1:] for line in lines if line.startswith('VARNAMES:'))
    return [dict(zip(columns, line.split()[1:])) for line in lines if line.startswith(prefix + ':')]


def write_csv(path, rows):
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)


def main():
    assert not (OUT / 'summary.json').exists()
    protocol_path = OUT / 'fit-protocol-shortprefix.json'
    protocol = json.loads(protocol_path.read_text())
    author_path = MODERN / 'author-optical-FITOPT000.FITRES.gz'
    author = {r['CID']: r for r in table(author_path, 'SN')}
    hashes = {str(protocol_path.relative_to(ROOT)): sha(protocol_path), str(author_path.relative_to(ROOT)): sha(author_path)}
    fits, index, masks = [], {}, {}
    fields = ('DLMAG', 'AV', 'STRETCH', 'PKMJD', 'FITCHI2')
    for job in protocol['jobs']:
        if job['stage'] != 'full':
            continue
        directory = (ROOT / job['nml']).parent
        execution = json.loads((directory / 'execution.json').read_text())
        assert execution['returncode'] == 0 and execution['protocol_sha256'] == sha(protocol_path)
        assert execution['log_sha256'] == sha(directory / 'fit.log')
        rows = table(directory / 'fit.FITRES.TEXT', 'SN')
        epochs = table(directory / 'fit.LCPLOT.TEXT', 'OBS')
        assert sorted(r['CID'] for r in rows) == protocol['cohort']
        for path in directory.glob('*'):
            if path.is_file():
                hashes[str(path.relative_to(ROOT))] = sha(path)
        for row in rows:
            cid = row['CID']
            accepted = [e for e in epochs if e['CID'] == cid and int(e['DATAFLAG']) == 1]
            mask = Counter((e['BAND'], e['MJD'], e['FLUXCAL']) for e in accepted)
            nfree = 4 if job['timing'] == 'free' else 3
            assert len(accepted) == int(row['NDOF']) + nfree
            covariance = np.diag([float(row[k]) ** 2 for k in ('STRETCHERR', 'AVERR', 'DLMAGERR')])
            for i, j, k in ((0, 1, 'COV_STRETCH_AV'), (0, 2, 'COV_STRETCH_DLMAG'), (1, 2, 'COV_AV_DLMAG')):
                covariance[i, j] = covariance[j, i] = float(row[k])
            record = {'CID': cid, 'law': job['law'], 'timing': job['timing'], 'arm': job['arm'],
                      'start': job['start'], 'ERRFLAG_FIT': int(row['ERRFLAG_FIT']),
                      **{k: float(row[k]) for k in fields}, 'NDOF': int(row['NDOF']),
                      'RV': float(row['RV']), 'RVERR': float(row['RVERR']), 'PKMJDINI': float(row['PKMJDINI']),
                      'PKMJDERR': float(row['PKMJDERR']), 'accepted': len(accepted),
                      'accepted_negative': sum(float(e['FLUXCAL']) < 0 for e in accepted),
                      'covariance_min_eigenvalue': float(np.linalg.eigvalsh(covariance)[0])}
            key = (job['law'], job['timing'], job['arm'], job['start'], cid)
            fits.append(record); index[key] = record; masks[key] = mask
    write_csv(OUT / 'fit-ledger.csv', fits)
    cases = []
    for law in (94, 99):
        for timing in ('free', 'fixed'):
            for arm in ('R', 'A', 'B', 'H'):
                for cid in protocol['cohort']:
                    keys = [(law, timing, arm, start, cid) for start in (1.0, .85, 1.15)]
                    if not all(k in index for k in keys):
                        continue
                    rs = [index[k] for k in keys]
                    cases.append({'CID': cid, 'law': law, 'timing': timing, 'arm': arm,
                                  'same_accepted_mask': all(masks[k] == masks[keys[0]] for k in keys[1:]),
                                  **{f'{k}_range': max(r[k] for r in rs) - min(r[k] for r in rs) for k in fields},
                                  'all_ERRFLAG_zero': all(r['ERRFLAG_FIT'] == 0 for r in rs)})
    write_csv(OUT / 'start-spreads.csv', cases)
    baseline = []
    for cid in protocol['cohort']:
        record = {'CID': cid}
        for law in (94, 99):
            row = index[(law, 'free', 'R', 1.0, cid)]
            for field in fields:
                record[f'old{law}_minus_author_{field}'] = row[field] - float(author[cid][field])
            record[f'old{law}_minus_author_NDOF'] = row['NDOF'] - int(author[cid]['NDOF'])
        modern = {r['CID']: r for r in table(MODERN / 'fits/full/free/R/fit.FITRES.TEXT', 'SN')}[cid]
        for field in fields:
            record[f'old99_minus_modern99_{field}'] = index[(99, 'free', 'R', 1.0, cid)][field] - float(modern[field])
        baseline.append(record)
    write_csv(OUT / 'author-and-version-baseline.csv', baseline)
    pairs = []
    for timing in ('free', 'fixed'):
        for start in (1.0, .85, 1.15):
            for cid in protocol['cohort']:
                for left, right in (('B', 'A'), ('H', 'R'), ('A', 'R')):
                    l = index[(94, timing, left, start, cid)]; r = index[(94, timing, right, start, cid)]
                    pairs.append({'CID': cid, 'timing': timing, 'start': start, 'contrast': left + '-' + right,
                                  **{field: l[field] - r[field] for field in ('DLMAG', 'AV', 'STRETCH', 'PKMJD')},
                                  'label': 'Algorithm-conditional difference, not validated physical correction'})
    write_csv(OUT / 'paired-start-responses.csv', pairs)
    primary = [c for c in cases if c['law'] == 94]
    summary = {'protocol_sha256': sha(protocol_path), 'fits': len(fits), 'primary_cases': len(primary),
               'all_ERRFLAG_zero': all(f['ERRFLAG_FIT'] == 0 for f in fits),
               'all_covariance_positive': all(f['covariance_min_eigenvalue'] > 0 for f in fits),
               'primary_same_masks': sum(c['same_accepted_mask'] for c in primary),
               'primary_distance_spreads_gt_0p001': sum(c['DLMAG_range'] > .001 for c in primary),
               'primary_distance_spreads_gt_0p01': sum(c['DLMAG_range'] > .01 for c in primary),
               'primary_max_distance_spread': max(c['DLMAG_range'] for c in primary),
               'most_unstable_primary': sorted(primary, key=lambda c: c['DLMAG_range'], reverse=True)[:5],
               'nominal_old94_vs_author_distance_abs_max': max(abs(r['old94_minus_author_DLMAG']) for r in baseline),
               'nominal_old94_vs_author_distance_abs_median': float(np.median([abs(r['old94_minus_author_DLMAG']) for r in baseline])),
               'no_preferred_start_or_physical_correction': True,
               'objective_limit': 'FITCHI2 excludes priors. Covariance and peak-prior centers vary across iterative native fits; same mask alone is not a common likelihood.',
               'source_sha256': sha(Path(__file__)), 'verified_output_hashes': hashes}
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (OUT / 'executed-summary-source.py').write_bytes(Path(__file__).read_bytes())
    print(json.dumps({k: v for k, v in summary.items() if k != 'verified_output_hashes'}, indent=2))


if __name__ == '__main__':
    main()
