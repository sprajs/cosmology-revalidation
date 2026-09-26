#!/usr/bin/env python3
"""Source-level FITRES/LCPLOT coherence check, independent of pilot refits."""
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'runs/research_2026_09_26/raisin_simulation_timing_pilot'
BASE = ROOT / 'runs/research_2026_09_26/raisin_simulation_execution_review/nominal_fitres'
FIT = BASE / 'nir.FITRES.gz'
PLOT = BASE / 'lcplot-feasibility/FITOPT000.LCPLOT.gz'

def rows(path, marker):
    with gzip.open(path, 'rt') as f:
        names = None
        for line in f:
            if line.startswith('VARNAMES:'):
                names = line.split()[1:]
            elif line.startswith(marker):
                yield dict(zip(names, line.split()[1:]))

fit = {}
for i, r in enumerate(rows(FIT, 'SN:')):
    if i >= 500:
        break
    fit[r['CID']] = r
plot = collections.defaultdict(list)
for r in rows(PLOT, 'OBS:'):
    if r['DATAFLAG'] == '1':
        plot[r['CID']].append(r)
assert len(fit) == len(plot) == 500
assert set(fit) == set(plot)
detail = []
for cid in fit:
    a = fit[cid]
    b = plot[cid]
    n = len(b)
    ndof = int(a['NDOF'])
    q = float(a['FITCHI2'])
    plot_q = sum(float(r['CHI2']) for r in b)
    detail.append({'CID': int(cid), 'FITRES_NDOF': ndof,
                   'LCPLOT_DATAFLAG1_rows': n, 'epoch_count_implied_if_one_free_DLMAG': ndof+1,
                   'count_match': n == ndof+1, 'FITRES_FITCHI2': q,
                   'LCPLOT_sum_printed_CHI2': plot_q, 'abs_Q_difference': abs(plot_q-q)})
result = {
    'sources': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in (FIT, PLOT)},
    'n_CIDs': 500,
    'count_matches': sum(r['count_match'] for r in detail),
    'count_mismatches': sum(not r['count_match'] for r in detail),
    'mismatch_CIDs': [r['CID'] for r in detail if not r['count_match']],
    'Q_printed_sum_within_0p1': sum(r['abs_Q_difference'] <= .1 for r in detail),
    'Q_printed_sum_over_1': sum(r['abs_Q_difference'] > 1 for r in detail),
    'first_20': detail[:20],
    'qualification': 'LCPLOT per-epoch CHI2 is rounded and may not sum to FITRES objective under covariance; count mismatch is independent of this qualification.'
}
with (OUT / 'archive-coherence-500.csv').open('w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(detail[0]))
    w.writeheader()
    w.writerows(detail)
(OUT / 'archive-coherence.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('first_20','sources','mismatch_CIDs')}, indent=2))
