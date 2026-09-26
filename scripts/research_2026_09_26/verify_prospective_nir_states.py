"""Independent engineering NIR checks against native photons and saved states."""
from pathlib import Path
from collections import Counter
import os
for key in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[key] = '1'
import argparse
import hashlib
import json
import numpy as np
from astropy.io import fits
from review_csp_native_states import parse


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('base', type=Path)
    ap.add_argument('out', type=Path)
    args = ap.parse_args()
    p = args.base
    args.out.mkdir(exist_ok=False)
    source = p/'readme-adapted/ledger/PTE'
    hpath = source/'PTE_HEAD.FITS'
    fpath = source/'PTE_PHOT.FITS'
    with fits.open(hpath) as f:
        head = f[1].data.copy()
    with fits.open(fpath) as f:
        phot = f[1].data.copy()
    wanted = {}
    for row in head:
        rows = phot[int(row['PTROBS_MIN'])-1:int(row['PTROBS_MAX'])]
        wanted[str(row['SNID']).strip()] = Counter(
            (str(r['BAND']).strip(), float(r['MJD']), float(r['FLUXCAL']),
             float(r['FLUXCALERR'])) for r in rows
            if str(r['BAND']).strip() in ['J', 'H'])
    assert len(wanted) == 8 and all(sum(v.values()) == 6 for v in wanted.values())
    null_head = p/'derived-start-only-bytes/null/PTE/PTE_HEAD.FITS'
    est_head = p/'derived-start-only-bytes/estimated/PTE/PTE_HEAD.FITS'
    assert null_head.read_bytes() == hpath.read_bytes()
    with fits.open(est_head) as f:
        est = {str(r['SNID']).strip(): float(r['PEAKMJD']) for r in f[1].data}
    truth = {str(r['SNID']).strip(): float(r['PEAKMJD']) for r in head}
    inputs = [Path(__file__), Path(__file__).with_name('review_csp_native_states.py'),
              hpath, fpath, null_head, est_head]
    results = {}
    total = 0
    for arm in ['true', 'estimated', 'null']:
        for suffix, nit in ([('12', 12)] if arm == 'null' else
                            [('12', 12), ('9', 9), ('minus', 12), ('plus', 12)]):
            label = 'NIR_'+arm+suffix
            log = p/'fits-start-only'/label/'native.log'
            inputs.append(log)
            records, entries = parse(log)
            assert len(records) == nit*8
            total += len(records)
            expected_peak = est if arm == 'estimated' else truth
            finals, checks = {}, []
            for cid in wanted:
                seq = [r for r in records if r['CID'] == cid]
                assert [r['iteration'] for r in seq] == list(range(1, nit+1))
                covs = []
                for r in seq:
                    actual = Counter((x['band'], x['values'][0], x['values'][4],
                                      x['values'][5]) for x in r['rows'])
                    assert actual == wanted[cid] and r['n'] == 6
                    x = np.array([a['values'] for a in r['rows']])
                    w = np.array(r['weights'])
                    w = w if r['cov'] else np.diag(w[:, 0])
                    assert np.isfinite(x).all() and np.isfinite(w).all()
                    assert np.max(abs(w-w.T)) <= max(1e-12, 1e-10*np.max(abs(w)))
                    chol = np.linalg.cholesky(w)
                    q = float(np.sum((chol.T@(x[:, 4]-x[:, 2]))**2))
                    o = r['objective']
                    gap = abs(q-(o[0]-o[1]-o[2]))
                    assert gap <= 1e-7 and o[4] == 1 and o[5] == 0
                    assert o[6] == expected_peak[cid]
                    assert -20 <= min(x[:, 1]) <= max(x[:, 1]) <= 70
                    alpha = float(x[:, 2]@w@x[:, 4]/(x[:, 2]@w@x[:, 2]))
                    assert alpha > 0
                    covs.append(np.linalg.inv(w))
                    checks.append(dict(CID=cid, iteration=r['iteration'], Q_gap=gap,
                                       fixed_C_D_gap=float(-2.5*np.log10(alpha))))
                a, b = seq[-2:]
                l = np.linalg.cholesky(covs[-2])
                dc = np.linalg.solve(l, covs[-1]-covs[-2])
                dc = np.linalg.solve(l, dc.T).T
                metric = float(np.linalg.norm(dc, 2))
                assert metric <= .001
                assert abs(b['objective'][3]-a['objective'][3]) <= .001
                assert abs(checks[-1]['fixed_C_D_gap']) <= .001
                finals[cid] = dict(D=b['objective'][3], peak=b['objective'][6],
                                   last_C_whitened_norm=metric)
            support = [s.split() for s in log.read_text().splitlines()
                       if s.startswith('PROSP_SUPPORT')]
            assert len(support) == 8 and {s[1] for s in support} == set(wanted)
            assert all(len(s) == 11 and int(s[2]) > 0 and
                       int(s[3]) == int(s[4]) == 0 for s in support)
            results[label] = dict(final=finals, checks=checks,
                                 first_D={e[0]: float(e[5]) for e in entries
                                          if int(e[1]) == 1})
    for arm in ['true', 'estimated']:
        base = results['NIR_'+arm+'12']
        for suffix in ['9', 'minus', 'plus']:
            case = results['NIR_'+arm+suffix]
            gap = max(abs(case['final'][c]['D']-base['final'][c]['D']) for c in wanted)
            assert gap <= .001
            case['max_final_D_difference_from_12'] = gap
            if suffix in ['minus', 'plus']:
                shift = -.2 if suffix == 'minus' else .2
                assert all(abs(case['first_D'][c]-base['first_D'][c]-shift) < 1e-4
                           for c in wanted)
    def tokens(name):
        return [s.split() for s in (p/'fits-start-only'/name/'native.log').read_text().splitlines()
                if s.startswith('CSP_') or s.startswith('PROSP_SUPPORT')]
    assert tokens('NIR_true12') == tokens('NIR_null12')
    result = dict(gate_pass=True, callbacks=total, exact_null_file_and_CSP=True,
                  pure_peak_and_fixed_photons=True, cases=results,
                  scope='Engineering validation only; no population bias or cosmology.')
    (args.out/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    (args.out/'manifest.json').write_text(json.dumps({str(x.resolve()):
        hashlib.sha256(x.read_bytes()).hexdigest() for x in inputs}, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'cases'}))


if __name__ == '__main__':
    main()
