"""Independent completed-state checks; never certifies a missing start test."""
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
    data = p / 'generation-v5/ledger/output/PTE'
    head_path, = data.glob('*HEAD.FITS*')
    phot_path, = data.glob('*PHOT.FITS*')
    with fits.open(head_path) as f:
        head = f[1].data.copy()
    with fits.open(phot_path) as f:
        phot = f[1].data.copy()
    wanted = {}
    for h in head:
        rows = phot[int(h['PTROBS_MIN'])-1:int(h['PTROBS_MAX'])]
        wanted[str(h['SNID']).strip()] = Counter(
            (str(r['BAND']).strip(), float(r['MJD']), float(r['FLUXCAL']),
             float(r['FLUXCALERR'])) for r in rows)
        assert len(rows) == 117
    assert len(wanted) == 8
    result = []
    inputs = [head_path, phot_path, Path(__file__),
              Path(__file__).with_name('review_csp_native_states.py')]
    last_states = {}
    for name, nit in [('joint12', 12), ('joint9', 9), ('joint_minus', 12)]:
        log = p / 'fits-restricted' / name / 'native.log'
        inputs.append(log)
        records, entries = parse(log)
        assert len(records) == nit*8
        assert {r['CID'] for r in records} == set(wanted)
        checks = []
        final = {}
        for cid in wanted:
            seq = [r for r in records if r['CID'] == cid]
            assert [r['iteration'] for r in seq] == list(range(1, nit+1))
            covariances = []
            for r in seq:
                actual = Counter((x['band'], x['values'][0], x['values'][4],
                                  x['values'][5]) for x in r['rows'])
                assert actual == wanted[cid] and r['n'] == 117
                x = np.array([a['values'] for a in r['rows']])
                w = np.array(r['weights'])
                w = w if r['cov'] else np.diag(w[:, 0])
                assert np.isfinite(x).all() and np.isfinite(w).all()
                assert np.max(abs(w-w.T)) <= max(1e-12, 1e-10*np.max(abs(w)))
                l = np.linalg.cholesky(w)
                q = float(np.sum((l.T @ (x[:, 4]-x[:, 2]))**2))
                o = r['objective']
                gap = abs(q-(o[0]-o[1]-o[2]))
                assert gap <= 1e-7 and o[4] == 1 and o[5] == 0
                assert -20 <= min(x[:, 1]) <= max(x[:, 1]) <= 70
                a = float(x[:, 2] @ w @ x[:, 4] / (x[:, 2] @ w @ x[:, 2]))
                assert a > 0
                covariances.append(np.linalg.inv(w))
                checks.append(dict(CID=cid, iteration=r['iteration'], Q_gap=gap,
                                   fixed_C_D_gap=float(-2.5*np.log10(a))))
            a, b = seq[-2:]
            l = np.linalg.cholesky(covariances[-2])
            dc = np.linalg.solve(l, covariances[-1]-covariances[-2])
            dc = np.linalg.solve(l, dc.T).T
            metric = float(np.linalg.norm(dc, 2))
            assert metric <= .001
            assert abs(b['objective'][3]-a['objective'][3]) <= .001
            assert abs(b['objective'][6]-a['objective'][6]) <= .01
            assert abs(b['objective'][6]-57707.80078125) <= 4
            assert abs(checks[-1]['fixed_C_D_gap']) <= .001
            final[cid] = dict(D=b['objective'][3], peak=b['objective'][6],
                              last_C_whitened_norm=metric)
        support = [x.split() for x in log.read_text().splitlines()
                   if x.startswith('PROSP_SUPPORT')]
        assert len(support) == 8
        assert {x[1] for x in support} == set(wanted)
        assert all(x[0] == 'PROSP_SUPPORT' and len(x) == 11 and
                   int(x[2]) > 0 and int(x[3]) == int(x[4]) == 0 for x in support)
        first_entries = {e[0]: dict(D=float(e[5]), peak=float(e[8]))
                         for e in entries if int(e[1]) == 1}
        assert set(first_entries) == set(wanted)
        result.append(dict(case=name, callbacks=len(records), checks=checks,
                           final=final, first_entries=first_entries))
        last_states[name] = final
    base = result[0]
    for case in result[1:]:
        case['entry_peak_delta_from_nominal'] = {
            c: case['first_entries'][c]['peak']-base['first_entries'][c]['peak']
            for c in wanted}
        case['max_final_distance_delta'] = max(
            abs(case['final'][c]['D']-base['final'][c]['D']) for c in wanted)
        case['max_final_peak_delta'] = max(
            abs(case['final'][c]['peak']-base['final'][c]['peak']) for c in wanted)
        assert case['max_final_distance_delta'] <= .001
        assert case['max_final_peak_delta'] <= .01
    out = dict(completed_state_checks_pass=True, callbacks=264, cases=result,
               actual_displaced_start_gate_pass=False,
               scope='Three completed joint cases only. The minus NML initializer '
               'does not produce the prescribed optimizer-entry displacement. '
               'No paired timing or cosmological effect is certified.')
    (args.out/'result.json').write_text(json.dumps(out, indent=2)+'\n')
    (args.out/'manifest.json').write_text(json.dumps({str(x.resolve()):
        hashlib.sha256(x.read_bytes()).hexdigest() for x in inputs}, indent=2)+'\n')
    print(json.dumps({k: v for k, v in out.items() if k != 'cases'}, indent=2))
    for case in result[1:]:
        print(json.dumps({k: v for k, v in case.items()
                          if k not in ['checks', 'final', 'first_entries']}))


if __name__ == '__main__':
    main()
