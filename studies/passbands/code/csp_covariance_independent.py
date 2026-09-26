#!/usr/bin/env python3
"""Independent read-only algebra audit of CSP native iteration covariance."""
from __future__ import annotations

import csv
import hashlib
import json
from math import log
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'runs/research_2026_09_26/csp_covariance_independent'
LOG = ROOT / 'runs/research_2026_09_26/csp_native_filter_response/cohort-execution/full/fits/nominal/fit.log'
K = 0.4 * log(10.0)


def parse(log: Path):
    states = {}
    rows = []
    active = None
    for line in log.read_text().splitlines():
        parts = line.split()
        if not parts:
            continue
        tag = parts[0]
        if tag == 'CSP_ROW:':
            rows.append({'CID': parts[1], 'iteration': int(parts[2]), 'index': int(parts[3]),
                         'epoch': int(parts[4]), 'band': parts[5], 'v': np.array([float(v) for v in parts[6:]])})
        elif tag == 'CSP_OBJECTIVE:':
            cid, it, n, cov = parts[1], int(parts[2]), int(parts[3]), parts[4] == 'T'
            assert len(rows) == n and all(r['CID'] == cid and r['iteration'] == it for r in rows)
            assert [r['index'] for r in rows] == list(range(1, n+1))
            active = {'CID': cid, 'iteration': it, 'cov': cov, 'n': n, 'rows': rows,
                      'objective': np.array([float(v) for v in parts[5:]]), 'weight_rows': []}
            assert len(active['objective']) == 10
            rows = []
        elif tag in ('CSP_WROW:', 'CSP_WDIAG:'):
            assert active is not None
            assert parts[1] == active['CID'] and int(parts[2]) == active['iteration']
            assert int(parts[3]) == len(active['weight_rows'])+1
            v = [float(x) for x in parts[4:]]
            assert len(v) == (active['n'] if active['cov'] else 1)
            active['weight_rows'].append(v)
            if len(active['weight_rows']) == active['n']:
                W = np.array(active['weight_rows'])
                active['W'] = W if active['cov'] else np.diag(W[:,0])
                key = (active['CID'], active['iteration'])
                assert key not in states
                states[key] = active
                active = None
    assert not rows and active is None
    return states


def main():
    expected = json.loads((OUT/'protocol.json').read_text())['inputs']
    for name, sha in expected.items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == sha
    states = parse(LOG)
    cids = sorted({cid for cid, _ in states})
    assert len(cids) == 42 and len(states) == 126
    records=[]
    for cid in cids:
        s1,s2,s3 = (states[(cid,j)] for j in (1,2,3))
        ids2 = [(r['epoch'],r['band'],r['v'][0]) for r in s2['rows']]
        ids3 = [(r['epoch'],r['band'],r['v'][0]) for r in s3['rows']]
        same_rows = ids2 == ids3
        a2=np.stack([r['v'] for r in s2['rows']]);a3=np.stack([r['v'] for r in s3['rows']])
        same_data=bool(np.array_equal(a2[:,4:6],a3[:,4:6]))
        result={'CID':cid,'n2':s2['n'],'n3':s3['n'],'same_rows':same_rows,'same_data_error':same_data,
                'D1':float(s1['objective'][3]),'D2':float(s2['objective'][3]),'D3':float(s3['objective'][3])}
        if not(same_rows and same_data):
            result['status']='row/data mismatch; no scaling or score calculated'
            records.append(result)
            continue
        W2,W3=s2['W'],s3['W']
        C2=np.linalg.inv(W2);C3=np.linalg.inv(W3)
        dataerr=a3[:,5]
        E=np.diag(dataerr**2)
        M2=C2-E
        scale=10**(-0.4*(result['D2']-result['D1']))
        C3_pred=E+scale**2*M2
        delta=C3_pred-C3
        invsqrt=1/dataerr
        model_eigs=np.linalg.eigvalsh(M2*invsqrt[:,None]*invsqrt[None,:])
        model_eigs3=np.linalg.eigvalsh((C3-E)*invsqrt[:,None]*invsqrt[None,:])
        # At state 3, model f and residual are evaluated at D3; covariance
        # was rebuilt before this fit at D2. Hold those coordinates distinct.
        f=a3[:,2];resid=a3[:,4]-f
        f_D=-K*f
        C_D=-2*K*(C3-E)
        Wr=W3@resid
        mean_score=-2*float(f_D@Wr)
        quadratic_score=-float(Wr@C_D@Wr)
        logdet_score=float(np.trace(W3@C_D))
        result.update(status='calculated',scale=scale,model_min_eigen_scaled=float(model_eigs[0]),
                      model_max_eigen_scaled=float(model_eigs[-1]),
                      final_model_min_eigen_scaled=float(model_eigs3[0]),
                      final_model_max_eigen_scaled=float(model_eigs3[-1]),
                      covariance_prediction_max_abs=float(np.max(abs(delta))),
                      covariance_prediction_relative_Frobenius=float(np.linalg.norm(delta)/np.linalg.norm(C3)),
                      covariance_prediction_whitened_max_abs=float(np.max(abs(delta*invsqrt[:,None]*invsqrt[None,:]))),
                      W2_C2_identity_max_abs=float(np.max(abs(W2@C2-np.eye(len(W2))))),
                      W3_C3_identity_max_abs=float(np.max(abs(W3@C3-np.eye(len(W3))))),
                      frozen_C_mean_term=mean_score, gaussian_cov_quadratic_term=quadratic_score,
                      gaussian_cov_logdet_term=logdet_score,
                      gaussian_total_score=mean_score+quadratic_score+logdet_score,
                      frozen_C_model_mean_Q=float(resid@W3@resid),
                      actual_model_mean_ratio_max_abs=float(np.max(abs(a3[:,2]/a2[:,2]-10**(-0.4*(result['D3']-result['D2']))))) if np.all(a2[:,2]!=0) else None)
        records.append(result)
    p=OUT/'per-object.csv'
    with p.open('w',newline='') as f:
        cols=list(dict.fromkeys(k for row in records for k in row))
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(records)
    numeric=[r for r in records if r['status']=='calculated']
    summary={'input_log':str(LOG.relative_to(ROOT)),'input_log_sha256':hashlib.sha256(LOG.read_bytes()).hexdigest(),
             'n_cids':len(cids),'n_calculated':len(numeric),'n_row_data_failures':len(records)-len(numeric),
             'max_covariance_prediction_relative_Frobenius':max(r['covariance_prediction_relative_Frobenius'] for r in numeric),
             'max_covariance_prediction_whitened_max_abs':max(r['covariance_prediction_whitened_max_abs'] for r in numeric),
             'min_model_eigen_scaled':min(r['model_min_eigen_scaled'] for r in numeric),
             'min_final_model_eigen_scaled':min(r['final_model_min_eigen_scaled'] for r in numeric),
             'score_extrema':{k:[min(r[k] for r in numeric),max(r[k] for r in numeric)] for k in ('frozen_C_mean_term','gaussian_cov_quadratic_term','gaussian_cov_logdet_term','gaussian_total_score')},
             'records':records}
    (OUT/'result.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    (OUT/'manifest.json').write_text(json.dumps({str(x.relative_to(ROOT)):hashlib.sha256(x.read_bytes()).hexdigest() for x in (OUT/'protocol.json',p,OUT/'result.json',Path(__file__))},indent=2)+'\n')


if __name__ == '__main__': main()
