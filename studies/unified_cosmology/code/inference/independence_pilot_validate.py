"""Independent saved-pilot density, prior-accounting and MH/hold audit."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import logsumexp
from scipy.stats import multivariate_normal, multivariate_t, norm
from independence_proposal import ROOT, digest, relative


def validate(proposal_folder, pilot_folder):
    proposal_folder, pilot_folder = Path(proposal_folder), Path(pilot_folder)
    result = json.loads((pilot_folder/'result.json').read_text())
    design = json.loads((pilot_folder/'design.json').read_text())
    for filename, expected in result['source_sha256'].items(): assert digest(ROOT/filename) == expected
    for filename, expected in result['artifact_sha256'].items(): assert digest(ROOT/filename) == expected
    manifest = json.loads((proposal_folder/'run-0.json').read_text())
    assert manifest['target_identity']['identity'] == result['target_identity']
    rows = [json.loads(line) for line in (pilot_folder/'records.jsonl').read_text().splitlines()]
    with np.load(proposal_folder/'proposal.npz') as data:
        names, mu, covariance = data['names'].tolist(), data['mean'], data['cov']
    with np.load(pilot_folder/'requests.npz') as data:
        assert data['names'].tolist() == names
        request, stored_q, logus = data['points'], data['logq'], data['log_uniform']
    assert len(rows) == len(request) == result['target_calls'] == 400
    # Standardize first: direct eigendecomposition of heterogeneous physical
    # units is an unnecessarily ill-conditioned independent density reference.
    scale = np.sqrt(np.diag(covariance)); correlation = covariance/np.outer(scale, scale)
    z = (request-mu)/scale
    independently = np.logaddexp(np.log(.9)+multivariate_normal.logpdf(z, np.zeros(len(names)), 1.05**2*correlation),
                                 np.log(.1)+multivariate_t.logpdf(z, np.zeros(len(names)), 2.4*correlation, df=5))-np.log(scale).sum()
    error = float(max(abs(independently-stored_q))); assert error < 1e-10
    current = None; occupied = []; closure = []; invalid = []; logweights = []
    for index, row in enumerate(rows):
        assert row['index'] == index and [row['point'][n] for n in names] == request[index].tolist()
        assert row['logq'] == stored_q[index] and row['log_uniform'] == logus[index]
        lp = 0.
        for name in names:
            value = row['point'][name]; p = manifest['target_identity']['configuration']['params'][name]['prior']
            if p.get('dist') == 'norm': lp += norm.logpdf(value, p['loc'], p['scale'])
            elif p['min'] <= value <= p['max']: lp -= np.log(p['max']-p['min'])
            else: lp = -np.inf; break
        assert row['prior_rejected'] == (not np.isfinite(lp))
        if row['finite']:
            assert all(v is not None for v in row['component_loglikes'].values())
            closure.append(abs(lp+sum(row['component_loglikes'].values())-row['target_logpost']))
            logweights.append(row['target_logpost']-independently[index])
        else:
            invalid.append({'index': index, 'prior_rejected': row['prior_rejected'],
                            'CPL_w_plus_wa': row['point'].get('w', -1)+row['point'].get('wa', 0)})
            logweights.append(-np.inf)
        if current is None: accept = row['finite']
        elif not row['finite']: accept = False
        else:
            ratio = min(0., row['target_logpost']-rows[current]['target_logpost']+independently[current]-independently[index])
            assert abs(ratio-row['log_acceptance']) < 1e-10
            accept = row['log_uniform'] < ratio
        assert accept == row['accepted']
        if accept: current = index
        assert row['occupied_candidate'] == current
        if current is not None: occupied.append(current)
    assert max(closure, default=0) < 1e-10
    weights = np.exp(np.asarray(logweights)-logsumexp(logweights))
    assert abs(1/(weights@weights)-result['iid_proposal_raw_importance_ESS']) < 1e-8
    unique, counts = np.unique(occupied, return_counts=True)
    assert len(unique) == result['unique_occupied_states'] and int(counts.max()) == result['longest_hold']
    return {'status': 'passed_independent_saved_pilot_audit', 'target_calls_checked': len(rows),
            'density_max_absolute_error': error, 'prior_plus_component_logdensity_error': max(closure, default=0),
            'all_MH_decisions_and_holding_counts_match': True, 'invalid_candidates': invalid,
            'target_identity': result['target_identity'],
            'source_sha256': digest(__file__), 'result_sha256': digest(pilot_folder/'result.json'),
            'design_sha256': digest(pilot_folder/'design.json'),
            'scope': 'Arithmetic and recorded-input identity only, no new target evaluations or posterior qualification.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('proposal_folder', type=Path); p.add_argument('pilot_folder', type=Path)
    a = p.parse_args(); r = validate(a.proposal_folder.resolve(), a.pilot_folder.resolve())
    (a.pilot_folder/'independent-audit.json').write_text(json.dumps(r, indent=2)+'\n')
    print(json.dumps({k:r[k] for k in ['status', 'target_calls_checked', 'density_max_absolute_error', 'prior_plus_component_logdensity_error']}, indent=2))
