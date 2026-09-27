"""A frozen proposal learned from a complete, audited anchored replacement."""
import argparse
from contextlib import ExitStack
from datetime import datetime, timezone
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
from scipy.linalg import solve_triangular
from scipy.special import logsumexp

from independence_proposal import FrozenMixture, ROOT, digest, relative

HERE = Path(__file__).resolve().parent
DESIGN = HERE/'anchored-proposal-design.json'
VALIDATION = ROOT/'studies/unified_cosmology/results/inference/anchored-proposal-validation.json'
KIND = 'anchored_bridge_weighted_v1'
SCHEMA = 'anchored-bridge-weighted-proposal-v1'
SCOPE = 'Numerical proposal only; no scientific posterior qualification follows from training.'


def training_settings(report):
    parent = report['parent_settings']
    return {**{k: parent[k] for k in ['model','evolution','calibration']},
            'sample': 'pantheon_shoes_anchored', 'native_accuracy': 1}


def sources():
    return {relative(p): digest(p) for p in [Path(__file__), DESIGN, HERE/'independence_proposal.py',
                                            HERE/'independence-sampling-design.json']}


def verify_hashes(mapping):
    for path, checksum in mapping.items():
        assert digest(ROOT/path) == checksum, 'Frozen proposal evidence changed: '+path


def partition(groups, chronology):
    groups, chronology = np.asarray(groups), np.asarray(chronology)
    assert groups.ndim == chronology.ndim == 1 and len(groups) == len(chronology) >= 2000
    assert set(groups) == {0, 1, 2, 3}
    assert np.isfinite(chronology).all() and np.array_equal(chronology, chronology.astype(int))
    mask = np.zeros(len(groups), dtype=bool)
    sizes = []
    for group in range(4):
        index = np.flatnonzero(groups == group)
        assert len(index) >= 500 and len(set(chronology[index])) == len(index)
        ordered = index[np.argsort(chronology[index], kind='stable')]
        split = int(.8*len(ordered))
        assert len(ordered)-split >= 100
        mask[ordered[:split]] = True
        sizes.append({'chain': group, 'selected': len(index), 'training': split, 'withheld': len(index)-split})
    return mask, sizes


def fit(points, logweights, groups, chronology, reference):
    points, logweights = np.asarray(points), np.asarray(logweights)
    assert points.shape == (len(logweights), reference.d) and logweights.ndim == 1
    assert len(points) == len(groups) == len(chronology)
    assert np.isfinite(points).all() and np.isfinite(logweights).all()
    mask, sizes = partition(groups, chronology)
    weights = np.exp(logweights[mask]-logsumexp(logweights[mask]))
    ess = float(1/np.sum(weights**2))
    assert ess >= max(20, 2*reference.d), 'Insufficient effective training support for this proposal construction.'
    white = solve_triangular(reference.L, (points[mask]-reference.mean).T, lower=True).T
    mean_white = weights@white
    centered = white-mean_white
    covariance_white = (centered*weights[:, None]).T@centered
    regularized = .95*covariance_white+.05*np.eye(reference.d)
    mean = reference.mean+reference.L@mean_white
    covariance = reference.L@regularized@reference.L.T
    covariance = (covariance+covariance.T)/2
    proposal = FrozenMixture(mean, covariance, reference.names)
    statistics = {'raw_training_weight_ESS': ess, 'largest_training_weight': float(weights.max()),
                  'training_points': int(mask.sum()), 'withheld_points': int((~mask).sum()), 'partition': sizes,
                  'training_chain_weight_fractions': {str(g): float(weights[np.asarray(groups)[mask] == g].sum()) for g in range(4)},
                  'empirical_white_covariance_eigenvalues': np.linalg.eigvalsh(covariance_white).tolist(),
                  'regularized_white_covariance_eigenvalues': np.linalg.eigvalsh(regularized).tolist(),
                  'regularization': .05, 'posterior_qualification': False}
    return proposal, mask, statistics


def read_reference(folder):
    folder = Path(folder).resolve()
    record_path = folder/'proposal.json'
    record = json.loads(record_path.read_text())
    assert record['source_sha256'] == digest(HERE/'independence_proposal.py')
    assert record['design_sha256'] == digest(HERE/'independence-sampling-design.json')
    for name, item in record['snapshot_files'].items():
        assert digest(folder/name) == item['sha256']
    return FrozenMixture.read(folder/'proposal.npz', record['proposal_sha256']), record


def validation_guard():
    value = json.loads(VALIDATION.read_text())
    assert value['status'] == 'passed_synthetic_anchored_proposal_validation'
    assert value['proposal_source_sha256'] == sources()
    verify_hashes(value['source_sha256'])
    return {relative(VALIDATION): digest(VALIDATION)}


def freeze(bridge_path, reference_folder, output):
    import anchored_bridge as bridge
    from exact_correction import verify_record
    bridge_path, reference_folder, output = [Path(p).resolve() for p in [bridge_path, reference_folder, output]]
    assert output.is_relative_to(ROOT/'.work') and not output.exists()
    validation = validation_guard()
    report = json.loads(bridge_path.read_text())
    assert report['status'] in {'qualified_conditional_anchored_bridge', 'insufficient_anchored_overlap_or_stability'}
    assert report['parent_qualified_under_declared_numerical_gates'] is True
    manifest = ROOT/report['cache_manifest_path']
    assert digest(manifest) == report['cache_manifest_sha256']
    cached = json.loads(manifest.read_text()); verify_hashes(cached)
    correction = ROOT/report['correction_summary_path']
    assert digest(correction) == report['correction_summary_sha256']
    correction_value = json.loads(correction.read_text())
    selection_path = ROOT/correction_value['selection_path']
    selection = json.loads(selection_path.read_text())
    assert len(selection['points']) >= 2000
    for index in range(len(selection['points'])):
        assert relative(manifest.parent/f'{index:05d}.json') in cached
    def forbidden(*args, **kwargs):
        raise RuntimeError('Proposal training may replay only a complete bridge, never evaluate missing data.')
    import camb
    import cobaya.model
    with ExitStack() as stack:
        for owner, name in [(bridge, 'background_densities'), (camb, 'get_background'),
                            (camb, 'get_results'), (camb, 'get_transfer_functions'), (cobaya.model, 'get_model')]:
            stack.enter_context(patch.object(owner, name, forbidden))
        replay = bridge.actual(selection_path.parent.parent, correction, manifest.parent)
    assert replay == report, 'Saved bridge does not reproduce from its current qualified parent.'
    reference, reference_record = read_reference(reference_folder)
    names = [k for k,v in report['native_target']['configuration']['params'].items() if isinstance(v,dict) and 'prior' in v]
    assert names == reference.names
    points, weights = [], []
    for index, point in enumerate(selection['points']):
        native_path = selection_path.parent/f'{index:05d}.json'
        native = json.loads(native_path.read_text()); verify_record(native)
        row = bridge.read_record(manifest.parent/f'{index:05d}.json', report['bridge_identity'], index, point, digest(native_path))
        assert row['status'] == 'finite_anchored_replacement'
        assert abs(row['source_SN_loglike_closure']) <= json.loads(bridge.DESIGN.read_text())['source_closure_absolute_tolerance']
        points.append([point[name] for name in names]); weights.append(row['target_logweight'])
    chronology = np.array([row['expanded_index'] for row in selection['locations']])
    groups = np.array(selection['groups']); points = np.array(points); weights = np.array(weights)
    proposal, mask, statistics = fit(points, weights, groups, chronology, reference)
    lineage_path = ROOT/report['lineage_path']
    lineage = json.loads(lineage_path.read_text())
    bindings = {}
    for mapping in [cached, lineage['qualified_parent_inputs'], report['source_sha256'], validation]:
        for path, checksum in mapping.items():
            if path in bindings:
                assert bindings[path] == checksum, 'Conflicting proposal lineage.'
            bindings[path] = checksum
    for p in [bridge_path, manifest, correction, selection_path, lineage_path, reference_folder/'proposal.json', reference_folder/'proposal.npz']:
        bindings[relative(p)] = digest(p)
    verify_hashes(bindings)
    output.mkdir(parents=True)
    (output/'reference-proposal.npz').write_bytes((reference_folder/'proposal.npz').read_bytes())
    np.savez(output/'training.npz', points=points, logweights=weights, groups=groups,
             chronology=chronology, training_mask=mask, names=names)
    np.savez(output/'proposal.npz', mean=proposal.mean, cov=proposal.covariance, names=names)
    record = {'schema': SCHEMA, 'proposal_kind': KIND, 'created_utc': datetime.now(timezone.utc).isoformat(),
              'source_sha256': sources(), 'input_sha256': bindings,
              'proposal_sha256': digest(output/'proposal.npz'),
              'snapshot_sha256': {name: digest(output/name) for name in ['training.npz','reference-proposal.npz']},
              'names': names, 'statistics': statistics,
              'training_target_identity': report['native_target']['identity'], 'training_settings': training_settings(report),
              'training_parent_settings': report['parent_settings'],
              'reference_training_target_identity': reference_record['parent_target_identity'],
              'reference_proposal_record_path': relative(reference_folder/'proposal.json'),
              'bridge_status': report['status'], 'bridge_report_path': relative(bridge_path),
              'scope': SCOPE,
              'design': json.loads(DESIGN.read_text())}
    (output/'proposal.json').write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
    verify_hashes(bindings)
    return record


def load(folder, info):
    folder = Path(folder).resolve()
    record_path = folder/'proposal.json'; record = json.loads(record_path.read_text())
    assert record['schema'] == SCHEMA and record['proposal_kind'] == KIND
    assert record['source_sha256'] == sources()
    assert record['design'] == json.loads(DESIGN.read_text())
    verify_hashes(record['input_sha256'])
    assert record['scope'] == SCOPE
    for key in ['bridge_report_path','reference_proposal_record_path']:
        assert record[key] in record['input_sha256']
    bridge_report = json.loads((ROOT/record['bridge_report_path']).read_text())
    assert record['training_target_identity'] == bridge_report['native_target']['identity']
    assert record['training_settings'] == training_settings(bridge_report)
    assert record['training_parent_settings'] == bridge_report['parent_settings']
    assert record['bridge_status'] == bridge_report['status']
    assert record['bridge_status'] in {'qualified_conditional_anchored_bridge', 'insufficient_anchored_overlap_or_stability'}
    assert bridge_report['parent_qualified_under_declared_numerical_gates'] is True
    reference_report = json.loads((ROOT/record['reference_proposal_record_path']).read_text())
    assert record['reference_training_target_identity'] == reference_report['parent_target_identity']
    assert record['snapshot_sha256']['reference-proposal.npz'] == reference_report['proposal_sha256']
    for name, checksum in record['snapshot_sha256'].items():
        assert digest(folder/name) == checksum
    assert set(record['snapshot_sha256']) == {'training.npz','reference-proposal.npz'}
    proposal = FrozenMixture.read(folder/'proposal.npz', record['proposal_sha256'])
    reference = FrozenMixture.read(folder/'reference-proposal.npz', record['snapshot_sha256']['reference-proposal.npz'])
    with np.load(folder/'training.npz', allow_pickle=False) as snapshot:
        assert set(snapshot) == {'points','logweights','groups','chronology','training_mask','names'}
        assert snapshot['names'].tolist() == reference.names == proposal.names == record['names']
        rebuilt, mask, statistics = fit(snapshot['points'], snapshot['logweights'], snapshot['groups'], snapshot['chronology'], reference)
        assert np.array_equal(snapshot['training_mask'], mask)
    assert statistics == record['statistics']
    assert np.array_equal(rebuilt.mean, proposal.mean) and np.array_equal(rebuilt.covariance, proposal.covariance)
    expected = [k for k,v in info['params'].items() if isinstance(v,dict) and 'prior' in v]
    assert expected == proposal.names
    assert not any(info['params'][n].get('periodic', False) for n in expected)
    return proposal, {'proposal_kind': KIND, 'proposal_file': str(folder/'proposal.npz'),
                      'proposal_sha256': record['proposal_sha256'], 'proposal_record_path': relative(record_path),
                      'proposal_record_sha256': digest(record_path),
                      'training_target_identity': record['training_target_identity'],
                      'training_settings': record['training_settings'], 'role': record['scope']}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bridge', type=Path, required=True)
    p.add_argument('--reference-proposal', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    result = freeze(args.bridge, args.reference_proposal, args.output)
    print(json.dumps({'scope': result['scope'], 'statistics': result['statistics']}))


if __name__ == '__main__':
    main()
