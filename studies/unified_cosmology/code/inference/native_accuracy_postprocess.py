"""Read-only summaries for an explicitly qualified accuracy-two target.

Pure aggregation kernels are shared with the original analyses. Their old
filesystem/configuration consumers are not called or monkeypatched. This
module performs no physical calculation and never changes a parent record.
"""
import argparse
import copy
import json
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

import probe_omission as omission
import quantile_precision as quantiles
from measurement_summary import PARAMETERS, weighted_fraction
from target_identity import canonical

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DESIGN = HERE/'native-accuracy-postprocess-design.json'


def arrays(points, derived, logweights, groups):
    """Pure numerical view; this function alone cannot qualify observations."""
    assert len(points) == len(derived) == len(logweights) == len(groups)
    names = sorted(set(points[0]) | set(derived[0]))
    assert all(set(p) == set(points[0]) and set(d) == set(derived[0])
               and not set(p).intersection(d) for p, d in zip(points, derived))
    values = {k: np.array([dict(p, **d)[k] for p, d in zip(points, derived)]) for k in names}
    lw = np.asarray(logweights, dtype=float)
    assert np.isfinite(lw).all() and all(np.isfinite(x).all() for x in values.values())
    weights = np.exp(lw-logsumexp(lw))
    return values, weights


def summaries(values, weights, groups, posterior):
    """Aggregate native2 arrays; priors fixed by the parent remain explicit."""
    names = [k for k in PARAMETERS if k in values]
    for k in names:
        assert abs(weights@values[k]-posterior[k]['mean']) < 1e-8
        assert np.max(abs(quantiles.weighted_quantiles(values[k], weights)
                          - posterior[k]['quantiles_025_16_50_84_975'])) <= 1e-12
    matrix = np.column_stack([values[k] for k in names])
    centered = matrix-weights@matrix
    signs = {k+'_below_zero': weighted_fraction(values[k] < 0, weights, groups)
             for k in ('q0', 'q05', 'q1', 'j0')}
    signs['accelerating_with_decreasing_scale_factor_acceleration_today'] = weighted_fraction(
        (values['q0'] < 0) & (values['j0'] < 0), weights, groups)
    correlations = {}
    for left, right in [('w','wa'), ('w','H0'), ('wa','H0'),
                        ('epsilon','w'), ('epsilon','wa'), ('q0','j0')]:
        if left not in values or right not in values:
            continue
        x, y = values[left]-weights@values[left], values[right]-weights@values[right]
        denominator = np.sqrt((weights@x**2)*(weights@y**2))
        correlations[left+','+right] = float(weights@(x*y)/denominator) if denominator else None
    return {'status':'qualified_conditional_native_accuracy2_summary',
            'posterior':{k:posterior[k] for k in names}, 'conditional_sign_fractions':signs,
            'weighted_correlations':correlations,
            'weighted_covariance':{'parameter_order':names,
                'matrix':((centered*weights[:,None]).T@centered).tolist(),
                'normalization':'Normalized raw posterior weights; no Bessel correction.'}}


def remove_components(point, derived, native_loglikes, native_logpost,
                      proposal_logpost, logweight, removed, expected, tolerance):
    """Exact native2/q density algebra, without fabricated proposal components."""
    assert set(native_loglikes) == set(expected)
    assert len(removed) == len(set(removed)) and set(removed) <= set(native_loglikes)
    assert np.isfinite([*native_loglikes.values(), native_logpost, proposal_logpost, logweight]).all()
    assert abs(logweight-native_logpost+proposal_logpost) <= tolerance
    removed_sum = float(sum(native_loglikes[k] for k in removed))
    target_logpost = float(native_logpost-removed_sum)
    target_weight = float(logweight-removed_sum)
    closure = float(target_logpost-proposal_logpost-target_weight)
    assert abs(closure) <= tolerance
    # This is only the documented arithmetic view expected by the old pure
    # summarizer. It is never written as an original native record or passed to
    # the old filesystem qualifier. All values explicitly originate at native2.
    view = {'point':point, 'derived':derived, 'log_weight':float(logweight)}
    row = {'status':'finite_component_removal', 'target_logweight':target_weight,
           'removed_native_loglike_sum':removed_sum,
           'target_weight_accounting_closure':closure,
           'target_logpost_in_inherited_normalization':target_logpost}
    return view, row


def chronology(points, groups, locations, proposal_logposts):
    groups = np.asarray(groups)
    assert len(groups) == len(points) == len(locations) == len(proposal_logposts) == 2000
    assert np.array_equal(groups, np.repeat(np.arange(4), 500))
    report = {}
    for group in range(4):
        ids = np.flatnonzero(groups == group)
        loc = [locations[i] for i in ids]
        assert len({r['chain'] for r in loc}) == 1
        assert all(isinstance(r[k], int) and not isinstance(r[k], bool) and r[k] >= 0
                   for r in loc for k in ('expanded_index','row'))
        expanded = np.array([r['expanded_index'] for r in loc])
        rows = np.array([r['row'] for r in loc])
        assert np.all(np.diff(expanded) >= 0) and np.all(np.diff(rows) >= 0)
        repeats = np.flatnonzero(np.diff(expanded) == 0)
        for j in repeats:
            left, right = int(ids[j]), int(ids[j+1])
            assert rows[j] == rows[j+1] and points[left] == points[right]
            assert proposal_logposts[left] == proposal_logposts[right]
        report[str(group)] = {'slots':500, 'unique_expanded_indices':int(len(np.unique(expanded))),
                              'repeated_slots':int(len(repeats))}
    return report


def sources():
    paths = [Path(__file__), DESIGN, HERE/'native_accuracy_postprocess_validate.py',
             HERE/'native_accuracy_qualify.py', HERE/'probe_omission.py', omission.DESIGN,
             omission.GATES, HERE/'measurement_summary.py', HERE/'quantile_precision.py',
             HERE/'luminosity_sensitivity.py', HERE/'target_identity.py']
    return {omission.relative(p):omission.digest(p) for p in paths}


def actual(work, summary_path, action):
    # This mandatory fresh typed boundary precedes all record/config access.
    from native_accuracy_qualify import QualifiedNativeAccuracy2, qualify_run
    source = sources()
    validation_path = ROOT/'studies/unified_cosmology/results/inference/native-accuracy-postprocess-validation.json'
    validation_sha = omission.digest(validation_path)
    validation = json.loads(validation_path.read_text())
    assert validation['status'] == 'passed_native_accuracy2_postprocess_validation'
    assert validation['source_sha256'] == source and validation['physical_calls'] == 0
    target = qualify_run(work, summary_path)
    assert isinstance(target, QualifiedNativeAccuracy2) and target.numerical_accuracy == 2
    omission.verify_hashes(target.input_sha256)
    records, points = target.records, target.selected_points
    groups = np.asarray(target.groups)
    derived = [r['native2_derived'] for r in records]
    values, weights = arrays(points, derived, target.logweights, groups)
    assert np.max(abs(weights-target.normalized_weights)) < 1e-14
    chronology_result = chronology(points, groups, target.locations, target.proposal_logposts)
    config = copy.deepcopy(target.native_configuration)
    for name in ('AccuracyBoost','lAccuracyBoost','lSampleBoost'):
        assert config['theory']['camb']['extra_args'][name] == 2
    if action == 'summary':
        result = summaries(values, weights, groups, target.posterior_summary['posterior'])
    elif action == 'quantiles':
        selected = {k:values[k] for k in PARAMETERS if k in values}
        result = quantiles.block_quantile_precision(selected, target.logweights, groups)
        for k in selected:
            assert np.max(abs(np.array(result['baseline_quantiles'][k])
                - target.posterior_summary['posterior'][k]['quantiles_025_16_50_84_975'])) <= 1e-12
    else:
        assert action in ('omit-sn','omit-bao','omit-lensing')
        omitted = action.removeprefix('omit-')
        design = json.loads(omission.DESIGN.read_text())
        gates = json.loads(omission.GATES.read_text())['overlap_gates']
        removed = design['alternatives'][omitted]['removed_components']
        assert set(config['likelihood']) == set(design['required_parent_components'])
        views, rows = [], []
        for i, row in enumerate(records):
            view, removed_row = remove_components(points[i], derived[i], row['native2_loglikes'],
                row['native2_logpost'], target.proposal_logposts[i], target.logweights[i],
                removed, config['likelihood'], design['logweight_accounting_absolute_tolerance'])
            views.append(view); rows.append(removed_row)
        sampled = {k for k,v in config['params'].items() if isinstance(v,dict) and 'prior' in v}
        unused = [k for k, case in [('epsilon','sn'),('A_fg','lensing')] if k in sampled and omitted == case]
        result = omission.summarize_omission(views, rows, groups, gates, unused)
        for component in removed:
            config['likelihood'].pop(component)
        result.update(omitted_probe=omitted, removed_components=removed,
                      all_parent_priors_and_dimensions_retained=True,
                      unused_auxiliary_coordinates=unused)
    identity_payload = {'native_accuracy':2, 'parent_numerical_target_identity':target.numerical_target_identity,
                        'configuration':canonical(config), 'action':action, 'source_sha256':source}
    result.update(numerical_accuracy=2, numerical_target_identity=omission.identity(identity_payload),
        target_definition=identity_payload, parent_numerical_target_identity=target.numerical_target_identity,
        parent_proposal_target_identity=target.parent_proposal_target_identity,
        settings={k:target.parent_settings[k] for k in ('model','evolution','sample','calibration')},
        chronology=chronology_result, physical_calls=0, source_sha256=source,
        input_sha256=dict(target.input_sha256, **{omission.relative(validation_path):validation_sha}),
        interpretation=json.loads(DESIGN.read_text())['reporting'])
    omission.verify_hashes(target.input_sha256)
    assert omission.digest(validation_path) == validation_sha
    assert sources() == source
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--summary', type=Path, required=True)
    p.add_argument('--action', choices=('summary','quantiles','omit-sn','omit-bao','omit-lensing'), required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    assert not a.output.exists(), 'Use a fresh report; retain previous results.'
    result = actual(a.work, a.summary, a.action)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'], 'action':a.action, 'numerical_accuracy':2}))


if __name__ == '__main__':
    main()
