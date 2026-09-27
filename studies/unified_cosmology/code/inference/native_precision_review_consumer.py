"""Verify a supplemental thermal precision report without computing cosmology.

An auxiliary closure repair is not a posterior qualification. Both original
and reviewed screen decisions remain visible in the return value.
"""
import argparse
from contextlib import ExitStack
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
from scipy.special import logsumexp

import native_posterior_precision as original
import native_precision_thermal_review as thermal
from exact_correction import verify_record
from measurement_summary import summarize_run
from target_identity import canonical

ROOT = original.ROOT
HERE = Path(__file__).resolve().parent
DESIGN = HERE/'native-precision-review-consumer-design.json'
THERMAL_VALIDATION = ROOT/'studies/unified_cosmology/results/inference/native-precision-thermal-review-validation.json'


def finite_tree(value):
    """Do not let a NaN in an unused report field escape identity replay."""
    if isinstance(value, dict):
        for child in value.values(): finite_tree(child)
    elif isinstance(value, list):
        for child in value: finite_tree(child)
    elif isinstance(value, (int, float)):
        assert np.isfinite(value), 'Nonfinite numerical report value.'


def read(path):
    value = json.loads(Path(path).read_text())
    finite_tree(value)
    return value


def hashes(mapping):
    assert isinstance(mapping, dict)
    for name, checksum in mapping.items():
        assert original.digest(ROOT/name) == checksum, 'Changed input: '+name


def merge(*mappings):
    out = {}
    for mapping in mappings:
        for name, checksum in mapping.items():
            assert name not in out or out[name] == checksum, 'Conflicting input identity.'
            out[name] = checksum
    return out


def _verify(screen_path, review_path):
    screen_path, review_path = Path(screen_path).resolve(), Path(review_path).resolve()
    assert screen_path != review_path
    screen, report = read(screen_path), read(review_path)
    assert report['status'] == 'completed_auxiliary_thermal_path_review'
    assert report['posterior_qualification'] is False and report['CMB_spectrum_calls'] == 0
    assert report['background_calls'] == 64
    assert report['original_screen_path'] == original.relative(screen_path)
    assert report['original_screen_sha256'] == original.digest(screen_path)
    plan_path = ROOT/screen['plan_path']
    assert original.digest(plan_path) == screen['plan_sha256']
    plan = original.sealed_read(plan_path); finite_tree(plan)
    assert plan['identity'] == original.identity({k:v for k,v in plan.items() if k not in ['identity','payload_sha256']})
    assert plan['design'] == read(original.DESIGN)
    assert len(plan['selected']) == plan['design']['points'] == 32
    correction_path = ROOT/plan['correction_summary_path']
    assert original.digest(correction_path) == plan['correction_summary_sha256']
    correction = read(correction_path)
    selection_path = ROOT/correction['selection_path']
    selection = read(selection_path)
    folder = selection_path.parent.parent
    # The real qualifier re-identifies every target asset and all native/chain
    # evidence. No synthetic bypass or external qualified=True flag is accepted.
    parent = summarize_run(folder, correction_path)
    assert parent['qualified_under_declared_numerical_gates'] is True
    assert parent['target_identity'] == plan['frozen_target']['identity'] == screen['frozen_target_identity']
    assert parent['input_sha256'] == plan['qualified_parent_inputs']
    assert read(folder/'run-0.json')['target_identity'] == plan['frozen_target']
    assert selection['settings'] == plan['settings']
    original.verify_bindings(plan)
    expected_source_names = {original.relative(p) for p in [Path(original.__file__), original.DESIGN,
        HERE/'measurement_summary.py', HERE/'exact_correction.py', HERE/'native_precision_audit.py',
        HERE/'likelihood.py', HERE/'target_identity.py']}
    assert set(plan['source_sha256']) == expected_source_names
    assert screen['source_sha256'] == plan['source_sha256']
    assert screen['declared_target_changed'] is False
    assert screen['weight_context_only'] == plan['weight_context_only']
    assert screen['numerical_controls'] == plan['design']['numerical_controls']
    nominal, high = original.native_configurations(plan['settings'], plan['design'])
    assert canonical(nominal) == plan['nominal_native_configuration']
    assert canonical(high) == plan['doubled_native_configuration']
    # Reconstruct chronological strata from the qualified parent, never from
    # the saved screen's likelihood outcomes.
    selected = original.select_indices(selection['groups'], selection['locations'], plan['design']['points_per_chain'])
    native_records = []
    for index, point in enumerate(selection['points']):
        path = selection_path.parent/f'{index:05d}.json'
        row = read(path); verify_record(row)
        assert row['status'] == 'finite' and row['point'] == point
        native_records.append(row)
    lw = np.array([r['log_weight'] for r in native_records]); weights = np.exp(lw-logsumexp(lw))
    for rebuilt, saved in zip(selected, plan['selected']):
        i = rebuilt['parent_index']; path = selection_path.parent/f'{i:05d}.json'
        rebuilt.update(point=selection['points'][i], native_record_path=original.relative(path),
            native_record_sha256=original.digest(path), original_location=selection['locations'][i],
            original_native_logweight=float(lw[i]), weight_in_full_qualified_parent=float(weights[i]))
        assert rebuilt == saved, 'Precision selection differs from fixed chronological strata.'
    expected_context = {'selected_weight_mass_in_full_parent': float(sum(x['weight_in_full_qualified_parent'] for x in selected)),
        'full_parent_chain_weight_masses': {str(c): float(weights[np.asarray(selection['groups']) == c].sum()) for c in range(4)},
        'use': 'Context only. Selection and numerical-screen centering are unweighted; no posterior reweighting.'}
    assert expected_context == plan['weight_context_only']
    manifest_path = ROOT/screen['record_manifest_path']
    assert original.digest(manifest_path) == screen['record_manifest_sha256']
    manifest = read(manifest_path); hashes(manifest)
    assert len(screen['points']) == len(report['point_reviews']) == 32
    expected_manifest = {}; adjusted_rows = []; comparisons = []
    review_design = read(thermal.DESIGN)
    assert report['interpretation'] == review_design['interpretation']
    for index, (row, selected, saved_review) in enumerate(zip(screen['points'], plan['selected'], report['point_reviews'])):
        assert row['status'] in {'finite_native_precision', 'failed_native_precision_checks'}, 'Incomplete native point cannot be repaired.'
        assert row['audit_index'] == index and row['native_point_evaluations'] == 1
        assert row['identity'] == plan['identity'] and row['plan_sha256'] == original.digest(plan_path)
        for key in ['point','parent_index','chain','native_record_sha256','original_native_logweight','weight_in_full_qualified_parent']:
            assert row[key] == selected[key]
        path = ROOT/row['record_path']; sealed = original.sealed_read(path); finite_tree(sealed)
        augmented = dict(sealed)
        for kind in ['record','log','spectrum']:
            name, checksum = row[kind+'_path'], row[kind+'_sha256']
            assert original.digest(ROOT/name) == checksum
            expected_manifest[name] = checksum
            if kind != 'spectrum': augmented.update({kind+'_path':name,kind+'_sha256':checksum})
        augmented['integrated_time_mismatch_warning_count'] = (ROOT/row['log_path']).read_text().count('mismatch in integrated times')
        assert row == augmented, 'Screen row does not reproduce sealed native output and log.'
        stored = native_records[selected['parent_index']]
        assert row['stored_declared_loglikes'] == stored['exact_loglikes']
        assert row['stored_declared_logpost'] == stored['exact_logpost']
        adjusted, rebuilt = thermal.reviewed_row(row, stored, saved_review['backgrounds'], plan['design'], review_design)
        rebuilt.update(audit_index=index, original_record_path=original.relative(path), original_record_sha256=original.digest(path))
        assert rebuilt == saved_review, 'Supplemental point review does not replay.'
        adjusted_rows.append(adjusted)
        comparisons.append({'audit_index': index, 'original_status': row['status'], 'original_failed_checks': row['failed_checks'],
                            'reviewed_status': adjusted['status'], 'reviewed_failed_checks': adjusted['failed_checks']})
    assert manifest == expected_manifest
    assert screen['known_native_point_evaluations'] == 32 and screen['incomplete_attempts_with_unknown_native_call_count'] == 0
    original_summary = original.summarize(screen['points'], plan['design'])
    assert original_summary == {key:screen[key] for key in original_summary}
    reviewed_summary = original.summarize(adjusted_rows, plan['design'])
    assert reviewed_summary == report['reviewed_screen']
    assert report['original_screen_status'] == original_summary['status']
    # Exactly reconstruct the supplemental producer's binding map; mere hash
    # verification of an incomplete user-supplied dictionary is insufficient.
    producer_paths = [screen_path, plan_path, correction_path, manifest_path, Path(thermal.__file__),
                      thermal.DESIGN, HERE/'expansion_history.py', Path(original.__file__)]
    expected_bindings = merge(parent['input_sha256'], manifest,
                             {original.relative(p):original.digest(p) for p in producer_paths})
    assert report['input_source_sha256'] == expected_bindings
    hashes(expected_bindings)
    validation = read(THERMAL_VALIDATION)
    assert validation['status'] == 'passed_background_only_thermal_review_validation'
    hashes(validation['source_sha256']); hashes(validation['input_sha256'])
    own = [Path(__file__), DESIGN, review_path, THERMAL_VALIDATION]
    inputs = merge(expected_bindings, plan['source_sha256'], validation['source_sha256'], validation['input_sha256'],
                   {original.relative(p):original.digest(p) for p in own})
    hashes(inputs)
    return {'schema':'verified-supplemental-native-precision-v1', 'status':'verified_saved_precision_review',
        'parent_target_identity': parent['target_identity'], 'parent_freshly_qualified': True,
        'original_screen_status': original_summary['status'], 'reviewed_screen_status': reviewed_summary['status'],
        'original_diagnostic_flags': original_summary.get('diagnostic_flags', []),
        'reviewed_diagnostic_flags': reviewed_summary.get('diagnostic_flags', []),
        'original_failures': original_summary.get('failures', []), 'reviewed_failures': reviewed_summary.get('failures', []),
        'reviewed_precision_screen_supported': reviewed_summary['status'] == 'no_large_variation_detected_on_fixed32',
        'original_screen': original_summary, 'reviewed_screen': reviewed_summary, 'point_statuses': comparisons,
        'original_screen_ref': {'path':original.relative(screen_path),'sha256':original.digest(screen_path)},
        'supplemental_review_ref': {'path':original.relative(review_path),'sha256':original.digest(review_path)},
        'input_sha256': inputs, 'posterior_qualification':False, 'posterior_reweighting_performed':False,
        'background_or_native_or_model_calls':0,
        'limits':read(DESIGN)['limits']}


def verify(screen_path, review_path):
    """Freshly qualify the parent and replay stored evidence, with theory forbidden."""
    import camb
    import cobaya.model
    def forbidden(*args, **kwargs):
        raise RuntimeError('This consumer may not construct models or calculate backgrounds/spectra.')
    with ExitStack() as stack:
        for owner, name in [(camb,'get_background'),(camb,'get_results'),(camb,'get_transfer_functions'),
                            (cobaya.model,'get_model'),(cobaya.model.Model,'__init__'),(thermal,'backgrounds')]:
            stack.enter_context(patch.object(owner,name,forbidden))
        return _verify(screen_path, review_path)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--screen',type=Path,required=True); parser.add_argument('--review',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(); assert not args.output.exists(), 'Fresh supplemental receipt only.'
    result=verify(args.screen,args.review)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['status','original_screen_status','reviewed_screen_status','reviewed_precision_screen_supported']}))


if __name__=='__main__': main()
