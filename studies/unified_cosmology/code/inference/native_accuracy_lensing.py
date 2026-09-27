"""Array-only Qu lensing replacements of a freshly qualified native2 target."""
from __future__ import annotations
import argparse
import copy
import importlib.metadata
import json
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

import joint_lensing_bridge as bridge
import native_accuracy_runtime as rt
import spectral_correction as capture
from exact_correction import verify_record
from measurement_summary import weighted_fraction
from native_accuracy_postprocess import chronology
from target_identity import canonical

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DESIGN = HERE/'native-accuracy-lensing-design.json'
VALIDATION = ROOT/'studies/unified_cosmology/results/inference/native-accuracy-lensing-validation.json'


def sources():
    paths = [Path(__file__), DESIGN, HERE/'native_accuracy_lensing_validate.py',
             HERE/'native_accuracy_qualify.py', HERE/'native_accuracy_runtime.py',
             HERE/'native_accuracy_postprocess.py', HERE/'joint_lensing_bridge.py',
             bridge.DESIGN, bridge.GATES, bridge.AUDIT, HERE/'spectral_correction.py',
             HERE/'exact_correction.py', HERE/'measurement_summary.py',
             HERE/'luminosity_sensitivity.py', HERE/'target_identity.py',
             HERE/'native_posterior_precision.py', bridge.EXTERNAL/'fast_lensing.py']
    return {rt.relative(p):rt.digest(p) for p in paths}


def source_constants(record, old_factors):
    """Metadata-only replay; never read an accuracy1 spectrum for new lensing."""
    path = ROOT/record['parent_native_path']
    assert rt.digest(path) == record['parent_native_sha256']
    parent = json.loads(path.read_text()); verify_record(parent)
    assert parent['status'] == 'finite'
    assert parent['index'] == record['index'] and parent['point'] == record['point']
    assert parent['exact_logpost'] == record['native1_logpost']
    assert parent['exact_loglikes'] == record['native1_loglikes']
    assert parent['log_weight'] == record['native1_logweight']
    assert parent['derived'] == record['native1_derived']
    sidecar_path = capture.sidecar_path(path)
    row = json.loads(sidecar_path.read_text()); verify_record(row)
    assert row['status'] == 'captured_from_parent_native_evaluation'
    assert row['parent_native_record_sha256'] == record['parent_native_sha256']
    assert row['parent_payload_sha256'] == parent['payload_sha256']
    assert row['index'] == parent['index'] and row['point'] == parent['point']
    assert row['correction_identity'] == parent['target_identity']
    assert row['source_sha256'] == capture.dependencies()
    assert row['CLIPY_NOJAX'] == '1' and row['additional_native_evaluations'] == 0
    assert row['spectral_units'] == 'raw C_l: TT/EE/TE/BB FIRASmuK2, pp dimensionless potential'
    constants = row['source_normalized_gaussian_constants']
    assert set(constants) == set(old_factors) and np.isfinite(list(constants.values())).all()
    hashes = rt.merge(row['source_sha256'], {rt.relative(path):record['parent_native_sha256'],
                         rt.relative(sidecar_path):rt.digest(sidecar_path)})
    return constants, hashes


def configuration(target, inherited):
    config = copy.deepcopy(target.native_configuration)
    assert set(config['theory']) == {'camb'}
    assert all(config['theory']['camb']['extra_args'][k] == 2
               for k in ('AccuracyBoost','lAccuracyBoost','lSampleBoost'))
    assert all(k in config['likelihood'] for k in inherited['old_factors'])
    assert config['params']['A_fg']['prior'] == inherited['auxiliary_parameter']['normalized_prior']
    act = config['likelihood'][inherited['old_factors'][0]]
    assert act['variant'] == 'actplanck_baseline' and not act['lens_only']
    assert act['apply_hartlap'] and act['nsims_planck'] == 400
    assert config['likelihood']['SPT2023_lensing']['clear_internal_priors']
    return config


def density(record, new, constants, old_factors, expected_components, tolerance):
    """Exact replacement algebra from native2, not a fabricated old record."""
    likes = record['native2_loglikes']
    assert set(likes) == set(expected_components)
    assert np.isfinite([*likes.values(), record['native2_logpost'],
                       record['proposal_logpost'], record['log_weight']]).all()
    assert abs(record['native2_logpost']-record['proposal_logpost']-record['log_weight']) <= tolerance
    assert abs(record['native2_logpost']-sum(likes.values())-sum(record['native2_logpriors'])) <= tolerance
    # Pure bridge kernel's documented arithmetic view; never serialized as an
    # old exact-correction record or supplied to its filesystem consumer.
    value = bridge.bridge_density({'exact_loglikes':likes, 'log_weight':record['log_weight']},
                                  new, constants, old_factors)
    target_post = record['native2_logpost']+value['lens_kernel_logratio']
    closure = target_post-record['proposal_logpost']-value['target_logweight']
    assert np.isfinite([*value.values(), target_post, closure]).all() and abs(closure) <= tolerance
    value.update(target_logpost_in_inherited_normalization=float(target_post),
                 target_weight_accounting_closure=float(closure))
    return value


def summarize(records, rows, groups, gates):
    view = [{'point':r['point'], 'derived':r['native2_derived'], 'log_weight':r['log_weight']} for r in records]
    result = bridge.summarize_variant(view, rows, groups, gates)
    qualified = result.get('qualified_under_declared_numerical_gates', False)
    if not qualified:
        result.update(posterior=None, weighted_covariance=None, conditional_sign_fractions=None)
        return result
    weights = np.exp(np.array([r['density']['target_logweight'] for r in rows])-logsumexp(
        [r['density']['target_logweight'] for r in rows]))
    values = {k:np.array([r['native2_derived'][k] for r in records]) for k in ('q0','q05','q1','j0')}
    signs = {k+'_below_zero':weighted_fraction(v < 0, weights, groups) for k,v in values.items()}
    signs['accelerating_with_decreasing_scale_factor_acceleration_today'] = weighted_fraction(
        (values['q0'] < 0) & (values['j0'] < 0), weights, groups)
    result['conditional_sign_fractions'] = signs
    assert 'A_fg' not in result['posterior'] and 'A_fg' not in result['weighted_covariance']['parameter_order']
    return result


def completion_guard(inputs, generated, source, versions, validation_sha):
    rt.verify_hashes(inputs); rt.verify_hashes(generated)
    assert rt.digest(VALIDATION) == validation_sha and sources() == source
    assert all(importlib.metadata.version(k) == v for k,v in versions.items()), 'Frozen environment changed.'


def actual(work, summary_path, cache):
    from native_accuracy_qualify import QualifiedNativeAccuracy2, qualify_run
    source = sources(); validation_sha = rt.digest(VALIDATION)
    validation = json.loads(VALIDATION.read_text())
    assert validation['status'] == 'passed_native_accuracy2_lensing_validation'
    assert validation['source_sha256'] == source and validation['physical_calls'] == 0
    # Runtime guard also prevents an accidental future physical route.
    with rt.guards_without_physics():
        target = qualify_run(work, summary_path)
        assert isinstance(target, QualifiedNativeAccuracy2) and target.numerical_accuracy == 2
        rt.verify_hashes(target.input_sha256)
        versions = dict(target.plan['evidence']['versions'])
        assert versions and all(importlib.metadata.version(k) == v for k,v in versions.items())
        design = json.loads(DESIGN.read_text()); inherited = json.loads(bridge.DESIGN.read_text())
        assert design['variants'] == inherited['variants']
        config = configuration(target, inherited)
        records = target.records; groups = np.asarray(target.groups)
        assert len(records) == len(target.selected_points) == 2000
        chronology_result = chronology(target.selected_points, groups, target.locations, target.proposal_logposts)
        constants = None; inputs = dict(target.input_sha256)
        for i, record in enumerate(records):
            assert record['index'] == i and record['point'] == target.selected_points[i]
            assert record['group'] == groups[i] and record['status'] == 'finite_native_accuracy2'
            assert record['numerical_target_identity'] == target.numerical_target_identity
            assert record['log_weight'] == target.logweights[i] and record['proposal_logpost'] == target.proposal_logposts[i]
            assert record['spectrum_path'] in target.input_sha256
            assert target.input_sha256[record['spectrum_path']] == record['spectrum_sha256']
            current, hashes = source_constants(record, inherited['old_factors'])
            inputs = rt.merge(inputs, hashes)
            if constants is None: constants = current
            else: assert constants == current, 'Source Gaussian constants changed between slots.'
        audit = json.loads(bridge.AUDIT.read_text())
        assert audit['status'] == 'released_likelihood_evaluation_reproduced_no_inference'
        assets = rt.merge(audit['input_sha256'], audit['source_sha256'])
        inputs = rt.merge(inputs, assets, {rt.relative(VALIDATION):validation_sha})
        rt.verify_hashes(inputs)
        descriptions = bridge.target_descriptions(config, inherited, assets)
        plan = {'schema':'native-accuracy2-joint-lensing-plan-v1', 'design':design,
            'parent_numerical_target_identity':target.numerical_target_identity,
            'parent_proposal_target_identity':target.parent_proposal_target_identity,
            'parent_settings':target.parent_settings, 'native_configuration':canonical(config),
            'scientific_versions':versions,
            'variant_target_descriptions':descriptions, 'source_sha256':source, 'input_sha256':inputs,
            'source_normalized_gaussian_constants':constants,
            'parent_summary_path':rt.relative(summary_path), 'parent_summary_sha256':rt.digest(summary_path),
            'chronology':chronology_result, 'slots':2000}
        plan['identity'] = rt.identity(plan)
        cache = Path(cache).resolve()
        assert cache.is_relative_to(ROOT/'.work') and not cache.exists(), 'Use a fresh ignored cache.'
        cache.mkdir(parents=True)
        rt.write_new(cache/'plan.json', plan)
        results = {}; generated = {rt.relative(cache/'plan.json'):rt.digest(cache/'plan.json')}
        for variant in design['variants']:
            definition = {'numerical_accuracy':2, 'parent_numerical_target_identity':target.numerical_target_identity,
                          'released_target':descriptions[variant], 'source_sha256':source}
            target_id = rt.identity(definition)
            response = bridge.JointResponse(variant)
            directory = cache/variant; directory.mkdir(); rows = []; hashes = {}
            for i, record in enumerate(records):
                try:
                    path = ROOT/record['spectrum_path']
                    assert rt.digest(path) == record['spectrum_sha256']
                    with np.load(path, allow_pickle=False) as spectra:
                        capture.validate_spectra(spectra)
                        new = response.evaluate(spectra)
                    value = density(record, new, constants, inherited['old_factors'], config['likelihood'],
                                    design['density_absolute_tolerance'])
                    row = {'status':'finite_bridge', 'density':value, 'new_lensing':new}
                except Exception as error:
                    row = {'status':'failed_bridge', 'error':repr(error)}
                row.update(index=i, point=record['point'], numerical_target_identity=target_id,
                    parent_numerical_target_identity=target.numerical_target_identity,
                    spectrum_path=record['spectrum_path'], spectrum_sha256=record['spectrum_sha256'],
                    native_CAMB_calls=0, background_calls=0)
                path = directory/f'{i:05d}.json'; rt.write_new(path, row)
                hashes[rt.relative(path)] = rt.digest(path); rows.append(row)
            manifest = directory/'records.json'
            rt.write_new(manifest, {'numerical_target_identity':target_id, 'file_sha256':hashes})
            generated = rt.merge(generated, hashes, {rt.relative(manifest):rt.digest(manifest)})
            result = summarize(records, rows, groups, json.loads(bridge.GATES.read_text())['overlap_gates'])
            result.update(numerical_accuracy=2, numerical_target_identity=target_id,
                          target_definition=definition, records_path=rt.relative(manifest), records_sha256=rt.digest(manifest))
            results[variant] = result
            del response
        completion_guard(inputs, generated, source, versions, validation_sha)
        return {'schema':'native-accuracy2-joint-lensing-result-v1',
            'status':'separate_native_accuracy2_lensing_sensitivities', 'numerical_accuracy':2,
            'parent_numerical_target_identity':target.numerical_target_identity,
            'parent_proposal_target_identity':target.parent_proposal_target_identity,
            'settings':{k:target.parent_settings[k] for k in ('model','evolution','sample','calibration')},
            'variants':results, 'native_CAMB_calls':0, 'background_calls':0, 'model_construction_calls':0,
            'slots':2000, 'chronology':chronology_result, 'source_sha256':source,
            'input_sha256':inputs, 'generated_sha256':generated, 'plan_path':rt.relative(cache/'plan.json'),
            'plan_sha256':rt.digest(cache/'plan.json'), 'limitations':design['limitations']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('work','summary','cache','output'): parser.add_argument('--'+key, type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists(), 'Preserve all prior results.'
    result = actual(args.work.resolve(), args.summary.resolve(), args.cache)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rt.write_new(args.output, result)
    print(json.dumps({'status':result['status'], 'variants':{k:v['status'] for k,v in result['variants'].items()},
                      'native_CAMB_calls':0}))


if __name__ == '__main__': main()
