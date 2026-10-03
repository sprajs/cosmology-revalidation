#!/usr/bin/env python3
"""One bounded released active46 box median comparison; observational use blocked.

Shared admission/receipt and structural FITS helpers retain their own historical
requests. This distinct controller performs no production Gaussian mathematics.
"""
import argparse
from datetime import datetime, timezone
from decimal import Decimal, localcontext
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import shutil
import subprocess
import sys
import tarfile
import time

ROOT = Path(__file__).resolve().parents[2]
FOLDER = Path(__file__).resolve().parent
REQUEST_SHA256 = '664b7054f1d2c3d23ae7af170c163e8d717213e2ccf05e5945e1971b56635de8'
CONTRACT_SHA256 = 'ee885213f48cdb04d3b6d0bada41fea3d65524feb9f33b0d3bda985645aa0741'
ACTIVE = list(range(44)) + [45, 46]
COMPONENT_SIZES = [2593, 55, 339, 143, 354] + [1] * 8
SHARED_PATHS = ('experiments/lcdm-campaign/controller.py',
                'experiments/released-ladder/controller.py', 'scripts/packet.py')
LOADED_SOURCE_HASHES = {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                        for name in SHARED_PATHS + (str(Path(__file__).resolve().relative_to(ROOT)),)}
TREE_BYTES = 2097152
STORE_BYTES = 268435456
PI = Decimal('3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679')


def module(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


receipt = module('released_box_receipt', SHARED_PATHS[0])
lineage = module('released_box_lineage', SHARED_PATHS[1])
parse, load, sha256, within = receipt.parse, receipt.load, receipt.sha256, receipt.within
keys, exact, finite, count = receipt.keys, receipt.exact, receipt.finite, receipt.count
decimal_string = receipt.decimal_string
ComparisonFailure = receipt.ComparisonFailure

NATIVE_FIELDS = ('schema_version interface_id contract_sha256 request_sha256 sdk_identity target '
                 'source_identities arithmetic policy resources stages payload_bounds result_status '
                 'result_numerical_status method_id enclosure_scope result_stage availability completion '
                 'box_diagnostics normalizations median work output_complete accepted qualification').split()
FLAGS = ('gaussian_completion endpoint_margins rectangle_enclosure normalization_enclosures '
         'quantile_enclosure endpoint_cdf_enclosures').split()
COMPLETION = ('mean variance named_variance_original46 minimum_quadratic '
              'reported_postcast_profile_quadratic log_design_precision_determinant '
              'source_covariance_log_determinant maximum_variance_sensitivity profile_stationarity').split()
DIAGNOSTICS = ('standardized_lower_margins standardized_upper_margins excluded_mass_upper '
               'box_probability log_box_probability log_prior_volume').split()
NORMALIZATIONS = ('log_relative_box_integral log_prior_normalized_relative_evidence '
                  'log_observation_normalized_evidence').split()
STAGES = ('last_stage error input_complete gaussian_status gaussian_numerical_status design_status '
          'design_numerical_status design_rank equilibrated_triangular_condition_inf '
          'equilibrated_transpose_triangular_condition_inf design_method box_preparation_status '
          'box_preparation_numerical_status box_preparation_numerical_status_scope design_consumed_by_box').split()
REFERENCE_FIELDS = ('schema_version interface_id status accepted error blockers identities target settings runtime_before '
                    'runtime_after source_structure diagnostics coarse fine refinement gates work qualification').split()
REFERENCE_SCALARS = ('q_min q_postcast q_postcast_unrounded_residual q_postcast_residual_rounding_allowance '
                     'logdet_H logdet_C excluded_mass_upper reduced_forward_sensitivity').split()
REFERENCE_VECTORS = ('mean variance endpoint_lower_margins endpoint_upper_margins').split()
REFERENCE_INTERVALS = ('box_probability log_box_probability log_prior_volume log_Z_relative_box '
                       'log_prior_normalized_relative_evidence log_observation_normalized_evidence median_original46').split()
REFERENCE_GATES = ('source_identity runtime_identity original_input_empirical_completion original_input_empirical_covariance_logdet '
                   'finite_box_tail reference_refinement original_input_empirical_normalization original46_median_reference').split()
REFINEMENT_SCALARS = ('q_min q_postcast logdet_H logdet_C log_prior_volume log_Z_relative_box '
                      'log_prior_normalized_relative_evidence log_observation_normalized_evidence median_original46').split()
ORIGINAL_SOURCE_SHA256 = {
    'MCMC_utils.py': 'e6ec3d83a9b126d7772ec6dd0d1b58ca757b820acd853fb89218c379cf873841',
    'lstsq_results.txt': '37d2d423d06b6a2100c47578eb9b1c566a575caf6ade0ed61f6b9586d2a35c95',
    'allc_shoes_ceph_topantheonwt6.0_112221.fits': 'a42778672d25df7a559bd2949b1e412b99b2e35dd9a36895f6b38828be019172',
    'alll_shoes_ceph_topantheonwt6.0_112221.fits': '9a2ce872dd20ed4fdf5005ce62a805d3eefbbb0b20b0c614013e7c0094009db8',
    'ally_shoes_ceph_topantheonwt6.0_112221.fits': '10bb034ca3fff53f6625c9fe47ebd054b37ef45d1140af1c766831cfa6871433',
    'run_mcmc.py': '00546b9de2d73cd6315c0333b5f7b07a8b948f9e73056a88e76abe0249464c85',
}
REFERENCE_SETTINGS = {
    'digits': [80, 120], 'correction_steps': [1, 3], 'diagnostic_digits': 36,
    'postcast_source_residual_digits': 2048, 'postcast_correction_steps': 3,
    'observations': 3492, 'original_columns': 47, 'active_columns': 46,
    'covariance_route': 'original-C-LAPACK-Cholesky-with-longdouble-residual-refinement/v1',
    'completion_peers': ['LAPACK-gesdd-SVD', 'LAPACK-pivoted-economic-QR'],
    'logdet_C_peer': 'original-C-component-LAPACK-eigvalsh-evr/v1',
    'reduced_route': 'stdlib-Decimal-LDLT-original-X-transpose-Cinv-X/v1',
    'tail_route': 'correlation-valid-union-Chernoff-Decimal/v1',
    'median_route': 'correlation-valid-union-bound-original46/v1', 'reference_refinement_fraction': '0.05',
    'allocations': {'coefficient_absolute': '1e-8', 'coefficient_relative': '1e-9', 'variance_absolute': '2e-12',
                    'variance_relative': '2e-12', 'quadratic_absolute': '1e-7', 'quadratic_relative': '1e-10',
                    'log_absolute': '1e-7', 'median_absolute': '5e-10'},
    'resources': {'jobs': 1, 'threads': 1, 'cpu_seconds': 900, 'wall_seconds': 900,
                  'address_limit_bytes': 3221225472, 'output_limit_bytes': 1048576},
    'maximum_excluded_mass': '0.1', 'maximum_covariance_forward_sensitivity': '1e-10',
    'minimum_longdouble_mantissa_bits': 64,
}


def reference_vector(value, label, positive=False, length=46):
    if not isinstance(value, list) or len(value) != length:
        raise ValueError(label + ' ordered length differs')
    result = [decimal_string(v, label) for v in value]
    if positive and any(v <= 0 for v in result):
        raise ValueError(label + ' requires positive values')
    return result


def reference_target(request):
    return {'active_original_indices': ACTIVE, 'fixed_coordinate': {'original_index': 44, 'value': 0.0},
            'marginal': {'original_index': 46, 'active_index': 45, 'cumulative_probability': 0.5},
            'lower_binary64_hex': [row['lower_binary64_hex'] for row in request['literal_support']],
            'upper_binary64_hex': [row['upper_binary64_hex'] for row in request['literal_support']],
            'parameter_measure': 'original-active46-Lebesgue-fixed44-point-mass-zero/v1'}


def reference_residual(record, steps, length=47, digits=None):
    array_keys = ('residual_infinity_by_rhs', 'backward_error_by_rhs', 'estimated_forward_error_by_rhs')
    fields = ['correction_steps', 'inverse_infinity_norm_estimate', 'component_estimates'] + list(array_keys)
    if digits is not None:
        fields.append('completion_digits')
    keys(record, fields, 'source residual diagnostic')
    exact(record['correction_steps'], steps, 'ordered source correction steps')
    if digits is not None:
        exact(record['completion_digits'], digits, 'ordered postcast completion precision')
    values = {key: reference_vector(record[key], key, length=length) for key in array_keys}
    if any(value < 0 for array in values.values() for value in array):
        raise ValueError('negative source residual diagnostic')
    components = record['component_estimates']
    if not isinstance(components, list) or len(components) != len(COMPONENT_SIZES):
        raise ValueError('source residual component coverage differs')
    component_values, inverse_norms = [], []
    for index, (component, rows) in enumerate(zip(components, COMPONENT_SIZES)):
        keys(component, ['component_index', 'rows', 'inverse_infinity_norm_estimate'] + list(array_keys), 'component residual')
        exact(component['component_index'], index, 'ordered source component')
        exact(component['rows'], rows, 'source component row count')
        arrays = {key: reference_vector(component[key], key, length=length) for key in array_keys}
        if any(value < 0 for array in arrays.values() for value in array):
            raise ValueError('negative source component residual diagnostic')
        norm = decimal_string(component['inverse_infinity_norm_estimate'], 'component inverse norm')
        if norm <= 0:
            raise ValueError('component inverse norm domain differs')
        component_values.append(arrays)
        inverse_norms.append(norm)
    for key in array_keys:
        if values[key] != [max(component[key][j] for component in component_values) for j in range(length)]:
            raise ValueError('source residual diagnostic differs from earned component maxima')
    if decimal_string(record['inverse_infinity_norm_estimate'], 'inverse norm') != max(inverse_norms):
        raise ValueError('source residual inverse norm differs from component maxima')


def check_reference(output, request, identity, require_complete=True):
    keys(output, REFERENCE_FIELDS, 'reference')
    exact(output['schema_version'], 1, 'reference schema')
    exact(output['interface_id'], 'released-fixed44-box-reference/v1', 'reference interface')
    if type(output['accepted']) is not bool:
        raise ValueError('reference accepted requires boolean')
    complete = output['accepted']
    exact(output['status'], 'accepted' if complete else 'refused', 'reference status')
    if complete:
        exact(output['error'], None, 'reference error')
        exact(output['blockers'], [], 'reference blockers')
    else:
        if output['error'] is not None and (not isinstance(output['error'], str) or not 1 <= len(output['error']) <= 4096):
            raise ValueError('reference refusal error domain differs')
        if not isinstance(output['blockers'], list) or not 1 <= len(output['blockers']) <= 64 or any(
                not isinstance(value, str) or not 1 <= len(value) <= 4096 for value in output['blockers']):
            raise ValueError('reference refusal causes missing or malformed')
    expected = {key: identity[key] for key in ('request_sha256', 'reference_input_sha256', 'reference_script_sha256')}
    expected.update({'contract_sha256': CONTRACT_SHA256, 'lineage_sha256': request['lineage']['packet_lineage_sha256'],
                     'constrained_sha256': request['lineage']['constrained_sha256'], 'original_source_sha256': ORIGINAL_SOURCE_SHA256,
                     **{'canonical_' + item['name'][0] + '_sha256': item['sha256'] for item in request['transport']['canonical_inputs']}})
    identities = output['identities']
    if complete:
        exact(identities, expected, 'reference identities')
    else:
        if not isinstance(identities, dict) or not set(identities).issubset(expected):
            raise ValueError('partial reference identity fields differ')
        for key, value in identities.items():
            if key == 'original_source_sha256':
                if not isinstance(value, dict) or not set(value).issubset(ORIGINAL_SOURCE_SHA256):
                    raise ValueError('partial reference source identities differ')
                for name, digest in value.items():
                    exact(digest, ORIGINAL_SOURCE_SHA256[name], 'partial reference source/' + name)
            else:
                exact(value, expected[key], 'partial reference identity/' + key)
    if complete or output['target'] is not None:
        exact(output['target'], reference_target(request), 'reference target')
    exact(output['settings'], REFERENCE_SETTINGS, 'reference settings')
    for key in ('runtime_before', 'runtime_after'):
        if complete or output[key] is not None:
            exact(output[key], identity['runtime'], 'reference/' + key)
    if complete or output['source_structure'] is not None:
        exact(output['source_structure'], {'inventory_sha256': 'e41d7f45428056e8c05e2952dbac54020009422e9a9e4c583aa4158e7eb2c2b1',
                                      'component_count': 13, 'maximum_component_size': 2593, 'maximum_general_component_size': 2593,
                                      'component_sizes': COMPONENT_SIZES,
                                      'dense_factor_multiply_add_pairs_per_precision': 2913157849}, 'original covariance structure')
    qualification = output['qualification']
    keys(qualification, ('scope', 'original_input_certificate', 'conditional_enclosure_scope', 'native_outputs_consumed',
                         'observational_qualification', 'limits'), 'reference qualification')
    for key, value in (('scope', 'empirical-original-input-comparison/v1'), ('original_input_certificate', False),
                       ('conditional_enclosure_scope', 'conditional-Gaussian-union-mass-and-median-only/v1'),
                       ('native_outputs_consumed', False), ('observational_qualification', 'blocked-original-contract-source-gaps')):
        exact(qualification[key], value, 'reference qualification/' + key)
    if not isinstance(qualification['limits'], list) or not qualification['limits'] or any(
            not isinstance(v, str) or not 1 <= len(v) <= 4096 for v in qualification['limits']):
        raise ValueError('reference limitations missing or malformed')
    keys(output['gates'], REFERENCE_GATES, 'reference gates')
    if any(type(value) is not bool for value in output['gates'].values()):
        raise ValueError('reference gates require booleans')
    if complete:
        exact(output['gates'], dict.fromkeys(REFERENCE_GATES, True), 'reference gates')
    work = output['work']
    keys(work, ('elapsed_seconds', 'cpu_user_seconds', 'cpu_system_seconds', 'maximum_rss_kib', 'covariance_factorizations',
                'covariance_eigenvalue_components', 'wide_residual_evaluations', 'correction_steps_completed',
                'reduced_decimal_factorizations', 'postcast_covariance_solves', 'postcast_wide_residual_evaluations',
                'postcast_correction_steps_completed', 'original_input_rows', 'native_outputs_consumed'), 'reference work')
    for key in ('elapsed_seconds', 'cpu_user_seconds', 'cpu_system_seconds'):
        if finite(work[key], key) < 0 or finite(work[key], key) > 901:
            raise ValueError('reference time budget exceeded')
    if finite(work['cpu_user_seconds'], 'CPU user') + finite(work['cpu_system_seconds'], 'CPU system') > 901:
        raise ValueError('reference total CPU budget exceeded')
    for key in ('maximum_rss_kib', 'covariance_factorizations', 'covariance_eigenvalue_components', 'wide_residual_evaluations',
                'correction_steps_completed', 'reduced_decimal_factorizations'):
        count(work[key], sys.maxsize, key)
    for key, value in (('covariance_factorizations', 1), ('covariance_eigenvalue_components', 13),
                       ('wide_residual_evaluations', 4), ('correction_steps_completed', 3), ('reduced_decimal_factorizations', 2),
                       ('postcast_covariance_solves', 2), ('postcast_wide_residual_evaluations', 8),
                       ('postcast_correction_steps_completed', 6)):
        if complete:
            exact(work[key], value, 'successful reference work/' + key)
        else:
            count(work[key], value, 'partial reference work/' + key)
    if work['maximum_rss_kib'] * 1024 > request['resources']['address_limit_bytes']:
        raise ValueError('reference memory budget exceeded')
    exact(work['original_input_rows'], 3492, 'reference rows')
    exact(work['native_outputs_consumed'], False, 'reference independence')
    diagnostics = output['diagnostics']
    diagnostic_fields = ('ancestry', 'original_input_LAPACK_SVD', 'original_input_LAPACK_pivoted_QR',
                         'covariance_condition_estimate', 'wide_residual_history', 'postcast_residual_history', 'logdet_C_peer')
    if complete:
        keys(diagnostics, diagnostic_fields, 'reference diagnostics')
    elif (not isinstance(diagnostics, dict) or not {'ancestry', 'wide_residual_history', 'postcast_residual_history'}.issubset(diagnostics)
          or not set(diagnostics).issubset(diagnostic_fields)):
        raise ValueError('partial reference diagnostic fields differ')
    if not isinstance(diagnostics['ancestry'], str) or not 1 <= len(diagnostics['ancestry']) <= 4096:
        raise ValueError('reference ancestry missing')
    if 'covariance_condition_estimate' in diagnostics:
        condition = diagnostics['covariance_condition_estimate']
        keys(condition, ('norm_infinity', 'reciprocal_condition_estimate', 'inverse_infinity_norm_estimate',
                         'estimation_scope', 'components'), 'reference covariance condition')
        exact(condition['estimation_scope'], 'empirical-LAPACK-pocon-not-certified-bound/v1', 'covariance estimator scope')
        for key in ('norm_infinity', 'reciprocal_condition_estimate', 'inverse_infinity_norm_estimate'):
            if decimal_string(condition[key], key) <= 0:
                raise ValueError('reference covariance condition domain differs')
        if decimal_string(condition['reciprocal_condition_estimate'], 'reciprocal condition') > 1:
            raise ValueError('reciprocal condition exceeds one')
        components = condition['components']
        if not isinstance(components, list) or len(components) != len(COMPONENT_SIZES):
            raise ValueError('covariance condition component coverage differs')
        for index, (component, rows) in enumerate(zip(components, COMPONENT_SIZES)):
            keys(component, ('component_index', 'rows', 'norm_infinity', 'reciprocal_condition_estimate',
                             'inverse_infinity_norm_estimate'), 'component covariance condition')
            exact(component['component_index'], index, 'ordered covariance component')
            exact(component['rows'], rows, 'covariance component row count')
            for key in ('norm_infinity', 'reciprocal_condition_estimate', 'inverse_infinity_norm_estimate'):
                if decimal_string(component[key], key) <= 0:
                    raise ValueError('component covariance condition domain differs')
            if decimal_string(component['reciprocal_condition_estimate'], 'component reciprocal condition') > 1:
                raise ValueError('component reciprocal condition exceeds one')
    for name in ('original_input_LAPACK_SVD', 'original_input_LAPACK_pivoted_QR'):
        if name not in diagnostics:
            continue
        peer = diagnostics[name]
        keys(peer, ('mean', 'variance', 'q_min', 'logdet_H', 'rank', 'condition2'), 'original-input peer')
        reference_vector(peer['mean'], name)
        reference_vector(peer['variance'], name, positive=True)
        exact(peer['rank'], 46, 'independent full rank')
        if decimal_string(peer['q_min'], name) < 0 or decimal_string(peer['condition2'], name) < 1:
            raise ValueError('original-input peer diagnostic domain differs')
        decimal_string(peer['logdet_H'], name)
    residuals = diagnostics['wide_residual_history']
    if not isinstance(residuals, list) or len(residuals) > 4 or (complete and len(residuals) != 4):
        raise ValueError('original-C residual refinement missing')
    for steps, residual in enumerate(residuals):
        reference_residual(residual, steps)
    postcast = diagnostics['postcast_residual_history']
    if not isinstance(postcast, list) or len(postcast) > 8 or (complete and len(postcast) != 8):
        raise ValueError('original-C postcast residual refinement missing')
    for index, residual in enumerate(postcast):
        reference_residual(residual, index % 4, length=1, digits=80 if index < 4 else 120)
    exact(work['wide_residual_evaluations'], len(residuals), 'earned basis residual evaluations')
    exact(work['postcast_wide_residual_evaluations'], len(postcast), 'earned postcast residual evaluations')
    if 'logdet_C_peer' in diagnostics:
        determinant = diagnostics['logdet_C_peer']
        keys(determinant, ('cholesky', 'eigenvalue', 'absolute_difference', 'component_count'), 'source determinant peer')
        for key in ('cholesky', 'eigenvalue', 'absolute_difference'):
            decimal_string(determinant[key], key)
        exact(determinant['component_count'], 13, 'source determinant components')
    for name, digits, corrections in (('coarse', 80, 1), ('fine', 120, 3)):
        result = output[name]
        if result is None and not complete:
            continue
        keys(result, ['digits', 'correction_steps'] + REFERENCE_VECTORS + REFERENCE_SCALARS + REFERENCE_INTERVALS, name)
        exact(result['digits'], digits, name + '/digits')
        exact(result['correction_steps'], corrections, name + '/corrections')
        for key in REFERENCE_VECTORS:
            reference_vector(result[key], name + '/' + key, positive=key == 'variance')
        for key in REFERENCE_SCALARS:
            value = decimal_string(result[key], name + '/' + key)
            if key in ('q_min', 'q_postcast', 'q_postcast_unrounded_residual', 'q_postcast_residual_rounding_allowance',
                       'excluded_mass_upper', 'reduced_forward_sensitivity') and value < 0:
                raise ValueError('negative reference quadratic/mass')
        for key in REFERENCE_INTERVALS:
            interval(result[key], name + '/' + key, highprecision=True)
        a, b = interval(result['box_probability'], name, highprecision=True)
        if not 0 < a <= b <= 1:
            raise ValueError('reference box probability domain differs')
        excluded = decimal_string(result['excluded_mass_upper'], name + '/excluded mass')
        if not 0 < excluded < Decimal(REFERENCE_SETTINGS['maximum_excluded_mass']):
            raise ValueError('reference excluded mass domain differs')
        if interval(result['log_box_probability'], name, highprecision=True)[1] > 0:
            raise ValueError('reference log probability exceeds zero')
        a, b = interval(result['median_original46'], name, highprecision=True)
        low, high = support(request)
        if a < finite(low[45], 'lower median support') or b > finite(high[45], 'upper median support'):
            raise ValueError('reference median outside literal support')
    refinement = output['refinement']
    if refinement is not None or complete:
        keys(refinement, ['estimation_scope', 'mean_absolute_errors', 'variance_absolute_errors'] +
         [key + '_absolute_error' for key in REFINEMENT_SCALARS], 'reference refinement')
        exact(refinement['estimation_scope'], 'empirical-original-input-calibration-not-certificate/v1', 'refinement scope')
        for key in ('mean_absolute_errors', 'variance_absolute_errors'):
            if any(v < 0 for v in reference_vector(refinement[key], key)):
                raise ValueError('negative reference error estimate')
        for key in REFINEMENT_SCALARS:
            if decimal_string(refinement[key + '_absolute_error'], key) < 0:
                raise ValueError('negative reference error estimate')
    if require_complete and not complete:
        raise ValueError('reference output is partial/refused; no numerical acceptance')
    return complete


def compare(native, reference, request, identity):
    check_reference(reference, request, identity)
    checks = []
    def gate(label, error, budget, **evidence):
        passed = error <= budget
        checks.append({'label': label, 'error_exact': str(error), 'budget_exact': str(budget), 'passed': passed, **evidence})
        return passed
    def discrepancy(values):
        return max(values, default=Decimal(0))
    def endpoint_error(a, b):
        return max(abs(a[0] - b[0]), abs(a[1] - b[1]))
    with localcontext() as context:
        context.prec = 180
        fine, coarse, aggregate = reference['fine'], reference['coarse'], reference['refinement']
        allocations = request['comparison_allocations']
        fraction = finite(allocations['reference_refinement_fraction'], 'refinement fraction')
        peers = [reference['diagnostics'][name] for name in ('original_input_LAPACK_SVD', 'original_input_LAPACK_pivoted_QR')]
        for key, allocation, aggregate_key in (('mean', 'coefficient', 'mean_absolute_errors'), ('variance', 'variance', 'variance_absolute_errors')):
            for index, original in enumerate(ACTIVE):
                value = decimal_string(fine[key][index], key)
                budget = finite(allocations[allocation + '_absolute'], key) + finite(allocations[allocation + '_relative'], key) * abs(value)
                measured = discrepancy([abs(value - decimal_string(peer[key][index], key)) for peer in peers + [coarse]])
                reported = decimal_string(aggregate[aggregate_key][index], key)
                gate(f'reference/{key}/original{original}/aggregate-covers-measured', measured, reported)
                gate(f'reference/{key}/original{original}/refinement', reported, fraction * budget)
                difference = abs(finite(native['completion'][key][index], key) - value)
                gate(f'native/{key}/original{original}', difference + reported, budget,
                     discrepancy_exact=str(difference), reference_error_exact=str(reported),
                     reference_exact=str(value), native_exact=str(finite(native['completion'][key][index], key)))
        for key, native_key in (('q_min', 'minimum_quadratic'), ('q_postcast', 'reported_postcast_profile_quadratic')):
            value = decimal_string(fine[key], key)
            budget = finite(allocations['quadratic_absolute'], key) + finite(allocations['quadratic_relative'], key) * abs(value)
            measured = abs(value - decimal_string(coarse[key], key))
            if key == 'q_min':
                measured = max(measured, *(abs(value - decimal_string(peer['q_min'], key)) for peer in peers))
            reported = decimal_string(aggregate[key + '_absolute_error'], key)
            if key == 'q_postcast':
                for name, result in (('coarse', coarse), ('fine', fine)):
                    rounding = decimal_string(result['q_postcast_residual_rounding_allowance'], 'postcast rounding allowance')
                    unrounded = decimal_string(result['q_postcast_unrounded_residual'], 'unrounded postcast residual')
                    gate('reference/' + name + '/postcast-rounding-allowance', abs(decimal_string(result[key], key) - unrounded), rounding)
                gate('reference/q_postcast/aggregate-covers-rounding',
                     decimal_string(fine['q_postcast_residual_rounding_allowance'], 'postcast rounding allowance'), reported)
            gate('reference/' + key + '/aggregate-covers-measured', measured, reported)
            gate('reference/' + key + '/refinement', reported, fraction * budget)
            difference = abs(finite(native['completion'][native_key], key) - value)
            gate('native/' + key, difference + reported, budget,
                 discrepancy_exact=str(difference), reference_error_exact=str(reported))
        log_budget = finite(allocations['log_relative_box_integral_absolute'], 'log allocation')
        for key, native_key in (('logdet_H', 'log_design_precision_determinant'), ('logdet_C', 'source_covariance_log_determinant')):
            value = decimal_string(fine[key], key)
            measured = abs(value - decimal_string(coarse[key], key))
            if key == 'logdet_H':
                measured = max(measured, *(abs(value - decimal_string(peer['logdet_H'], key)) for peer in peers))
            else:
                det = reference['diagnostics']['logdet_C_peer']
                peer_difference = abs(decimal_string(det['cholesky'], key) - decimal_string(det['eigenvalue'], key))
                measured = max(measured, peer_difference)
                gate('reference/logdet_C/reported-peer-difference', peer_difference, decimal_string(det['absolute_difference'], key))
            reported = decimal_string(aggregate[key + '_absolute_error'], key)
            gate('reference/' + key + '/aggregate-covers-measured', measured, reported)
            gate('reference/' + key + '/refinement', reported, fraction * log_budget)
            difference = abs(finite(native['completion'][native_key], key) - value)
            gate('native/' + key, difference + reported, log_budget,
                 discrepancy_exact=str(difference), reference_error_exact=str(reported))
        for key in ('log_prior_volume', 'log_Z_relative_box', 'log_prior_normalized_relative_evidence', 'log_observation_normalized_evidence', 'median_original46'):
            result = interval(fine[key], key, highprecision=True)
            measured = max(endpoint_error(result, interval(coarse[key], key, highprecision=True)), (result[1] - result[0]) / 2)
            reported = decimal_string(aggregate[key + '_absolute_error'], key)
            budget = finite(allocations['median_reference_refinement_absolute'], key) if key == 'median_original46' else fraction * log_budget
            gate('reference/' + key + '/aggregate-covers-measured', measured, reported)
            gate('reference/' + key + '/refinement', reported, budget)
            if key == 'log_prior_volume':
                actual = interval(native['box_diagnostics']['log_prior_volume'], key)
            elif key == 'median_original46':
                actual = interval(native['median']['quantile'], key)
                width = finite(allocations['requested_original46_median_width'], key)
                gate('reference/median_original46/width', result[1] - result[0], width)
            else:
                actual = interval(native['normalizations']['log_relative_box_integral' if key == 'log_Z_relative_box' else key], key)
            native_budget = finite(allocations['requested_original46_median_width'], key) if key == 'median_original46' else log_budget
            difference = endpoint_error(actual, result)
            gate('native/' + key, difference + reported, native_budget,
                 discrepancy_exact=str(difference), reference_error_exact=str(reported))
        gate('reference/original-C/reduced-forward-sensitivity',
             decimal_string(fine['reduced_forward_sensitivity'], 'reduced sensitivity'),
             Decimal(REFERENCE_SETTINGS['maximum_covariance_forward_sensitivity']))
        if not all(check['passed'] for check in checks):
            raise ComparisonFailure('unchanged released-box numerical comparison/refinement allocation exceeded', checks)
    return checks


def child(*arguments, **keywords):
    previous = os.environ.get('BLIS_NUM_THREADS')
    os.environ['BLIS_NUM_THREADS'] = '1'
    try:
        return receipt.child(*arguments, **keywords)
    finally:
        if previous is None:
            os.environ.pop('BLIS_NUM_THREADS', None)
        else:
            os.environ['BLIS_NUM_THREADS'] = previous


def verify_runtime(value, executable):
    keys(value, ('python', 'packages', 'thread_environment', 'inventory', 'inventory_sha256',
                 'numeric_configuration'), 'reference runtime')
    python = value['python']
    keys(python, ('version', 'implementation', 'executable'), 'reference Python')
    if not isinstance(python['version'], str) or not python['version'] or python['implementation'] != 'CPython':
        raise ValueError('reference Python version/implementation differs')
    path = Path(executable).absolute()
    exact(python['executable'], {'path': str(path), 'resolved_path': str(path.resolve()),
                                 'bytes': path.stat().st_size, 'sha256': sha256(path)}, 'reference executable')
    packages = value['packages']
    keys(packages, ('numpy', 'scipy', 'decimal'), 'reference packages')
    required = set()
    for name in ('numpy', 'scipy'):
        item = packages[name]
        keys(item, ('version', 'module_root'), name)
        if not isinstance(item['version'], str) or not re.fullmatch(r'\d+\.\d+\.\d+(?:[A-Za-z0-9.+-]*)', item['version']):
            raise ValueError('package version domain differs')
        exact(item['version'], '2.5.3' if name == 'numpy' else '1.18.1', 'frozen reference/' + name)
        root = Path(item['module_root'])
        if not root.is_absolute() or str(root.resolve()) != str(root) or not root.is_dir():
            raise ValueError('package root domain differs')
        required.update(str(p.resolve()) for p in root.rglob('*.so'))
        libraries = root.parent / (name + '.libs')
        if libraries.is_dir():
            required.update(str(p.resolve()) for p in libraries.rglob('*') if p.is_file())
    keys(packages['decimal'], ('precision_policy', 'module_path'), 'decimal package')
    exact(packages['decimal']['precision_policy'], [80, 120], 'Decimal precision')
    required.add(str(Path(packages['decimal']['module_path']).resolve()))
    exact(value['thread_environment'], dict.fromkeys(('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
                                                     'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'), '1'), 'single-thread runtime')
    inventory = value['inventory']
    if not isinstance(inventory, list) or not 1 <= len(inventory) <= 8192:
        raise ValueError('runtime inventory quota differs')
    paths, total = [], 0
    for item in inventory:
        keys(item, ('path', 'bytes', 'sha256'), 'runtime artifact')
        path = Path(item['path'])
        if not path.is_absolute() or str(path.resolve()) != str(path) or path.is_symlink() or not path.is_file():
            raise ValueError('runtime artifact path/type differs')
        size = count(item['bytes'], 268435456, 'runtime artifact bytes')
        total += size
        if total > 1073741824 or path.stat().st_size != size or sha256(path) != item['sha256']:
            raise ValueError('runtime artifact bytes/hash/quota differs')
        paths.append(str(path))
    if paths != sorted(set(paths)) or not required.issubset(paths):
        raise ValueError('complete runtime extension/library inventory differs')
    exact(value['inventory_sha256'], hashlib.sha256(json.dumps(inventory, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
          'runtime inventory identity')
    if not isinstance(value['numeric_configuration'], str) or not 1 <= len(value['numeric_configuration']) <= 65536:
        raise ValueError('numeric configuration domain differs')
    configuration = parse(value['numeric_configuration'])
    keys(configuration, ('show_config', 'threadpools', 'decimal'), 'numeric configuration')
    if not isinstance(configuration['show_config'], str) or not configuration['show_config']:
        raise ValueError('runtime backend configuration missing')
    pools = configuration['threadpools']
    if not isinstance(pools, list) or not 1 <= len(pools) <= 64:
        raise ValueError('runtime backend inventory missing or unbounded')
    pool_paths = []
    for pool in pools:
        if not isinstance(pool, dict) or not {'filepath', 'num_threads'}.issubset(pool):
            raise ValueError('runtime backend fields missing')
        exact(pool['num_threads'], 1, 'actual backend thread count')
        if not isinstance(pool['filepath'], str) or pool['filepath'] not in paths:
            raise ValueError('actual backend library missing from runtime inventory')
        pool_paths.append(pool['filepath'])
    if len(set(pool_paths)) != len(pool_paths):
        raise ValueError('duplicate actual runtime backend')
    keys(configuration['decimal'], ('version', 'libmpdec_version', 'rounding'), 'Decimal runtime')
    for key in ('version', 'libmpdec_version'):
        if not isinstance(configuration['decimal'][key], str) or not configuration['decimal'][key]:
            raise ValueError('Decimal runtime version missing')
    exact(configuration['decimal']['rounding'], 'ROUND_HALF_EVEN', 'Decimal rounding policy')


def request_bytes(blob):
    if len(blob) > 65536 or hashlib.sha256(blob).hexdigest() != REQUEST_SHA256:
        raise ValueError('exact reviewed box request SHA-256 differs')
    request = parse(blob.decode('utf-8'))
    exact(request['transport']['active_original_indices'], ACTIVE, 'request active order')
    exact(request['resources']['jobs'], 1, 'jobs')
    exact(request['resources']['threads'], 1, 'threads')
    for key, value in (('source_tree_bytes', TREE_BYTES), ('source_archive_bytes', 4194304),
                       ('source_file_bytes', 1048576), ('attempt_store_bytes', STORE_BYTES)):
        exact(request['resources'][key], value, 'controller bound/' + key)
    return request


def check_packet(folder):
    packet, _, identities = receipt.read_packet(folder)
    exact(packet['id'], 'released-box', 'packet ID')
    exact(packet['origin']['kind'], 'development_smoke', 'conditional control origin')
    exact(packet['status'], 'blocked', 'observational packet status')
    exact(packet['execution'], None, 'generic packet execution')
    return packet, identities


def support(request):
    return ([float.fromhex(row['lower_binary64_hex']) for row in request['literal_support']],
            [float.fromhex(row['upper_binary64_hex']) for row in request['literal_support']])


def nullable_enum(value, maximum, label):
    if value is not None:
        count(value, maximum, label)


def vector(value, label, positive=False):
    if not isinstance(value, list) or len(value) != 46:
        raise ValueError(label + ' requires all ordered46 values')
    result = [finite(v, label) for v in value]
    if positive and any(v <= 0 for v in result):
        raise ValueError(label + ' requires strictly positive values')
    return result


def interval(value, label, highprecision=False):
    keys(value, ('lower', 'upper'), label)
    number = decimal_string if highprecision else finite
    low, high = number(value['lower'], label), number(value['upper'], label)
    if low > high:
        raise ValueError(label + ' interval reversed')
    return low, high


def check_native(output, request, guide_sha256, require_complete=True):
    """Admit explicit partial schemas too; only complete results earn acceptance."""
    keys(output, NATIVE_FIELDS, 'native')
    exact(output['schema_version'], 1, 'native schema')
    exact(output['interface_id'], request['consumer_interface_id'], 'native interface')
    exact(output['contract_sha256'], CONTRACT_SHA256, 'native original contract')
    exact(output['request_sha256'], REQUEST_SHA256, 'native request')
    sdk = {**request['sdk_identity'], 'gaussian_box_guide_sha256': guide_sha256,
           'verification_owner': 'distinct bounded controller; echoed compile-time declarations are not self-certified'}
    exact(output['sdk_identity'], sdk, 'native SDK')
    low, high = support(request)
    target = {'observations': 3492, 'original_columns': 47, 'active_columns': 46,
              'active_original_indices': ACTIVE,
              'fixed_coordinates': [{'original_index': 44, 'value': 0,
                                     'measure': 'point mass outside active46 Lebesgue'}],
              'parameter_measure': request['transport']['parameter_measure'],
              'prior_identity': 'lstsq_results.txt@c447f0fea703fcd0fff57de5000947b5ca81286b;SHA256=' +
                                request['lineage']['immutable_prior_sha256'] + ';literal closed endpoints',
              'requested_marginal': {'original_index': 46, 'active_index': 45,
                                     'cumulative_probability': 0.5, 'unit': '5log10(H0/[1km/s/Mpc])'},
              'support_lower': low, 'support_upper': high, 'endpoints': 'closed literal binary64',
              'row_order': 'released-row-0 through released-row-3491; unchanged original FITS order'}
    exact(output['target'], target, 'native target')
    original = request['lineage']
    canonical = request['transport']['canonical_inputs']
    sources = {'source_revision': original['source_revision'],
               'packet_source_revision': original['packet_source_revision'],
               'lineage_sha256': original['packet_lineage_sha256'],
               'constrained_sha256': original['constrained_sha256'],
               **{'canonical_' + item['name'][0] + '_sha256': item['sha256'] for item in canonical},
               'hash_admission_owner': 'controller before/after source and canonical byte verification'}
    exact(output['source_identities'], sources, 'native source identities')
    exact(output['arithmetic'], {'covariance_id': 'F02/longdouble-cpu/v1', 'qr_id': 'longdouble-cpu/v1',
                               'double_mantissa_bits': 53, 'long_double_mantissa_bits': 64,
                               'long_double_max_exponent': 16384, 'round_to_nearest': True,
                               'binary64_ieee': True, 'little_endian': True}, 'native arithmetic')
    policy = {**request['numerical_policy']['design'],
              **{k: v for k, v in request['numerical_policy']['box'].items()
                 if k != 'maximum_total_cdf_nodes_scope'}}
    exact(output['policy'], policy, 'native producer policy')
    limits = request['resources']
    exact(output['resources'], {'jobs': 1, 'threads': 1,
                              **{k: limits[k] for k in ('address_limit_bytes', 'output_limit_bytes', 'native_cpu_seconds')},
                              'wall_timeout_owner': 'controller'}, 'native resources')
    stage = output['stages']
    keys(stage, STAGES, 'native stages')
    if stage['last_stage'] not in ('input_admission', 'gaussian_preparation', 'design_preparation',
                                   'box_preparation', 'box_evaluation', 'native_complete',
                                   'native_refusal_or_output_postcondition'):
        raise ValueError('native last stage differs')
    if stage['error'] is not None and (not isinstance(stage['error'], str) or not 1 <= len(stage['error']) <= 256):
        raise ValueError('native error domain differs')
    for key in ('input_complete', 'design_consumed_by_box'):
        if type(stage[key]) is not bool:
            raise ValueError('native stage flags require boolean')
    for key in ('gaussian_status', 'design_status', 'box_preparation_status'):
        nullable_enum(stage[key], 5, key)
    for key in ('gaussian_numerical_status', 'design_numerical_status'):
        nullable_enum(stage[key], 8, key)
    nullable_enum(stage['design_rank'], 3, 'design rank')
    for key in ('equilibrated_triangular_condition_inf', 'equilibrated_transpose_triangular_condition_inf'):
        if stage[key] is not None and finite(stage[key], key) < 1:
            raise ValueError('condition domain differs')
    exact(stage['design_method'], 'retained-whitened-pivoted-householder-qr/v1', 'design method')
    exact(stage['box_preparation_numerical_status'], None, 'unexposed box numerical status')
    exact(stage['box_preparation_numerical_status_scope'], 'not exposed by current API', 'box status scope')
    keys(output['payload_bounds'], ('gaussian_preparation', 'design_preparation', 'box_preparation', 'box_evaluation'), 'payload bounds')
    for key, value in output['payload_bounds'].items():
        if value is not None:
            count(value, sys.maxsize, key)
    for key, maximum in (('result_status', 5), ('result_numerical_status', 8), ('result_stage', 6)):
        nullable_enum(output[key], maximum, key)
    for key in ('output_complete', 'accepted'):
        if type(output[key]) is not bool:
            raise ValueError(key + ' requires boolean')
    exact(output['qualification'], 'native conditional completion/enclosure only; independent original-input comparison unassessed; observational qualification blocked', 'native qualification')
    flags = output['availability']
    result = output['result_stage'] is not None
    if result:
        count(output['result_status'], 5, 'assessed result status')
        count(output['result_numerical_status'], 8, 'assessed numerical status')
        exact(output['method_id'], 'retained-qr-rational-tail-box-enclosure/v1', 'box method')
        exact(output['enclosure_scope'], request['qualification']['enclosure_scope'], 'box enclosure scope')
        keys(flags, FLAGS, 'availability')
        for index, flag in enumerate(FLAGS, 1):
            exact(flags[flag], output['result_stage'] >= index, 'causal availability/' + flag)
    else:
        for key in ('availability', 'result_status', 'result_numerical_status', 'method_id', 'enclosure_scope'):
            exact(output[key], None, 'unassessed/' + key)
        flags = dict.fromkeys(FLAGS, False)
    for key, flag in (('completion', 'gaussian_completion'), ('box_diagnostics', 'endpoint_margins'),
                      ('normalizations', 'normalization_enclosures'), ('median', 'quantile_enclosure')):
        if not flags[flag]:
            exact(output[key], None, 'unavailable/' + key)
        elif not isinstance(output[key], dict):
            raise ValueError('available group missing: ' + key)
    if flags['gaussian_completion']:
        completion = output['completion']
        keys(completion, COMPLETION, 'Gaussian completion')
        vector(completion['mean'], 'unboxed mean')
        vector(completion['variance'], 'unboxed variance', positive=True)
        exact(completion['named_variance_original46'], completion['variance'][45], 'named variance')
        for key in COMPLETION[2:]:
            value = finite(completion[key], key)
            if key in ('minimum_quadratic', 'reported_postcast_profile_quadratic',
                       'maximum_variance_sensitivity', 'profile_stationarity') and value < 0:
                raise ValueError('negative completion diagnostic ' + key)
    if flags['endpoint_margins']:
        diagnostics = output['box_diagnostics']
        keys(diagnostics, DIAGNOSTICS, 'box diagnostics')
        for key in DIAGNOSTICS[:2]:
            if not isinstance(diagnostics[key], list) or len(diagnostics[key]) != 46:
                raise ValueError('all46 endpoint margins required')
            for value in diagnostics[key]:
                interval(value, key)
        if finite(diagnostics['excluded_mass_upper'], 'excluded mass') < 0:
            raise ValueError('negative excluded mass')
        interval(diagnostics['log_prior_volume'], 'prior volume')
        for key in ('box_probability', 'log_box_probability'):
            if flags['rectangle_enclosure']:
                interval(diagnostics[key], key)
            else:
                exact(diagnostics[key], None, 'unavailable/' + key)
    if flags['normalization_enclosures']:
        keys(output['normalizations'], NORMALIZATIONS, 'normalizations')
        for key in NORMALIZATIONS:
            interval(output['normalizations'][key], key)
    if flags['quantile_enclosure']:
        median = output['median']
        keys(median, ('original_index', 'active_index', 'cumulative_probability', 'quantile',
                      'lower_endpoint_cdf', 'upper_endpoint_cdf'), 'median')
        for key, value in (('original_index', 46), ('active_index', 45), ('cumulative_probability', 0.5)):
            exact(median[key], value, 'median/' + key)
        a, b = interval(median['quantile'], 'median quantile')
        if a < finite(low[45], 'lower support') or b > finite(high[45], 'upper support'):
            raise ValueError('median outside literal support')
        for key in ('lower_endpoint_cdf', 'upper_endpoint_cdf'):
            if flags['endpoint_cdf_enclosures']:
                a, b = interval(median[key], key)
                if a < 0 or b > 1:
                    raise ValueError('conditional CDF probability domain differs')
                if output['accepted'] and ((key == 'lower_endpoint_cdf' and a > Decimal('0.5')) or
                                           (key == 'upper_endpoint_cdf' and b < Decimal('0.5'))):
                    raise ValueError('endpoint CDF contradicts reported median bracket')
            else:
                exact(median[key], None, 'unavailable/' + key)
    keys(output['work'], ('cdf_node_evaluations', 'cdf_evaluations', 'bisections'), 'work')
    for key, value in output['work'].items():
        if result:
            count(value, sys.maxsize, key)
        else:
            exact(value, None, 'unassessed work/' + key)
    if output['accepted'] != output['output_complete']:
        raise ValueError('native acceptance/completion disagree')
    if output['accepted']:
        for key in ('gaussian_status', 'design_status', 'box_preparation_status',
                    'gaussian_numerical_status', 'design_numerical_status'):
            exact(stage[key], 0, 'success/' + key)
        exact(stage['design_rank'], 3, 'success rank')
        exact(stage['input_complete'], True, 'success input')
        exact(stage['design_consumed_by_box'], True, 'success ownership')
        exact(stage['error'], None, 'success error')
        exact(stage['last_stage'], 'native_complete', 'success last stage')
        exact(output['result_status'], 0, 'success result')
        exact(output['result_numerical_status'], 0, 'success numerical')
        exact(output['result_stage'], 6, 'success complete stage')
        for key, value in output['payload_bounds'].items():
            count(value, policy['maximum_payload_bytes'], 'success payload/' + key)
            if value == 0:
                raise ValueError('zero success payload bound')
        for key in ('maximum_variance_sensitivity', 'profile_stationarity'):
            if finite(completion[key], key) > finite(policy['maximum_forward_sensitivity'], key):
                raise ValueError('forward sensitivity budget exceeded')
        with localcontext() as context:
            context.prec = 128
            e = finite(diagnostics['excluded_mass_upper'], 'excluded mass')
            probability = interval(diagnostics['box_probability'], 'box probability')
            log_probability = interval(diagnostics['log_box_probability'], 'log probability')
            if not 0 < e < 1 or not 0 < probability[0] <= probability[1] <= 1 or log_probability[1] > 0:
                raise ValueError('box probability/log domain differs')
            if probability[0] > 1 - e or probability[1] != 1 or log_probability[1] != 0:
                raise ValueError('reported union-bound mass identities differ')
            if log_probability[1] - log_probability[0] > finite(policy['maximum_log_probability_width'], 'log width'):
                raise ValueError('log probability allocation exceeded')
            a, b = interval(output['median']['quantile'], 'median')
            if b - a > finite(policy['maximum_quantile_width'], 'median width'):
                raise ValueError('median allocation exceeded')
            for index, (a, b, mean, variance) in enumerate(zip(low, high, completion['mean'], completion['variance'])):
                sigma = finite(variance, 'variance').sqrt()
                expected = ((finite(mean, 'mean') - finite(a, 'support')) / sigma,
                            (finite(b, 'support') - finite(mean, 'mean')) / sigma)
                for key, wanted in zip(DIAGNOSTICS[:2], expected):
                    margin = interval(diagnostics[key][index], key)
                    if not margin[0] <= wanted <= margin[1]:
                        raise ValueError('reported completion/endpoint margin identity differs')
            qmin = finite(completion['minimum_quadratic'], 'minimum quadratic')
            qcast = finite(completion['reported_postcast_profile_quadratic'], 'postcast quadratic')
            allocations = request['comparison_allocations']
            if abs(qmin - qcast) > finite(allocations['quadratic_absolute'], 'q allocation') + \
                    finite(allocations['quadratic_relative'], 'q allocation') * abs(qmin):
                raise ValueError('minimum/postcast quadratic discrepancy exceeds frozen allocation')
            volume = interval(diagnostics['log_prior_volume'], 'prior volume')
            exact_volume = sum((finite(b, 'support') - finite(a, 'support')).ln()
                               for a, b in zip(low, high))
            if not volume[0] <= exact_volume <= volume[1]:
                raise ValueError('literal support volume identity differs')
            log2pi = (2 * PI).ln()
            unboxed = -qmin / 2 + 23 * log2pi - finite(completion['log_design_precision_determinant'], 'design determinant') / 2
            relative = interval(output['normalizations'][NORMALIZATIONS[0]], NORMALIZATIONS[0])
            prior = interval(output['normalizations'][NORMALIZATIONS[1]], NORMALIZATIONS[1])
            observation = interval(output['normalizations'][NORMALIZATIONS[2]], NORMALIZATIONS[2])
            expected_relative = (unboxed + log_probability[0], unboxed + log_probability[1])
            expected_prior = (relative[0] - volume[1], relative[1] - volume[0])
            observation_offset = (finite(completion['source_covariance_log_determinant'], 'source determinant') + 3492 * log2pi) / 2
            expected_observation = (prior[0] - observation_offset, prior[1] - observation_offset)
            for actual, expected, label in ((relative, expected_relative, 'relative integral'),
                                             (prior, expected_prior, 'prior evidence'),
                                             (observation, expected_observation, 'observation evidence')):
                if actual[1] < expected[0] or expected[1] < actual[0]:
                    raise ValueError('normalization identity differs: ' + label)
        for key, maximum in (('cdf_node_evaluations', policy['maximum_total_cdf_nodes']),
                             ('cdf_evaluations', policy['maximum_cdf_evaluations']),
                             ('bisections', 2 * policy['maximum_bisections_per_inverse'])):
            count(output['work'][key], maximum, 'success work/' + key)
    elif require_complete:
        raise ValueError('native output is partial/refused; no numerical acceptance')
    return output['accepted']


def file_identity(path, maximum, expected_sha256=None, expected_bytes=None):
    blob = receipt.read_log(path, maximum)
    identity = {'path': str(path), 'bytes': len(blob), 'sha256': hashlib.sha256(blob).hexdigest()}
    if expected_sha256 is not None and identity['sha256'] != expected_sha256:
        raise ValueError('changed exact input: ' + str(path))
    if expected_bytes is not None and len(blob) != expected_bytes:
        raise ValueError('changed input byte count: ' + str(path))
    return blob, identity


def raw_sources(directory, manifest):
    sources, identities = {}, []
    for item in manifest['sources']:
        path = directory / item['name']
        if path.is_symlink() or path.resolve().parent != directory.resolve():
            raise ValueError('raw source path/type differs')
        blob, identity = file_identity(path, item['bytes'], item['sha256'], item['bytes'])
        sources[item['name']] = blob
        identities.append({**item, **identity})
    return sources, identities


def canonical_identity(directory, request):
    result = []
    for item in request['transport']['canonical_inputs']:
        _, identity = file_identity(directory / item['name'], item['bytes'], item['sha256'], item['bytes'])
        result.append({**item, **identity})
    return result


def engine_identity(source, artifacts, request):
    sdk = request['sdk_identity']
    head = receipt.git(source, 'rev-parse', 'HEAD')
    status = receipt.git(source, 'status', '--porcelain', '--untracked-files=all')
    if head != sdk['engine_revision'] or status:
        raise ValueError('clean pinned GaussianBox engine source required')
    paths = {'manifest': artifacts / 'build/build-manifest-release.json',
             'archive': artifacts / 'build/native-release/libirred_core.a',
             'cli': artifacts / 'target/release/irred'}
    hashes = {key: sha256(path) for key, path in paths.items()}
    for key, digest in hashes.items():
        exact(digest, sdk[key + '_sha256'], 'actual engine/' + key)
    build = load(paths['manifest'])
    for key, value in (('git_head', head), ('git_status', ''), ('build_id', sdk['build_id']),
                       ('profile', 'release'), ('backend', 'portable_cpu')):
        exact(build[key], value, 'build/' + key)
    content = {k: v for k, v in build.items() if k not in ('build_id', 'git_head', 'git_status')}
    exact(hashlib.sha256(json.dumps(content, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
          build['build_id'], 'reconstructed build identity')
    patterns = ('src/**/*.rs', 'tests/**/*.rs', 'cpp/**/*.cpp', 'cpp/**/*.h', 'cpp/**/*.hpp',
                'cpp/**/*.inc', 'cpp/**/*.cmake', 'schema/*.json', 'tools/*.py')
    inventory = {str(path.relative_to(source)) for pattern in patterns for path in source.glob(pattern)}
    inventory.update(('Cargo.toml', 'Cargo.lock', 'build.rs', 'cpp/CMakeLists.txt'))
    if inventory != set(build['sources']):
        raise ValueError('complete GaussianBox engine inventory differs')
    exact(len(inventory), sdk['source_inventory_count'], 'complete SDK source count')
    for name, digest in build['sources'].items():
        path = source / name
        if path.is_symlink() or not path.is_file() or not path.resolve().is_relative_to(source.resolve()):
            raise ValueError('engine source type/path differs: ' + name)
        exact(sha256(path), digest, 'engine source/' + name)
    tools = {'compiler': Path('/usr/bin/c++').resolve(),
             'standard_library': Path(build['standard_library']).resolve()}
    expected = {'compiler': build['compiler_executable_digest'],
                'standard_library': build['standard_library_digest']}
    for name, digest in build['tool_executable_digests'].items():
        path = shutil.which(name)
        if path is None:
            raise ValueError('missing build tool: ' + name)
        tools[name] = Path(path).resolve()
        expected[name] = digest
    tools = {name: {'path': str(path), 'bytes': path.stat().st_size, 'sha256': sha256(path)}
             for name, path in tools.items()}
    for name in tools:
        exact(tools[name]['sha256'], expected[name], 'actual build tool/' + name)
    result = {**hashes, 'source_path': str(source), 'artifact_path': str(artifacts), 'head': head,
              'status': status, 'build_id': build['build_id'], 'sources': build['sources'], 'tools': tools}
    exact(result['tools']['compiler']['sha256'], sdk['compiler_executable_sha256'], 'compiler identity')
    exact(result['tools']['standard_library']['sha256'], sdk['standard_library_sha256'], 'standard library identity')
    result['gaussian_box_header_sha256'] = sha256(source / 'cpp/include/irred/gaussian_box.hpp')
    exact(result['gaussian_box_header_sha256'], sdk['gaussian_box_header_sha256'], 'GaussianBox header')
    result['gaussian_box_guide_sha256'] = sha256(source / 'docs/gaussian-box.md')
    exact(result['gaussian_box_guide_sha256'], hashlib.sha256(subprocess.check_output(
          ['git', '-C', str(source), 'show', sdk['engine_revision'] + ':docs/gaussian-box.md'], timeout=30)).hexdigest(),
          'committed GaussianBox guide')
    return result


def snapshot_source(store, identity, record):
    """New packet's reviewed2MiB bound; historical thermal source limit is intact."""
    if identity['status']:
        raise ValueError('final attempts require clean committed Reproducible source')
    child(['git', '-C', str(ROOT), 'archive', '--format=tar', '--output=' + str(store / 'source.tar'), identity['head']],
          store, 'source', {'source_seconds': 30, 'memory_bytes': 1073741824,
                           'output_bytes': 4194304}, record=record)
    destination = store / 'source'
    destination.mkdir()
    seen, total = set(), 0
    with tarfile.open(store / 'source.tar', 'r:') as archive:
        for member in archive:
            path = within(destination, member.name)
            if member.isdir():
                path.mkdir(parents=True, exist_ok=True)
                continue
            total += member.size
            if not member.isfile() or member.size > 1048576 or total > TREE_BYTES or member.name in seen:
                raise ValueError('committed snapshot type/size/duplicate refused')
            blob = archive.extractfile(member).read(member.size + 1)
            if len(blob) != member.size or hashlib.sha256(blob).hexdigest() != identity['files'].get(member.name):
                raise ValueError('committed snapshot bytes differ')
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(blob)
            seen.add(member.name)
    if seen != set(identity['files']):
        raise ValueError('committed snapshot inventory differs')
    for relative in SHARED_PATHS + (str(Path(__file__).resolve().relative_to(ROOT)),):
        if relative not in seen:
            raise ValueError('controller/helper source is not committed')
    record['source_archive_sha256'] = sha256(store / 'source.tar')
    return destination


def terminal_outputs(store, record, errors):
    drift, stdout = receipt.verify_outputs(store, record)
    errors.extend(drift)
    receipt.ingest_outputs(store, record, stdout)
    record['output_sha256'] = {}
    total = 0
    for path in store.rglob('*'):
        if path.is_symlink():
            errors.append('attempt symlink refused: ' + str(path.relative_to(store)))
        elif path.is_file():
            name = str(path.relative_to(store))
            try:
                total += path.stat().st_size
                record['output_sha256'][name] = sha256(path)
            except Exception as error:
                errors.append(name + ' terminal output identity: ' + str(error))
    if total > STORE_BYTES:
        errors.append('attempt storage quota exceeded')
    record['output_bytes'] = total
    for label, operation in record.get('subprocesses', {}).items():
        for suffix in ('out', 'err'):
            name = label + '.' + suffix
            if operation.get(suffix + '_sha256') is not None and record['output_sha256'].get(name) != operation[suffix + '_sha256']:
                errors.append(name + ' raw log drift: terminal manifest differs')
    if errors:
        record['status'] = 'failed'
        record.setdefault('error', 'terminal identity failure')
        record['gates']['numerical'] = 'not accepted: identity failure'


def execute(args):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,95}', args.name):
        raise ValueError('attempt name must be one bounded path component')
    store = within(ROOT, 'results/released-box/' + args.name)
    store.mkdir(parents=True, exist_ok=False)
    record = {'schema_version': 1, 'status': 'started', 'request_sha256': REQUEST_SHA256,
              'started_utc': datetime.now(timezone.utc).isoformat(),
              'gates': {'execution': 'unassessed', 'numerical': 'unassessed',
                        'inference': 'conditional released active46 box median only',
                        'interpretation': 'observational qualification blocked: unchanged source gaps'},
              'paths': {'engine_source': str(args.engine_source), 'engine_artifacts': str(args.engine_artifacts),
                        'original_sources': str(args.sources)},
              'controller_python': {'version': sys.version, 'path': str(Path(sys.executable).absolute()),
                                    'resolved_path': str(Path(sys.executable).resolve()),
                                    'sha256': sha256(Path(sys.executable))},
              'controller_bounds': {'source_tree_bytes': TREE_BYTES, 'source_archive_bytes': 4194304,
                                    'source_file_bytes': 1048576, 'attempt_store_bytes': STORE_BYTES}}
    started = time.monotonic()
    usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    source_before = engine_before = raw_before = canonical_before = runtime = identity = None
    snapshot_hashes, packet_hashes, immutable = {}, {}, {}
    request = None
    try:
        source_before = receipt.source_identity()
        record['source_before'] = source_before
        if any(source_before['files'].get(name) != digest for name, digest in LOADED_SOURCE_HASHES.items()):
            raise ValueError('loaded controller/helper bytes differ from committed source identity')
        record['loaded_source_hashes'] = LOADED_SOURCE_HASHES
        record['packet'], record['packet_identities'] = check_packet(FOLDER)
        request = request_bytes(receipt.read_log(FOLDER / 'request.json', 65536))
        record['request'] = request
        packet_hashes = {p.name: sha256(p) for p in FOLDER.iterdir() if p.is_file() and not p.is_symlink()}
        if len(packet_hashes) > 8:
            raise ValueError('new packet exceeds eight source files')
        record['packet_hashes'] = packet_hashes
        snapshot = snapshot_source(store, source_before, record)
        packet = snapshot / FOLDER.relative_to(ROOT)
        for name, digest in packet_hashes.items():
            blob, _ = file_identity(packet / name, 1048576, digest)
            (store / name).write_bytes(blob)
        snapshot_hashes = {str(p.relative_to(store)): sha256(p) for p in store.rglob('*') if p.is_file()}
        record['snapshot_hashes'] = snapshot_hashes
        original = snapshot / 'experiments/released-ladder'
        for name, digest in (('lineage.json', request['lineage']['packet_lineage_sha256']),
                             ('constrained.json', request['lineage']['constrained_sha256'])):
            file_identity(original / name, 65536, digest)
        contract_blob, _ = file_identity(store / 'contract.snapshot.json', 65536, CONTRACT_SHA256)
        contract = parse(contract_blob.decode())
        literal = [{k: row[k] for k in ('original_index', 'lower_binary64_hex', 'upper_binary64_hex')}
                   for row in contract['all_original_prior_rows'] if row['original_index'] != 44]
        exact(request['literal_support'], literal, 'literal original contract support')
        exact(sha256(store / 'bounds.hpp'), request['bounds_header_sha256'], 'bounds transport source')
        manifest = load(original / 'lineage.json')
        lineage.validate_lineage(manifest)
        lineage.validate_constrained(load(original / 'constrained.json'))
        sources, raw_before = raw_sources(args.sources, manifest)
        record['raw_inputs_before'] = raw_before
        record['source_audit'] = lineage.source_audit(sources, manifest, store)
        record['constrained_source_audit'] = lineage.constrained_source_audit(store, sources)
        canonical_before = canonical_identity(store, request)
        record['canonical_inputs_before'] = canonical_before
        del sources
        engine_before = engine_identity(args.engine_source, args.engine_artifacts, request)
        record['engine_before'] = engine_before
        resources = request['resources']
        if len(str(store).encode()) > resources['input_directory_maximum_bytes']:
            raise ValueError('canonical directory path byte quota exceeded')
        limits = {'memory_bytes': resources['address_limit_bytes'], 'output_bytes': resources['output_limit_bytes'],
                  'compile_seconds': resources['compile_seconds'], 'native_seconds': resources['native_wall_seconds'],
                  'reference_seconds': resources['reference_seconds'], 'runtime_seconds': 30, 'discovery_seconds': 30}
        record['subprocess_limits'] = limits
        discovery = parse(child([str(args.engine_artifacts / 'target/release/irred'), 'describe', '--json'],
                                store, 'discovery', limits, record=record))
        record['discovery'] = discovery
        if (discovery.get('product') != 'Irreducible' or discovery['build']['build_id'] != request['sdk_identity']['build_id'] or
                discovery['build']['git_head'] != request['sdk_identity']['engine_revision'] or discovery['build']['git_status'] != ''):
            raise ValueError('compiled discovery identity differs')
        record['native_route'] = request['native_route']
        sdk = request['sdk_identity']
        definitions = dict(zip(('ENGINE_REVISION', 'BUILD_ID', 'MANIFEST_SHA256', 'ARCHIVE_SHA256', 'CLI_SHA256', 'HEADER_SHA256'),
                               (sdk['engine_revision'], sdk['build_id'], sdk['manifest_sha256'], sdk['archive_sha256'],
                                sdk['cli_sha256'], sdk['gaussian_box_header_sha256'])))
        command = ['/usr/bin/c++', '-std=c++20', '-O2', '-Wall', '-Wextra', '-Wpedantic', '-fno-fast-math', '-ffp-contract=off',
                   *['-DRELEASED_BOX_' + key + '="' + value + '"' for key, value in definitions.items()],
                   str(store / 'consumer.cpp'), '-I', str(args.engine_source / 'cpp/include'),
                   str(args.engine_artifacts / 'build/native-release/libirred_core.a'), '-o', str(store / 'consumer')]
        child(command, store, 'compile', limits, record=record)
        immutable['consumer'] = sha256(store / 'consumer')
        record['consumer_sha256'] = immutable['consumer']
        native = parse(child([str(store / 'consumer'), str(store)], store, 'native', limits, record=record))
        record['native'] = native
        check_native(native, request, engine_before['gaussian_box_guide_sha256'])
        record['gates']['execution'] = 'passed complete native conditional box result'
        reference_input = {'schema_version': 1, 'request_sha256': REQUEST_SHA256,
                           'canonical_directory': str(store), 'original_source_directory': str(args.sources),
                           'source_root': str(snapshot)}
        (store / 'reference-input.json').write_text(json.dumps(reference_input, allow_nan=False) + '\n')
        immutable['reference-input.json'] = sha256(store / 'reference-input.json')
        runtime = parse(child([str(args.reference_python), str(store / 'reference.py'), '--runtime-fingerprint'],
                              store, 'runtime', limits, record=record))
        record['reference_runtime'] = runtime
        verify_runtime(runtime, args.reference_python)
        identity = {'runtime': runtime, 'reference_input_sha256': immutable['reference-input.json'],
                    'reference_script_sha256': sha256(store / 'reference.py'), 'request_sha256': REQUEST_SHA256}
        reference = parse(child([str(args.reference_python), str(store / 'reference.py'), str(store / 'reference-input.json')],
                                store, 'reference', limits, record=record))
        record['reference'] = reference
        record['comparisons'] = compare(native, reference, request, identity)
        record['gates']['numerical'] = 'passed empirical original-input completion, box integral and median/refinement comparisons; not a certificate'
        record['status'] = 'completed'
    except Exception as error:
        record['status'] = 'failed'
        record['error'] = str(error)
        if isinstance(error, ComparisonFailure):
            record['comparisons'] = error.checks
        if 'native' in record and record['gates']['execution'] == 'unassessed':
            record['gates']['execution'] = 'failed native subprocess/admission'
    finally:
        errors, stdout = receipt.verify_outputs(store, record)
        receipt.ingest_outputs(store, record, stdout)
        if 'native' in record and request is not None and engine_before is not None and record['status'] == 'failed':
            try:
                record['native_output_complete'] = check_native(
                    record['native'], request, engine_before['gaussian_box_guide_sha256'], require_complete=False)
                record['native_output_schema_admitted'] = True
            except Exception as error:
                record['native_partial_admission_error'] = str(error)
        if 'reference' in record and identity is not None and record['status'] == 'failed':
            try:
                record['reference_output_complete'] = check_reference(
                    record['reference'], request, identity, require_complete=False)
                record['reference_output_schema_admitted'] = True
            except Exception as error:
                record['reference_partial_admission_error'] = str(error)
        checks = []
        if source_before is not None:
            checks.append(('source_after', receipt.source_identity, source_before))
        if engine_before is not None:
            checks.append(('engine_after', lambda: engine_identity(args.engine_source, args.engine_artifacts, request), engine_before))
        if raw_before is not None:
            checks.append(('raw_inputs_after', lambda: raw_sources(args.sources, manifest)[1], raw_before))
        if canonical_before is not None:
            checks.append(('canonical_inputs_after', lambda: canonical_identity(store, request), canonical_before))
        if packet_hashes:
            checks.append(('packet_hashes_after', lambda: {p.name: sha256(p) for p in FOLDER.iterdir()
                                                          if p.is_file() and not p.is_symlink()}, packet_hashes))
        for key, read, expected in checks:
            try:
                record[key] = read()
                if record[key] != expected:
                    errors.append(key + ' identity changed')
            except Exception as error:
                errors.append(key + ' final verification: ' + str(error))
        if runtime is not None:
            try:
                verify_runtime(runtime, args.reference_python)
            except Exception as error:
                errors.append('runtime final verification: ' + str(error))
        for name, digest in {**snapshot_hashes, **immutable}.items():
            try:
                if sha256(store / name) != digest:
                    errors.append(name + ' immutable snapshot changed')
            except Exception as error:
                errors.append(name + ' immutable snapshot verification: ' + str(error))
        terminal_outputs(store, record, errors)
        record['integrity_errors'] = errors
        record['elapsed_seconds'] = time.monotonic() - started
        record['finished_utc'] = datetime.now(timezone.utc).isoformat()
        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        record['children_resources'] = {'cpu_user_seconds': usage.ru_utime - usage_before.ru_utime,
                                        'cpu_system_seconds': usage.ru_stime - usage_before.ru_stime,
                                        'process_lifetime_maximum_rss_kib': usage.ru_maxrss, 'jobs': 1, 'threads': 1}
        (store / 'record.json').write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')
        for path in store.rglob('*'):
            if path.is_file() and not path.is_symlink():
                path.chmod(0o555 if path.name == 'consumer' else 0o444)
    print(json.dumps({'status': record['status'], 'record': str(store / 'record.json'),
                      'gates': record['gates'], 'error': record.get('error')}))
    return 0 if record['status'] == 'completed' else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('engine-source', 'sources', 'reference-python'):
        parser.add_argument('--' + name, required=True, type=Path)
    parser.add_argument('--engine-artifacts', type=Path)
    parser.add_argument('--name', required=True)
    args = parser.parse_args()
    args.engine_source = args.engine_source.resolve()
    args.engine_artifacts = (args.engine_artifacts or args.engine_source).resolve()
    args.sources = args.sources.resolve()
    args.reference_python = args.reference_python.absolute()
    sys.exit(execute(args))
