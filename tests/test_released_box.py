"""Synthetic admission/receipt controls; these do not qualify released science."""
import copy
from decimal import Decimal, localcontext
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

FOLDER = Path(__file__).resolve().parent
ROOT = next(path for path in FOLDER.parents if (path / 'experiments/released-ladder/controller.py').is_file())
PUBLIC = ROOT / 'experiments/released-box/controller.py'
CONTROLLER = PUBLIC if PUBLIC.is_file() else FOLDER / 'controller.py'
spec = importlib.util.spec_from_file_location('box_controller_test', CONTROLLER)
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)
REQUEST = c.request_bytes((CONTROLLER.parent / 'request.json').read_bytes())
GUIDE = '7db2a06a5ef95c07729d0dbec4acd7d06e3b0eeb0678137e3d5c06e1d1a35571'


def enclosure(value, radius=0):
    value = float(value)
    return {'lower': math.nextafter(value - radius, -math.inf),
            'upper': math.nextafter(value + radius, math.inf)}


def native_fixture():
    """A schema fixture with internally consistent scalar identities, no data fit."""
    r = REQUEST
    low, high = c.support(r)
    means = [(a + b) / 2 for a, b in zip(low, high)]
    variances = [0.00001] * 46
    with localcontext() as context:
        context.prec = 128
        volume = sum((c.finite(b, 'support') - c.finite(a, 'support')).ln() for a, b in zip(low, high))
        log2pi = (2 * c.PI).ln()
        relative = -1 + 23 * log2pi
        prior = relative - volume
        observation = prior - 1746 * log2pi
        margins = [[enclosure((c.finite(mean, 'mean') - c.finite(a, 'support')) / c.finite(var, 'variance').sqrt()),
                    enclosure((c.finite(b, 'support') - c.finite(mean, 'mean')) / c.finite(var, 'variance').sqrt())]
                   for a, b, mean, var in zip(low, high, means, variances)]
    o = {'schema_version': 1, 'interface_id': r['consumer_interface_id'], 'contract_sha256': c.CONTRACT_SHA256,
         'request_sha256': c.REQUEST_SHA256,
         'sdk_identity': {**r['sdk_identity'], 'gaussian_box_guide_sha256': GUIDE,
                          'verification_owner': 'distinct bounded controller; echoed compile-time declarations are not self-certified'},
         'target': {'observations': 3492, 'original_columns': 47, 'active_columns': 46,
                    'active_original_indices': c.ACTIVE,
                    'fixed_coordinates': [{'original_index': 44, 'value': 0, 'measure': 'point mass outside active46 Lebesgue'}],
                    'parameter_measure': r['transport']['parameter_measure'],
                    'prior_identity': 'lstsq_results.txt@c447f0fea703fcd0fff57de5000947b5ca81286b;SHA256=' + r['lineage']['immutable_prior_sha256'] + ';literal closed endpoints',
                    'requested_marginal': {'original_index': 46, 'active_index': 45, 'cumulative_probability': 0.5, 'unit': '5log10(H0/[1km/s/Mpc])'},
                    'support_lower': low, 'support_upper': high, 'endpoints': 'closed literal binary64',
                    'row_order': 'released-row-0 through released-row-3491; unchanged original FITS order'},
         'source_identities': {'source_revision': r['lineage']['source_revision'], 'packet_source_revision': r['lineage']['packet_source_revision'],
                               'lineage_sha256': r['lineage']['packet_lineage_sha256'], 'constrained_sha256': r['lineage']['constrained_sha256'],
                               **{'canonical_' + item['name'][0] + '_sha256': item['sha256'] for item in r['transport']['canonical_inputs']},
                               'hash_admission_owner': 'controller before/after source and canonical byte verification'},
         'arithmetic': {'covariance_id': 'F02/longdouble-cpu/v1', 'qr_id': 'longdouble-cpu/v1', 'double_mantissa_bits': 53,
                        'long_double_mantissa_bits': 64, 'long_double_max_exponent': 16384, 'round_to_nearest': True,
                        'binary64_ieee': True, 'little_endian': True},
         'policy': {**r['numerical_policy']['design'], **{k: v for k, v in r['numerical_policy']['box'].items() if k != 'maximum_total_cdf_nodes_scope'}},
         'resources': {'jobs': 1, 'threads': 1, **{k: r['resources'][k] for k in ('address_limit_bytes', 'output_limit_bytes', 'native_cpu_seconds')},
                       'wall_timeout_owner': 'controller'},
         'stages': {'last_stage': 'native_complete', 'error': None, 'input_complete': True,
                    'gaussian_status': 0, 'gaussian_numerical_status': 0, 'design_status': 0, 'design_numerical_status': 0,
                    'design_rank': 3, 'equilibrated_triangular_condition_inf': 2.0, 'equilibrated_transpose_triangular_condition_inf': 2.0,
                    'design_method': 'retained-whitened-pivoted-householder-qr/v1', 'box_preparation_status': 0,
                    'box_preparation_numerical_status': None, 'box_preparation_numerical_status_scope': 'not exposed by current API',
                    'design_consumed_by_box': True},
         'payload_bounds': dict.fromkeys(('gaussian_preparation', 'design_preparation', 'box_preparation', 'box_evaluation'), 1024),
         'result_status': 0, 'result_numerical_status': 0, 'result_stage': 6,
         'method_id': 'retained-qr-rational-tail-box-enclosure/v1',
         'enclosure_scope': r['qualification']['enclosure_scope'], 'availability': dict.fromkeys(c.FLAGS, True),
         'completion': {'mean': means, 'variance': variances, 'named_variance_original46': variances[-1],
                        'minimum_quadratic': 2.0, 'reported_postcast_profile_quadratic': 2.0,
                        'log_design_precision_determinant': 0.0, 'source_covariance_log_determinant': 0.0,
                        'maximum_variance_sensitivity': 1e-12, 'profile_stationarity': 1e-12},
         'box_diagnostics': {'standardized_lower_margins': [m[0] for m in margins], 'standardized_upper_margins': [m[1] for m in margins],
                             'excluded_mass_upper': 1e-12, 'box_probability': {'lower': math.nextafter(1 - 1e-12, -math.inf), 'upper': 1.0},
                             'log_box_probability': {'lower': -2e-12, 'upper': 0.0}, 'log_prior_volume': enclosure(volume, 2e-11)},
         'normalizations': dict(zip(c.NORMALIZATIONS, (enclosure(relative, 2e-11), enclosure(prior, 2e-11), enclosure(observation, 2e-11)))),
         'median': {'original_index': 46, 'active_index': 45, 'cumulative_probability': 0.5, 'quantile': enclosure(means[-1], 2e-9),
                    'lower_endpoint_cdf': {'lower': 0.499999, 'upper': 0.5}, 'upper_endpoint_cdf': {'lower': 0.5, 'upper': 0.500001}},
         'work': {'cdf_node_evaluations': 1000, 'cdf_evaluations': 20, 'bisections': 10}, 'output_complete': True, 'accepted': True,
         'qualification': 'native conditional completion/enclosure only; independent original-input comparison unassessed; observational qualification blocked'}
    return copy.deepcopy(o)


def reference_fixture(native=None):
    native = native or native_fixture()
    scalar = lambda value: str(c.finite(value, 'fixture'))
    bounds = lambda value: {key: scalar(value[key]) for key in ('lower', 'upper')}
    result = {'digits': 120, 'correction_steps': 3,
              'mean': [scalar(v) for v in native['completion']['mean']],
              'variance': [scalar(v) for v in native['completion']['variance']],
              'q_min': '2', 'q_postcast': '2', 'logdet_H': '0', 'logdet_C': '0', 'excluded_mass_upper': '1e-12',
              'q_postcast_unrounded_residual': '2', 'q_postcast_residual_rounding_allowance': '0',
              'reduced_forward_sensitivity': '1e-12',
              'endpoint_lower_margins': [scalar(v['lower']) for v in native['box_diagnostics']['standardized_lower_margins']],
              'endpoint_upper_margins': [scalar(v['upper']) for v in native['box_diagnostics']['standardized_upper_margins']],
              **{key: bounds(native['box_diagnostics'][key]) for key in ('box_probability', 'log_box_probability', 'log_prior_volume')},
              'log_Z_relative_box': bounds(native['normalizations']['log_relative_box_integral']),
              **{key: bounds(native['normalizations'][key]) for key in ('log_prior_normalized_relative_evidence', 'log_observation_normalized_evidence')},
              'median_original46': bounds(enclosure(native['completion']['mean'][-1], 1e-10))}
    coarse = copy.deepcopy(result)
    coarse.update(digits=80, correction_steps=1)
    peer = {key: result[key] for key in ('mean', 'variance', 'q_min', 'logdet_H')}
    peer.update(rank=46, condition2='2')
    runtime = {'synthetic_verified_runtime': True}
    identity = {'runtime': runtime, 'request_sha256': c.REQUEST_SHA256,
                'reference_input_sha256': 'a' * 64, 'reference_script_sha256': 'b' * 64}
    original = REQUEST['lineage']
    def residual(steps, length=47, digits=None):
        arrays = {key: ['0'] * length for key in ('residual_infinity_by_rhs', 'backward_error_by_rhs', 'estimated_forward_error_by_rhs')}
        record = {'correction_steps': steps, **arrays, 'inverse_infinity_norm_estimate': '1',
                  'component_estimates': [{'component_index': index, 'rows': rows, 'inverse_infinity_norm_estimate': '1', **arrays}
                                         for index, rows in enumerate(c.COMPONENT_SIZES)]}
        if digits is not None:
            record['completion_digits'] = digits
        return record
    reference = {'schema_version': 1, 'interface_id': 'released-fixed44-box-reference/v1', 'status': 'accepted',
                 'accepted': True, 'error': None, 'blockers': [],
                 'identities': {key: identity[key] for key in ('request_sha256', 'reference_input_sha256', 'reference_script_sha256')},
                 'target': c.reference_target(REQUEST), 'settings': c.REFERENCE_SETTINGS,
                 'runtime_before': runtime, 'runtime_after': runtime,
                 'source_structure': {'inventory_sha256': 'e41d7f45428056e8c05e2952dbac54020009422e9a9e4c583aa4158e7eb2c2b1',
                                      'component_count': 13, 'maximum_component_size': 2593, 'maximum_general_component_size': 2593,
                                      'component_sizes': [2593, 55, 339, 143, 354] + [1] * 8,
                                      'dense_factor_multiply_add_pairs_per_precision': 2913157849},
                 'diagnostics': {'ancestry': 'synthetic schema fixture; no released computation',
                                 'original_input_LAPACK_SVD': copy.deepcopy(peer), 'original_input_LAPACK_pivoted_QR': copy.deepcopy(peer),
                                 'covariance_condition_estimate': {'norm_infinity': '2', 'reciprocal_condition_estimate': '0.5',
                                                                    'inverse_infinity_norm_estimate': '1',
                                                                    'estimation_scope': 'empirical-LAPACK-pocon-not-certified-bound/v1',
                                                                    'components': [{'component_index': index, 'rows': rows, 'norm_infinity': '2',
                                                                                    'reciprocal_condition_estimate': '0.5', 'inverse_infinity_norm_estimate': '1'}
                                                                                   for index, rows in enumerate(c.COMPONENT_SIZES)]},
                                 'wide_residual_history': [residual(steps) for steps in range(4)],
                                 'postcast_residual_history': [residual(steps, length=1, digits=digits) for digits in (80, 120) for steps in range(4)],
                                 'logdet_C_peer': {'cholesky': '0', 'eigenvalue': '0', 'absolute_difference': '0', 'component_count': 13}},
                 'coarse': coarse, 'fine': result,
                 'refinement': {'estimation_scope': 'empirical-original-input-calibration-not-certificate/v1',
                                'mean_absolute_errors': ['1e-12'] * 46, 'variance_absolute_errors': ['1e-16'] * 46,
                                **{key + '_absolute_error': '3e-10' if key == 'median_original46' else '1e-10'
                                   for key in c.REFINEMENT_SCALARS}},
                 'gates': dict.fromkeys(c.REFERENCE_GATES, True),
                 'work': {'elapsed_seconds': 1.0, 'cpu_user_seconds': 1.0, 'cpu_system_seconds': 0.0, 'maximum_rss_kib': 1024,
                          'covariance_factorizations': 1, 'covariance_eigenvalue_components': 13, 'wide_residual_evaluations': 4,
                          'correction_steps_completed': 3, 'reduced_decimal_factorizations': 2, 'original_input_rows': 3492,
                          'postcast_covariance_solves': 2, 'postcast_wide_residual_evaluations': 8, 'postcast_correction_steps_completed': 6,
                          'native_outputs_consumed': False},
                 'qualification': {'scope': 'empirical-original-input-comparison/v1', 'original_input_certificate': False,
                                   'conditional_enclosure_scope': 'conditional-Gaussian-union-mass-and-median-only/v1',
                                   'native_outputs_consumed': False, 'observational_qualification': 'blocked-original-contract-source-gaps',
                                   'limits': ['synthetic fixture; empirical method only']}}
    reference['identities'].update({'contract_sha256': c.CONTRACT_SHA256, 'lineage_sha256': original['packet_lineage_sha256'],
                                    'constrained_sha256': original['constrained_sha256'], 'original_source_sha256': c.ORIGINAL_SOURCE_SHA256,
                                    **{'canonical_' + item['name'][0] + '_sha256': item['sha256'] for item in REQUEST['transport']['canonical_inputs']}})
    return copy.deepcopy(reference), identity


class BoxControllerTests(unittest.TestCase):
    def test_request_change_is_not_admitted(self):
        blob = (CONTROLLER.parent / 'request.json').read_bytes()
        c.request_bytes(blob)
        for changed in (blob + b' ', blob.replace(b'"jobs": 1', b'"jobs": 2'), blob.replace(b'original44', b'original45')):
            with self.subTest(changed=changed[-20:]), self.assertRaises(ValueError):
                c.request_bytes(changed)

    def test_closed_native_success_fixture(self):
        self.assertTrue(c.check_native(native_fixture(), REQUEST, GUIDE))
        overlapping = native_fixture()
        overlapping['median']['lower_endpoint_cdf'] = {'lower': 0.49, 'upper': 0.51}
        overlapping['median']['upper_endpoint_cdf'] = {'lower': 0.49, 'upper': 0.51}
        self.assertTrue(c.check_native(overlapping, REQUEST, GUIDE))

    def test_reference_complete_fixture_and_exact_refinement(self):
        native = native_fixture()
        reference, identity = reference_fixture(native)
        self.assertTrue(all(check['passed'] for check in c.compare(native, reference, REQUEST, identity)))
        with localcontext() as context:
            context.prec = 180
            reference['coarse']['mean'][0] = str(Decimal(reference['fine']['mean'][0]) + Decimal('1e-40'))
        checks = c.compare(native, reference, REQUEST, identity)
        check = next(row for row in checks if row['label'] == 'reference/mean/original0/aggregate-covers-measured')
        self.assertEqual(Decimal(check['error_exact']), Decimal('1e-40'))

    def test_false_reference_and_unearned_refinement_refuse(self):
        native = native_fixture()
        paths = [(('accepted',), False), (('gates', 'original_input_empirical_completion'), False),
                 (('settings', 'digits'), [60, 90]), (('target', 'active_original_indices'), list(range(46))),
                 (('runtime_after',), {'changed': True}), (('coarse',), None),
                 (('fine', 'mean'), ['0'] * 45), (('fine', 'variance'), ['NaN'] * 46),
                 (('fine', 'excluded_mass_upper'), '0'), (('coarse', 'excluded_mass_upper'), '0.1'),
                 (('fine', 'log_box_probability'), {'lower': '0', 'upper': '0.0001'}),
                 (('fine', 'median_original46'), {'lower': '-100', 'upper': '-99'}),
                 (('work', 'cpu_user_seconds'), 900.5),
                 (('refinement', 'mean_absolute_errors'), ['0'] * 46),
                 (('refinement', 'log_Z_relative_box_absolute_error'), '1e-6')]
        for path, value in paths:
            reference, identity = reference_fixture(native)
            if path == ('refinement', 'mean_absolute_errors'):
                reference['coarse']['mean'][0] = str(Decimal(reference['fine']['mean'][0]) + Decimal('1e-9'))
            if path == ('work', 'cpu_user_seconds'):
                reference['work']['cpu_system_seconds'] = 1.0
            target = reference
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = value
            with self.subTest(path=path), self.assertRaises(ValueError):
                c.compare(native, reference, REQUEST, identity)

    def test_component_order_residual_work_and_postcast_rounding_are_admitted(self):
        changes = [(('diagnostics', 'covariance_condition_estimate', 'components', 0, 'rows'), 2592),
                   (('diagnostics', 'wide_residual_history', 0, 'component_estimates', 0, 'component_index'), 1),
                   (('diagnostics', 'wide_residual_history', 0, 'component_estimates', 0, 'estimated_forward_error_by_rhs'), ['1'] * 47),
                   (('diagnostics', 'postcast_residual_history', 4, 'completion_digits'), 80),
                   (('diagnostics', 'postcast_residual_history', 4, 'residual_infinity_by_rhs'), ['0'] * 47),
                   (('work', 'postcast_covariance_solves'), 1), (('work', 'postcast_wide_residual_evaluations'), 7),
                   (('fine', 'q_postcast_unrounded_residual'), '2.000001'),
                   (('fine', 'q_postcast_residual_rounding_allowance'), '1e-8')]
        for path, value in changes:
            native = native_fixture()
            reference, identity = reference_fixture(native)
            target = reference
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = value
            with self.subTest(path=path), self.assertRaises(ValueError):
                c.compare(native, reference, REQUEST, identity)

    def test_refused_reference_retains_earned_groups_without_acceptance(self):
        reference, identity = reference_fixture()
        reference.update(accepted=False, status='refused', error='synthetic fine solve refusal',
                         blockers=['synthetic fine solve refusal'], fine=None, refinement=None)
        reference['gates'] = dict.fromkeys(c.REFERENCE_GATES, False)
        reference['gates'].update(source_identity=True, runtime_identity=True)
        reference['diagnostics']['postcast_residual_history'] = reference['diagnostics']['postcast_residual_history'][:4]
        reference['work'].update(reduced_decimal_factorizations=1, postcast_covariance_solves=1,
                                 postcast_wide_residual_evaluations=4, postcast_correction_steps_completed=3)
        earned = copy.deepcopy(reference['coarse'])
        self.assertFalse(c.check_reference(reference, REQUEST, identity, require_complete=False))
        self.assertEqual(reference['coarse'], earned)
        with self.assertRaises(ValueError):
            c.compare(native_fixture(), reference, REQUEST, identity)
        admission = copy.deepcopy(reference)
        admission.update(target=None, identities={'request_sha256': c.REQUEST_SHA256,
                                                   'original_source_sha256': {'MCMC_utils.py': c.ORIGINAL_SOURCE_SHA256['MCMC_utils.py']}},
                         source_structure=None, coarse=None, runtime_before=None, runtime_after=None)
        admission['diagnostics'] = {'ancestry': 'source admission refusal fixture', 'wide_residual_history': [], 'postcast_residual_history': []}
        for key in ('covariance_factorizations', 'covariance_eigenvalue_components', 'wide_residual_evaluations',
                    'correction_steps_completed', 'reduced_decimal_factorizations', 'postcast_covariance_solves',
                    'postcast_wide_residual_evaluations', 'postcast_correction_steps_completed'):
            admission['work'][key] = 0
        self.assertFalse(c.check_reference(admission, REQUEST, identity, require_complete=False))
        admission['identities']['unknown_source'] = 'a' * 64
        with self.assertRaises(ValueError):
            c.check_reference(admission, REQUEST, identity, require_complete=False)

    def test_native_discrepancy_and_reference_error_share_total_allocation(self):
        families = [('mean', ('completion', 'mean', 0), 'mean_absolute_errors', 0),
                    ('variance', ('completion', 'variance', 0), 'variance_absolute_errors', 0),
                    ('q_min', ('completion', 'minimum_quadratic'), 'q_min_absolute_error', None),
                    ('q_postcast', ('completion', 'reported_postcast_profile_quadratic'), 'q_postcast_absolute_error', None),
                    ('logdet_H', ('completion', 'log_design_precision_determinant'), 'logdet_H_absolute_error', None),
                    ('logdet_C', ('completion', 'source_covariance_log_determinant'), 'logdet_C_absolute_error', None),
                    ('log_prior_volume', ('box_diagnostics', 'log_prior_volume'), 'log_prior_volume_absolute_error', None),
                    ('log_Z_relative_box', ('normalizations', 'log_relative_box_integral'), 'log_Z_relative_box_absolute_error', None),
                    ('log_prior_normalized_relative_evidence', ('normalizations', 'log_prior_normalized_relative_evidence'),
                     'log_prior_normalized_relative_evidence_absolute_error', None),
                    ('log_observation_normalized_evidence', ('normalizations', 'log_observation_normalized_evidence'),
                     'log_observation_normalized_evidence_absolute_error', None),
                    ('median_original46', ('median', 'quantile'), 'median_original46_absolute_error', None)]
        for key, path, error_key, index in families:
            native = native_fixture()
            reference, identity = reference_fixture(native)
            baseline = c.compare(native, reference, REQUEST, identity)
            label = 'native/' + key + ('/original0' if index is not None else '')
            budget = Decimal(next(row for row in baseline if row['label'] == label)['budget_exact'])
            with localcontext() as context:
                context.prec = 180
                reported = Decimal('0.04') * budget
                if index is None:
                    reference['refinement'][error_key] = str(reported)
                else:
                    reference['refinement'][error_key][index] = str(reported)
                target = native
                for part in path[:-1]:
                    target = target[part]
                if isinstance(target[path[-1]], dict):
                    value = reference['fine'][key]
                    target[path[-1]] = {edge: float(Decimal(value[edge]) + Decimal('0.98') * budget)
                                       for edge in ('lower', 'upper')}
                else:
                    target[path[-1]] = float(c.finite(target[path[-1]], key) + Decimal('0.98') * budget)
            with self.subTest(key=key), self.assertRaises(c.ComparisonFailure) as raised:
                c.compare(native, reference, REQUEST, identity)
            row = next(row for row in raised.exception.checks if row['label'] == label)
            self.assertLess(Decimal(row['discrepancy_exact']), budget)
            self.assertGreater(Decimal(row['error_exact']), budget)
            self.assertFalse(row['passed'])

    def test_packet_remains_blocked_conditional_control(self):
        packet = {'id': 'released-box', 'origin': {'kind': 'development_smoke'}, 'status': 'blocked', 'execution': None}
        with patch.object(c.receipt, 'read_packet', return_value=(packet, None, {'packet': 'a' * 64})):
            self.assertEqual(c.check_packet(Path('fixture'))[0], packet)
        for field, value in (('id', 'released-ladder'), ('origin', {'kind': 'prospector_candidate'}),
                             ('status', 'runnable'), ('execution', {'operation': 'guessed'})):
            changed = {**packet, field: value}
            with self.subTest(field=field), patch.object(c.receipt, 'read_packet', return_value=(changed, None, {})), self.assertRaises(ValueError):
                c.check_packet(Path('fixture'))

    def test_axis_policy_status_and_numeric_changes_refuse(self):
        def change(path, value):
            o = native_fixture()
            target = o
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = value
            return o
        changes = [(('target', 'active_original_indices'), list(range(46))),
                   (('target', 'support_lower'), list(reversed(c.support(REQUEST)[0]))),
                   (('target', 'fixed_coordinates'), [{'original_index': 44, 'value': 1, 'measure': 'point mass outside active46 Lebesgue'}]),
                   (('median', 'active_index'), 46), (('median', 'cumulative_probability'), True),
                   (('policy', 'maximum_total_cdf_nodes'), 2097153), (('resources', 'threads'), 2),
                   (('stages', 'design_rank'), 2), (('result_numerical_status',), 4), (('accepted',), 1),
                   (('availability', 'endpoint_margins'), False), (('completion', 'mean'), [0.0] * 45),
                   (('completion', 'variance'), [float('nan')] * 46), (('completion', 'minimum_quadratic'), float('inf')),
                   (('completion', 'maximum_variance_sensitivity'), 1e-9), (('work', 'cdf_evaluations'), 257),
                   (('box_diagnostics', 'excluded_mass_upper'), 0.0),
                   (('box_diagnostics', 'log_prior_volume'), {'lower': 0.0, 'upper': 0.0}),
                   (('normalizations', 'log_relative_box_integral'), {'lower': 0.0, 'upper': 0.0}),
                   (('median', 'quantile'), {'lower': 9.0, 'upper': 9.0001}),
                   (('median', 'lower_endpoint_cdf'), {'lower': -1.0, 'upper': 0.5}),
                   (('median', 'lower_endpoint_cdf'), {'lower': 0.51, 'upper': 0.52}),
                   (('median', 'upper_endpoint_cdf'), {'lower': 0.48, 'upper': 0.49})]
        for path, value in changes:
            with self.subTest(path=path), self.assertRaises(ValueError):
                c.check_native(change(path, value), REQUEST, GUIDE)
        o = native_fixture()
        o['unknown'] = 0
        with self.assertRaises(ValueError):
            c.check_native(o, REQUEST, GUIDE)

    def test_partial_groups_are_retained_without_acceptance(self):
        o = native_fixture()
        original_completion = copy.deepcopy(o['completion'])
        o.update(result_status=4, result_numerical_status=8, result_stage=2, output_complete=False, accepted=False)
        o['stages']['last_stage'] = 'native_refusal_or_output_postcondition'
        o['availability'] = {flag: index <= 2 for index, flag in enumerate(c.FLAGS, 1)}
        o['normalizations'] = o['median'] = None
        o['box_diagnostics']['box_probability'] = o['box_diagnostics']['log_box_probability'] = None
        self.assertFalse(c.check_native(o, REQUEST, GUIDE, require_complete=False))
        self.assertEqual(o['completion'], original_completion)
        with self.assertRaises(ValueError):
            c.check_native(o, REQUEST, GUIDE)
        o['normalizations'] = dict.fromkeys(c.NORMALIZATIONS, {'lower': 0, 'upper': 0})
        with self.assertRaises(ValueError):
            c.check_native(o, REQUEST, GUIDE, require_complete=False)

    def test_changed_raw_and_canonical_bytes_are_refused(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            blob = b'original source'
            (directory / 'raw').write_bytes(blob)
            manifest = {'sources': [{'name': 'raw', 'bytes': len(blob), 'sha256': hashlib.sha256(blob).hexdigest(), 'role': 'test'}]}
            c.raw_sources(directory, manifest)
            (directory / 'raw').write_bytes(b'changed! source')
            with self.assertRaises(ValueError):
                c.raw_sources(directory, manifest)
            (directory / 'raw').unlink()
            (directory / 'outside').write_bytes(blob)
            (directory / 'raw').symlink_to(directory / 'outside')
            with self.assertRaises(ValueError):
                c.raw_sources(directory, manifest)
            canonical_blob = b'12345678'
            (directory / 'C.f64').write_bytes(canonical_blob)
            request = {'transport': {'canonical_inputs': [{'name': 'C.f64', 'bytes': 8,
                                                           'sha256': hashlib.sha256(canonical_blob).hexdigest()}]}}
            c.canonical_identity(directory, request)
            (directory / 'C.f64').write_bytes(b'92345678')
            with self.assertRaises(ValueError):
                c.canonical_identity(directory, request)

    def test_runtime_extension_and_library_closure_and_drift(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            for package in ('numpy', 'scipy'):
                (root / package).mkdir()
                (root / (package + '.libs')).mkdir()
                (root / package / 'extension.so').write_bytes(b'synthetic extension')
                (root / (package + '.libs') / 'blas.so').write_bytes(b'synthetic blas')
            executable = root / 'python'
            executable.write_bytes(b'synthetic interpreter')
            decimal = root / 'decimal.py'
            decimal.write_bytes(b'synthetic decimal source')
            inventory = [{'path': str(p), 'bytes': p.stat().st_size, 'sha256': c.sha256(p)}
                         for p in sorted(root.rglob('*'), key=str) if p.is_file()]
            runtime = {'python': {'version': 'synthetic CPython', 'implementation': 'CPython',
                                  'executable': {'path': str(executable), 'resolved_path': str(executable),
                                                 'bytes': executable.stat().st_size, 'sha256': c.sha256(executable)}},
                       'packages': {**{p: {'version': '2.5.3' if p == 'numpy' else '1.18.1', 'module_root': str(root / p)} for p in ('numpy', 'scipy')},
                                    'decimal': {'precision_policy': [80, 120], 'module_path': str(decimal)}},
                       'thread_environment': dict.fromkeys(('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
                                                            'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'), '1'),
                       'inventory': inventory,
                       'inventory_sha256': hashlib.sha256(json.dumps(inventory, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
                       'numeric_configuration': json.dumps({'show_config': 'synthetic library inventory fixture',
                           'threadpools': [{'filepath': str(root / 'numpy.libs/blas.so'), 'num_threads': 1},
                                           {'filepath': str(root / 'scipy.libs/blas.so'), 'num_threads': 1}],
                           'decimal': {'version': '1.70', 'libmpdec_version': '4.0.0', 'rounding': 'ROUND_HALF_EVEN'}})}
            c.verify_runtime(runtime, executable)
            for change in ('threads', 'library', 'rounding', 'empty', 'nonfinite'):
                changed = copy.deepcopy(runtime)
                configuration = json.loads(changed['numeric_configuration'])
                if change == 'threads':
                    configuration['threadpools'][0]['num_threads'] = 2
                elif change == 'library':
                    configuration['threadpools'][0]['filepath'] = str(root / 'unlisted.so')
                elif change == 'rounding':
                    configuration['decimal']['rounding'] = 'ROUND_DOWN'
                elif change == 'empty':
                    configuration['threadpools'] = []
                else:
                    configuration['threadpools'][0]['num_threads'] = float('nan')
                changed['numeric_configuration'] = json.dumps(configuration)
                with self.subTest(change=change), self.assertRaises(ValueError):
                    c.verify_runtime(changed, executable)
            omitted = copy.deepcopy(runtime)
            omitted['inventory'] = [i for i in omitted['inventory'] if not i['path'].endswith('extension.so')]
            omitted['inventory_sha256'] = hashlib.sha256(json.dumps(omitted['inventory'], sort_keys=True, separators=(',', ':')).encode()).hexdigest()
            with self.assertRaises(ValueError):
                c.verify_runtime(omitted, executable)
            (root / 'scipy.libs/blas.so').write_bytes(b'changed! blas!!')
            with self.assertRaises(ValueError):
                c.verify_runtime(runtime, executable)

    def test_complete_sdk_inventory_includes_cmake_and_changed_actual_bytes(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            (root / 'src').mkdir()
            (root / 'cpp/cmake').mkdir(parents=True)
            (root / 'cpp/include/irred').mkdir(parents=True)
            (root / 'docs').mkdir()
            (root / 'build/native-release').mkdir(parents=True)
            (root / 'target/release').mkdir(parents=True)
            sources = [root / 'src' / (str(i) + '.rs') for i in range(305)] + [
                root / 'Cargo.toml', root / 'Cargo.lock', root / 'build.rs', root / 'cpp/CMakeLists.txt',
                root / 'cpp/cmake/gaussian_box.cmake', root / 'cpp/include/irred/gaussian_box.hpp']
            for path in sources:
                path.write_bytes(b'fixed synthetic SDK source')
            guide = root / 'docs/gaussian-box.md'
            guide.write_bytes(b'fixed synthetic guide')
            library = root / 'build/stdlib.so'
            library.write_bytes(b'fixed synthetic standard library')
            archive = root / 'build/native-release/libirred_core.a'
            archive.write_bytes(b'fixed synthetic archive')
            cli = root / 'target/release/irred'
            cli.write_bytes(b'fixed synthetic CLI')
            compiler_digest = c.sha256(Path('/usr/bin/c++').resolve())
            build = {'profile': 'release', 'backend': 'portable_cpu',
                     'sources': {str(path.relative_to(root)): c.sha256(path) for path in sources},
                     'compiler_executable_digest': compiler_digest, 'standard_library': str(library),
                     'standard_library_digest': c.sha256(library), 'tool_executable_digests': {}}
            build_id = hashlib.sha256(json.dumps(build, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
            build.update(build_id=build_id, git_head=REQUEST['sdk_identity']['engine_revision'], git_status='')
            manifest = root / 'build/build-manifest-release.json'
            manifest.write_text(json.dumps(build))
            request = copy.deepcopy(REQUEST)
            request['sdk_identity'].update(build_id=build_id, manifest_sha256=c.sha256(manifest), archive_sha256=c.sha256(archive),
                                           cli_sha256=c.sha256(cli), gaussian_box_header_sha256=c.sha256(sources[-1]),
                                           compiler_executable_sha256=compiler_digest, standard_library_sha256=c.sha256(library))
            def git(root, *args):
                return request['sdk_identity']['engine_revision'] if args[0] == 'rev-parse' else ''
            with patch.object(c.receipt, 'git', git), patch.object(c.subprocess, 'check_output', return_value=guide.read_bytes()):
                identity = c.engine_identity(root, root, request)
                self.assertEqual(len(identity['sources']), 311)
                self.assertIn('cpp/cmake/gaussian_box.cmake', identity['sources'])
                sources[-2].write_bytes(b'changed synthetic SDK source')
                with self.assertRaises(ValueError):
                    c.engine_identity(root, root, request)
                sources[-2].write_bytes(b'fixed synthetic SDK source')
                (root / 'cpp/include/unlisted.hpp').write_bytes(b'unlisted header')
                with self.assertRaises(ValueError):
                    c.engine_identity(root, root, request)
                (root / 'cpp/include/unlisted.hpp').unlink()
                archive.write_bytes(b'changed synthetic archive')
                with self.assertRaises(ValueError):
                    c.engine_identity(root, root, request)

    def test_failed_and_timeout_child_preserve_identity_and_partial_output(self):
        limits = {'memory_bytes': 1 << 30, 'output_bytes': 4096, 'native_seconds': 1}
        for timeout in (False, True):
            with self.subTest(timeout=timeout), tempfile.TemporaryDirectory() as name:
                store, record = Path(name), {}
                def run(command, **kwargs):
                    kwargs['stdout'].write(b'{"partial":true}')
                    kwargs['stderr'].write(b'refused\n')
                    if timeout:
                        raise subprocess.TimeoutExpired(command, 1)
                    return SimpleNamespace(returncode=3)
                with patch.object(c.receipt.subprocess, 'run', run), self.assertRaises(ValueError):
                    c.child(['synthetic'], store, 'native', limits, record=record)
                operation = record['subprocesses']['native']
                self.assertEqual(operation['timed_out'], timeout)
                self.assertEqual(operation['out_sha256'], hashlib.sha256(b'{"partial":true}').hexdigest())
                errors, stdout = c.receipt.verify_outputs(store, record)
                self.assertEqual(errors, [])
                c.receipt.ingest_outputs(store, record, stdout)
                self.assertEqual(record['native'], {'partial': True})

    def test_terminal_raw_drift_revokes_nominal_pass_without_overwriting(self):
        for suffix, replacement in (('out', b'{"changed":true}'), ('out', b'malformed'), ('out', None), ('err', b'changed error')):
            with self.subTest(suffix=suffix, replacement=replacement), tempfile.TemporaryDirectory() as name:
                store = Path(name)
                record = {'status': 'completed', 'gates': {'numerical': 'passed'}, 'native': {'original': True}, 'comparisons': [{'passed': True}],
                          'subprocesses': {'native': {'out_bytes': 17, 'out_sha256': hashlib.sha256(b'{"original":true}').hexdigest(),
                                                      'err_bytes': 0, 'err_sha256': hashlib.sha256(b'').hexdigest()}}}
                (store / 'native.out').write_bytes(b'{"original":true}')
                (store / 'native.err').write_bytes(b'')
                path = store / ('native.' + suffix)
                path.unlink() if replacement is None else path.write_bytes(replacement)
                errors = []
                c.terminal_outputs(store, record, errors)
                self.assertEqual(record['status'], 'failed')
                self.assertEqual(record['native'], {'original': True})
                self.assertEqual(record['comparisons'], [{'passed': True}])
                self.assertTrue(any('raw log drift' in error for error in errors))
                self.assertEqual(record['gates']['numerical'], 'not accepted: identity failure')

    def test_names_and_fresh_failed_attempts(self):
        with self.assertRaises(ValueError):
            c.execute(SimpleNamespace(name='../escape'))
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            args = SimpleNamespace(name='failed', engine_source=root, engine_artifacts=root, sources=root)
            with patch.object(c, 'ROOT', root), patch.object(c.receipt, 'source_identity', side_effect=ValueError('synthetic admission failure')):
                self.assertEqual(c.execute(args), 1)
                record = root / 'results/released-box/failed/record.json'
                self.assertIn('synthetic admission failure', json.loads(record.read_text())['error'])
                self.assertEqual(record.stat().st_mode & 0o222, 0)
                original = record.read_bytes()
                with self.assertRaises(FileExistsError):
                    c.execute(args)
                self.assertEqual(record.read_bytes(), original)


if __name__ == '__main__':
    unittest.main()
