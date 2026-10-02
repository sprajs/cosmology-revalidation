#!/usr/bin/env python3
"""One bounded conditional thermal/DESI comparison within the blocked campaign."""
import argparse
from datetime import datetime, timezone
from decimal import Decimal, localcontext
import hashlib
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
sys.path.insert(0, str(ROOT / 'scripts'))
from packet import load, parse, sha256, within, read_packet

REQUEST_SHA256 = '0e6fca28ec07187da6bcc91c0ae917d0b3caad7d16f91ceae9721eab2f937d15'
CANDIDATE_ORIGIN = {
    'kind': 'prospector_candidate', 'repository': 'sprajs/prospector',
    'revision': '023b3839b2fd4fd55e737fbf1415d3b3360c6ed8',
    'path': 'designs/candidate-lcdm-baseline-reference-audit.json',
    'snapshot': 'candidate.json',
    'sha256': 'b098fa4423e2a228ddf392e621ed0bff6748bc094b2bc49696e83033719f45a1',
}
METHOD = 'BAO/flat-thermal-FD-supplied-drag-normalized-ratio-density/v1'
MAPPING = 'physical-omega-h2-explicit-Kelvin-temperature-state-weight-SI/v1'
MODEL_ORDER = ['thermal-massless-supplied-drag', 'thermal-massive-FD-supplied-drag']
INPUT_IDENTITIES = [
    {'path': 'data/bao/desi_gaussian_bao_ALL_GCcomb_mean.txt', 'bytes': 472,
     'sha256': '9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585',
     'source': 'https://raw.githubusercontent.com/CobayaSampler/bao_data/bb0c1c9009dc76d1391300e169e8df38fd1096db/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt'},
    {'path': 'data/bao/desi_gaussian_bao_ALL_GCcomb_cov.txt', 'bytes': 2547,
     'sha256': '252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509',
     'source': 'https://raw.githubusercontent.com/CobayaSampler/bao_data/bb0c1c9009dc76d1391300e169e8df38fd1096db/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt'},
]
DENSITY_KEYS = ('quadratic', 'log_determinant', 'normalization', 'log_density')
WORK_KEYS = ('callbacks', 'preparation_callbacks', 'outer_callbacks', 'momentum_callbacks')
ARITHMETIC = {'id': 'F02/longdouble-cpu/v1', 'double_mantissa_bits': 53,
              'long_double_mantissa_bits': 64, 'long_double_max_exponent': 16384,
              'round_to_nearest': True}
POLICY_METADATA = {
    'momentum_method': 'scaled-adaptive-direct-momentum-exponential-tail/v1',
    'provider_maximum_depth': 30, 'provider_maximum_points': 4096,
    'provider_maximum_native_bytes': 16777216,
    'momentum_maximum_callbacks_per_evaluation': 200000,
    'momentum_maximum_total_callbacks': 500000000, 'momentum_maximum_depth': 30,
    'momentum_maximum_points': 4096, 'momentum_maximum_species': 16,
    'momentum_maximum_native_bytes': 16777216,
    'density_maximum_models': 2, 'density_maximum_queries': 13,
    'density_maximum_string_bytes': 8192, 'density_maximum_native_bytes': 16777216,
    'density_maximum_forward_sensitivity': 1e-8, 'density_maximum_projection_log_density_error': 1e-8,
    'preparation_maximum_queries': 13, 'preparation_maximum_matrix_elements': 169,
    'preparation_maximum_string_bytes': 8192, 'preparation_maximum_native_bytes': 1048576,
    'preparation_maximum_forward_sensitivity': 1e-8,
}
REFERENCE_METHOD = 'reference/direct-momentum-GL-direct-z-sqrt-a-scalar-LDLT/v1'
REFERENCE_CONSTANTS = 'SI2019-exact-h-c-kB-eV-IAU2012-AU-CODATA2018-G-fixed'
REFERENCE_SETTINGS = {
    name: dict(zip(('digits', 'momentum_nodes_per_panel', 'outer_nodes_per_panel', 'momentum_tail'), values))
    for name, values in (('coarse', (60, 32, 32, 128)), ('fine', (90, 48, 48, 160)),
                         ('momentum_refined', (90, 48, 32, 160)), ('outer_refined', (90, 32, 48, 128)))
}
DECIMAL_PATTERN = re.compile(r'-?(?:0|[1-9]\d*)(?:\.\d+)?(?:[eE][+-]?\d+)?\Z')


class ComparisonFailure(ValueError):
    def __init__(self, message, checks):
        super().__init__(message)
        self.checks = list(checks)


def keys(value, expected, label):
    if not isinstance(value, dict) or set(value) != set(expected):
        raise ValueError(label + ' closed fields differ')


def finite(value, label):
    try:
        admitted = type(value) in (int, float) and math.isfinite(value)
    except OverflowError:
        admitted = False
    if not admitted:
        raise ValueError(label + ' requires finite binary64 numbers')
    return Decimal.from_float(value) if type(value) is float else Decimal(value)


def count(value, maximum, label):
    if type(value) is not int or not 0 <= value <= maximum:
        raise ValueError(label + ' integer work/status domain differs')
    return value


def decimal_string(value, label):
    if not isinstance(value, str) or len(value) > 256 or not DECIMAL_PATTERN.fullmatch(value):
        raise ValueError(label + ' requires a finite decimal string')
    result = Decimal(value)
    if not result.is_finite() or abs(result.adjusted()) > 1000:
        raise ValueError(label + ' decimal domain differs')
    return result


def exact(value, expected, label):
    """JSON equality that does not admit True as 1 or discard numeric types."""
    if isinstance(expected, dict):
        keys(value, expected, label)
        for key in expected:
            exact(value[key], expected[key], label + '/' + key)
    elif isinstance(expected, list):
        if not isinstance(value, list) or len(value) != len(expected):
            raise ValueError(label + ' ordered length differs')
        for i, (actual, wanted) in enumerate(zip(value, expected)):
            exact(actual, wanted, label + '/' + str(i))
    elif type(expected) is float:
        if finite(value, label) != Decimal.from_float(expected):
            raise ValueError(label + ' exact binary64 value differs')
    elif type(value) is not type(expected) or value != expected:
        raise ValueError(label + ' exact value differs')


def git(root, *arguments):
    return subprocess.check_output(['git', '-C', str(root), *arguments], text=True,
                                   timeout=30).strip()


def validate_request_bytes(blob):
    if len(blob) > 65536 or hashlib.sha256(blob).hexdigest() != REQUEST_SHA256:
        raise ValueError('exact reviewed request SHA-256 differs; review a new design')
    request = parse(blob.decode('utf-8'))
    if request['resources']['jobs'] != 1 or request['resources']['threads'] != 1:
        raise ValueError('one job and one thread required')
    return request


def validate_admission(packet):
    exact(packet.get('origin'), CANDIDATE_ORIGIN, 'immutable candidate origin')
    if packet.get('id') != 'lcdm-campaign' or packet.get('status') != 'blocked' or packet.get('execution') is not None:
        raise ValueError('blocked full campaign and separate conditional control required')
    input_sources(packet)


def source_identity():
    head = git(ROOT, 'rev-parse', 'HEAD')
    status = git(ROOT, 'status', '--porcelain', '--untracked-files=all')
    paths = git(ROOT, 'ls-files', '-z').split('\0')
    paths = [name for name in paths if name]
    hashes = {name: sha256(within(ROOT, name)) if within(ROOT, name).is_file() else 'missing' for name in paths}
    return {'head': head, 'status': status, 'files': hashes}


def snapshot_source(store, identity, record):
    if identity['status']:
        raise ValueError('final attempts require clean committed Reproducible source')
    child(['git', '-C', str(ROOT), 'archive', '--format=tar',
           '--output=' + str(store / 'source.tar'), identity['head']],
          store, 'source', {'source_seconds': 30, 'memory_bytes': 1073741824,
                           'output_bytes': 2097152}, record=record)
    destination = store / 'source'
    destination.mkdir()
    seen = set()
    total = 0
    with tarfile.open(store / 'source.tar', 'r:') as archive:
        for member in archive:
            target = within(destination, member.name)
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            if not member.isfile() or member.size > 1048576 or member.name in seen:
                raise ValueError('committed snapshot type/size/duplicate rejected')
            total += member.size
            if total > 1048576:
                raise ValueError('committed source snapshot exceeds tree quota')
            blob = archive.extractfile(member).read(member.size + 1)
            if (len(blob) != member.size or hashlib.sha256(blob).hexdigest() !=
                    identity['files'].get(member.name) or within(ROOT, member.name).read_bytes() != blob):
                raise ValueError('working source differs from immutable committed snapshot')
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(blob)
            seen.add(member.name)
    if seen != set(identity['files']):
        raise ValueError('committed source snapshot inventory differs')
    record['source_archive_sha256'] = sha256(store / 'source.tar')
    return destination


def fingerprint(engine, identity, artifacts=None):
    artifacts = artifacts or engine
    head = git(engine, 'rev-parse', 'HEAD')
    status = git(engine, 'status', '--porcelain', '--untracked-files=all')
    if head != identity['revision'] or status:
        raise ValueError('clean pinned engine revision required')
    paths = {'manifest': artifacts / 'build/build-manifest-release.json',
             'archive': artifacts / 'build/native-release/libirred_core.a',
             'cli': artifacts / 'target/release/irred'}
    hashes = {key: sha256(path) for key, path in paths.items()}
    for key, digest in hashes.items():
        if digest != identity[key + '_sha256']:
            raise ValueError('changed engine ' + key)
    build = load(paths['manifest'])
    if (build['git_head'] != head or build['git_status'] != '' or
            build['build_id'] != identity['build_id'] or build['profile'] != 'release' or
            build['backend'] != 'portable_cpu'):
        raise ValueError('build/source identity mismatch')
    content = {key: value for key, value in build.items()
               if key not in ('build_id', 'git_head', 'git_status')}
    if hashlib.sha256(json.dumps(content, sort_keys=True, separators=(',', ':')).encode()).hexdigest() != build['build_id']:
        raise ValueError('build manifest identity does not reconstruct')
    patterns = ('src/**/*.rs', 'tests/**/*.rs', 'cpp/**/*.cpp', 'cpp/**/*.h',
                'cpp/**/*.hpp', 'cpp/**/*.inc', 'schema/*.json', 'tools/*.py')
    inventory = {str(path.relative_to(engine)) for pattern in patterns for path in engine.glob(pattern)}
    inventory.update(('Cargo.toml', 'Cargo.lock', 'build.rs', 'cpp/CMakeLists.txt'))
    if inventory != set(build['sources']):
        raise ValueError('complete engine source/header inventory differs')
    sources = {}
    for name, digest in build['sources'].items():
        path = within(engine, name)
        if sha256(path) != digest:
            raise ValueError('changed build source ' + name)
        sources[name] = digest
    tools = {'compiler': Path('/usr/bin/c++').resolve(),
             'standard_library': Path(build['standard_library']).resolve()}
    expected = {'compiler': build['compiler_executable_digest'],
                'standard_library': build['standard_library_digest']}
    for name, digest in build['tool_executable_digests'].items():
        executable = shutil.which(name)
        if executable is None:
            raise ValueError('missing build tool ' + name)
        tools[name] = Path(executable).resolve()
        expected[name] = digest
    tool_identity = {}
    for name, path in tools.items():
        digest = sha256(path)
        if digest != expected[name]:
            raise ValueError('changed compiler/standard library/tool ' + name)
        tool_identity[name] = {'path': str(path), 'sha256': digest}
    return {**hashes, 'source_path': str(engine), 'artifact_path': str(artifacts),
            'sources': sources, 'tools': tool_identity, 'head': head, 'status': status,
            'build_id': build['build_id']}


def input_sources(packet):
    if not isinstance(packet.get('inputs'), list) or len(packet['inputs']) != 2:
        raise ValueError('exact two ordered original DESI sources required')
    for item, identity in zip(packet['inputs'], INPUT_IDENTITIES):
        keys(item, (*identity, 'role', 'semantics'), 'input declaration')
        exact({key: item[key] for key in identity}, identity, 'exact input identity')
        if item['role'] != 'released_fitted_summary' or not isinstance(item['semantics'], str) or not item['semantics']:
            raise ValueError('exact released fitted summary roles/semantics required')
    return [dict(identity) for identity in INPUT_IDENTITIES]


def inputs(input_root, packet):
    identities = input_sources(packet)
    blobs = []
    for item in identities:
        path = within(input_root, item['path'])
        blob = path.read_bytes()
        if len(blob) != item['bytes'] or hashlib.sha256(blob).hexdigest() != item['sha256']:
            raise ValueError('changed exact input ' + item['path'])
        blobs.append(blob)
    kinds = {'DM_over_rs': 0, 'DH_over_rs': 1, 'DV_over_rs': 2}
    rows = []
    for line in blobs[0].decode('utf-8').splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        z, observed, kind = line.split()
        rows.append({'z': float(z), 'observed': float(observed), 'kind': kinds[kind]})
    covariance = [list(map(float, line.split())) for line in blobs[1].decode('utf-8').splitlines()
                  if line.strip() and not line.lstrip().startswith('#')]
    if len(rows) != 13 or len(covariance) != 13 or any(len(row) != 13 for row in covariance):
        raise ValueError('full13 axis dimensions required')
    for row in rows:
        if finite(row['z'], 'source z') <= 0 or finite(row['observed'], 'source observed') <= 0:
            raise ValueError('positive source axes required')
    for row in covariance:
        for value in row:
            finite(value, 'source covariance')
    return blobs, rows, covariance


def queries(rows):
    return [{**row, 'id': 'DESI-DR2-row-' + str(i)} for i, row in enumerate(rows)]


def transport(request, rows, covariance):
    lines = [str(len(request['models']))]
    for model in request['models']:
        lines.append(' '.join(str(model[key]) for key in
                              ('H0', 'omega_b', 'omega_c', 'Tcmb', 'omega_other', 'z_drag')) +
                     ' ' + str(len(model['species'])))
        for species in model['species']:
            lines.append(' '.join(str(species[key]) for key in ('mass_eV', 'temperature_K', 'g')))
    lines.append(str(len(rows)))
    lines.extend(f"{row['z']} {row['kind']} {row['observed']}" for row in rows)
    lines.extend(' '.join(map(str, row)) for row in covariance)
    return '\n'.join(lines) + '\n'


def check_native(native, request, rows):
    keys(native, ('schema_version', 'method', 'mapping', 'requested', 'batch_status',
                  'numerical_status', 'arithmetic', 'producer_policy', 'policy_metadata', 'model_order',
                  'queries', 'model_sources', 'slots', *WORK_KEYS, 'accepted'), 'native')
    for name, expected in (('schema_version', 1), ('method', METHOD), ('mapping', MAPPING),
                           ('requested', 7), ('batch_status', 0), ('numerical_status', 0),
                           ('arithmetic', ARITHMETIC), ('producer_policy', request['producer_policy']),
                           ('policy_metadata', POLICY_METADATA),
                           ('model_order', MODEL_ORDER), ('model_sources', request['models']),
                           ('queries', queries(rows)), ('accepted', True)):
        exact(native[name], expected, 'native/' + name)
    if not isinstance(native['slots'], list) or len(native['slots']) != 2:
        raise ValueError('complete two-model native slots required')
    totals = dict.fromkeys(WORK_KEYS, 0)
    shared_density = None
    for index, slot in enumerate(native['slots']):
        keys(slot, ('source_index', 'model_source', 'numerical_status', 'preparation_status', 'predictions',
                    'residuals', 'predictions_state', 'residuals_state', 'density_state',
                    'density_status', 'projection_estimate', 'density', *WORK_KEYS), 'native slot')
        exact(slot['source_index'], index, 'native source index')
        exact(slot['model_source'], request['models'][index], 'native slot source model')
        for name in ('numerical_status', 'preparation_status', 'density_status'):
            exact(slot[name], 0, 'native slot/' + name)
        for name in ('predictions_state', 'residuals_state', 'density_state'):
            exact(slot[name], {'availability': 1, 'status': 0, 'numerical_status': 0}, name)
        for name in ('predictions', 'residuals'):
            if not isinstance(slot[name], list) or len(slot[name]) != 13:
                raise ValueError('native full13 ' + name + ' required')
            for value in slot[name]:
                finite(value, 'native ' + name)
        for row, prediction, residual in zip(rows, slot['predictions'], slot['residuals']):
            if prediction <= 0 or residual != row['observed'] - prediction:
                raise ValueError('native prediction/residual binary64 consistency differs')
        projection = finite(slot['projection_estimate'], 'native projection')
        if not 0 <= projection <= Decimal(str(request['budgets']['projection_absolute'])):
            raise ValueError('native projection budget failed')
        for name in WORK_KEYS:
            totals[name] += count(slot[name], request['producer_policy']['maximum_callbacks_per_model'], name)
        if (slot['callbacks'] == 0 or slot['outer_callbacks'] == 0 or
                slot['callbacks'] != slot['outer_callbacks'] + slot['momentum_callbacks'] or
                slot['preparation_callbacks'] > slot['momentum_callbacks'] or
                (index == 0 and (slot['momentum_callbacks'] or slot['preparation_callbacks']))):
            raise ValueError('native work accounting differs')
        keys(slot['density'], DENSITY_KEYS, 'native density')
        density = {key: finite(slot['density'][key], 'native density/' + key) for key in DENSITY_KEYS}
        if density['quadratic'] < 0:
            raise ValueError('native nonnegative quadratic required')
        with localcontext() as context:
            context.prec = 180
            error = abs(density['log_density'] + sum(density[key] for key in DENSITY_KEYS[:3]) / 2)
            scale = max(Decimal(1), sum(abs(density[key]) for key in DENSITY_KEYS))
            if error > Decimal.from_float(16 * sys.float_info.epsilon) * scale:
                raise ValueError('native normalized density identity differs')
        normalization = 13 * math.log(2 * math.pi)
        if abs(slot['density']['normalization'] - normalization) > 16 * sys.float_info.epsilon * normalization:
            raise ValueError('native full13 normalization differs')
        shared = (slot['density']['log_determinant'], slot['density']['normalization'])
        if shared_density is not None and shared != shared_density:
            raise ValueError('native shared full13 covariance density differs')
        shared_density = shared
    for name in WORK_KEYS:
        if count(native[name], request['producer_policy']['maximum_callbacks_total'], 'batch/' + name) != totals[name]:
            raise ValueError('native batch work totals differ')


def compare(native, reference, request, rows, reference_identity=None):
    check_native(native, request, rows)
    validate_reference(reference, request, rows, reference_identity)
    checks = []
    budgets = request['budgets']
    with localcontext() as context:
        context.prec = 180
        def scalar(name, value, coarse, fine, absolute, relative, outer_control, momentum_control, tail):
            actual = finite(value, name)
            low = decimal_string(coarse, name + '/coarse')
            high = decimal_string(fine, name + '/fine')
            budget = Decimal(str(absolute)) + Decimal(str(relative)) * abs(high)
            refinement = abs(low - high)
            outer_refinement = abs(decimal_string(outer_control, name + '/outer control') - high)
            momentum_refinement = abs(decimal_string(momentum_control, name + '/momentum control') - high)
            tail_bound = decimal_string(tail, name + '/tail bound')
            reference_envelope = outer_refinement + momentum_refinement + tail_bound
            difference = abs(actual - high)
            check = {'name': name, 'native_binary64_exact': str(actual),
                     'reference_coarse': str(low), 'reference_fine': str(high),
                     'difference': str(difference), 'budget': str(budget),
                     'fraction': str(difference / budget), 'reference_refinement': str(refinement),
                     'reference_fraction': str(refinement / budget),
                     'outer_refinement': str(outer_refinement), 'momentum_tail_refinement': str(momentum_refinement),
                     'analytic_tail_bound': str(tail_bound), 'reference_envelope': str(reference_envelope)}
            checks.append(check)
            if refinement > Decimal(str(budgets['reference_fraction'])) * budget:
                raise ComparisonFailure('reference refinement failed ' + name, checks)
            if (outer_refinement > Decimal(str(budgets['reference_fraction'])) * budget or
                    momentum_refinement > Decimal(str(budgets['reference_fraction'])) * budget or
                    reference_envelope > Decimal(str(budgets['reference_fraction'])) * budget):
                raise ComparisonFailure('independent reference axis/tail allocation failed ' + name, checks)
            if difference > budget:
                raise ComparisonFailure('frozen native comparison failed ' + name, checks)
        for index, (slot, coarse, fine) in enumerate(zip(native['slots'], reference['coarse'], reference['fine'])):
            outer_control = reference['momentum_refined'][index]
            momentum_control = reference['outer_refined'][index]
            for row, (value, low, high) in enumerate(zip(slot['predictions'], coarse['predictions'], fine['predictions'])):
                scalar(f'model{index}/row{row}', value, low, high,
                       budgets['ratio_absolute'], budgets['ratio_relative'],
                       outer_control['predictions'][row], momentum_control['predictions'][row],
                       fine['tail_bounds']['predictions'][row])
            for key in DENSITY_KEYS:
                scalar(f'model{index}/{key}', slot['density'][key], coarse['density'][key],
                       fine['density'][key], budgets['density_absolute'], 0,
                       outer_control['density'][key], momentum_control['density'][key],
                       fine['tail_bounds'].get(key, '0'))
    if len(checks) != 34:
        raise ValueError('complete full13/density comparison coverage required')
    return checks


def validate_reference(reference, request, rows, identity=None):
    keys(reference, ('schema_version', 'method', 'mapping', 'constants', 'mpmath_version',
                     'model_order', 'models', 'query_order', 'queries', 'sources', 'request_sha256',
                     'request_canonical_sha256', 'reference_input_sha256', 'reference_script_sha256',
                     'settings', 'runtime_before', 'runtime_after', *REFERENCE_SETTINGS), 'reference')
    for name, wanted in (('schema_version', 1), ('mpmath_version', '1.3.0'),
                         ('method', REFERENCE_METHOD), ('mapping', MAPPING), ('constants', REFERENCE_CONSTANTS),
                         ('model_order', MODEL_ORDER), ('models', request['models']),
                         ('query_order', [row['id'] for row in queries(rows)]),
                         ('queries', rows), ('sources', INPUT_IDENTITIES),
                         ('request_sha256', REQUEST_SHA256), ('settings', REFERENCE_SETTINGS),
                         ('request_canonical_sha256', hashlib.sha256(json.dumps(request, sort_keys=True,
                                                                  separators=(',', ':')).encode()).hexdigest())):
        exact(reference[name], wanted, 'reference/' + name)
    for name in ('reference_input_sha256', 'reference_script_sha256'):
        if not isinstance(reference[name], str) or not re.fullmatch('[0-9a-f]{64}', reference[name]):
            raise ValueError('reference digest domain differs')
    if identity is not None:
        for name in ('reference_input_sha256', 'reference_script_sha256'):
            exact(reference[name], identity[name], 'reference/' + name)
        exact(reference['runtime_before'], identity['runtime'], 'reference runtime before')
    exact(reference['runtime_after'], reference['runtime_before'], 'reference runtime after')
    verify_runtime(reference['runtime_before'], Path(reference['runtime_before']['python']['executable']['path']))
    for key, expected in REFERENCE_SETTINGS.items():
        values = reference[key]
        if not isinstance(values, list) or len(values) != 2:
            raise ValueError('complete two-model reference required')
        for index, value in enumerate(values):
            keys(value, ('model', 'ruler_mpc', 'predictions', 'density', 'tail_bounds', *expected), 'reference model')
            exact(value['model'], MODEL_ORDER[index], 'reference model order')
            for name, wanted in expected.items():
                exact(value[name], wanted, 'reference model settings/' + name)
            if decimal_string(value['ruler_mpc'], 'reference ruler') <= 0:
                raise ValueError('positive reference ruler required')
            if not isinstance(value['predictions'], list) or len(value['predictions']) != 13:
                raise ValueError('reference full13 predictions required')
            for prediction in value['predictions']:
                if decimal_string(prediction, 'reference prediction') <= 0:
                    raise ValueError('positive reference prediction required')
            keys(value['density'], DENSITY_KEYS, 'reference density')
            density = {name: decimal_string(value['density'][name], 'reference density/' + name) for name in DENSITY_KEYS}
            with localcontext() as context:
                context.prec = 180
                if density['quadratic'] < 0:
                    raise ValueError('reference nonnegative quadratic required')
                error = abs(density['log_density'] + sum(density[name] for name in DENSITY_KEYS[:3]) / 2)
                scale = max(Decimal(1), sum(abs(number) for number in density.values()))
                if error > Decimal(10) ** (6 - expected['digits']) * scale:
                    raise ValueError('reference normalized density identity differs')
            bounds = value['tail_bounds']
            keys(bounds, ('scaled_expansion', 'relative_scaled_expansion', 'predictions', 'quadratic', 'log_density'),
                 'reference analytic tail bounds')
            if not isinstance(bounds['predictions'], list) or len(bounds['predictions']) != 13:
                raise ValueError('reference full13 analytic tail bounds required')
            for number in [bounds[name] for name in ('scaled_expansion', 'relative_scaled_expansion', 'quadratic', 'log_density')] + bounds['predictions']:
                if decimal_string(number, 'reference tail bound') < 0:
                    raise ValueError('nonnegative reference tail bounds required')
            if index == 0 and any(decimal_string(number, 'massless tail bound') != 0
                                  for number in bounds['predictions'] + [bounds[name] for name in
                                                                        ('scaled_expansion', 'relative_scaled_expansion', 'quadratic', 'log_density')]):
                raise ValueError('massless reference must have no momentum tail')


def child(command, store, label, limits, stdin=None, record=None):
    def bound():
        resource.setrlimit(resource.RLIMIT_AS, (limits['memory_bytes'], limits['memory_bytes']))
        resource.setrlimit(resource.RLIMIT_FSIZE, (limits['output_bytes'], limits['output_bytes']))
        resource.setrlimit(resource.RLIMIT_CPU, (math.ceil(limits[label + '_seconds']) + 1,) * 2)
    operation = {'command': command, 'status': 'started', 'returncode': None, 'timed_out': False}
    if record is not None:
        record.setdefault('subprocesses', {})[label] = operation
    try:
        with (store / (label + '.out')).open('xb') as out, (store / (label + '.err')).open('xb') as err:
            result = subprocess.run(command, cwd=ROOT, stdout=out, stderr=err, input=stdin,
                                    text=True, timeout=limits[label + '_seconds'], preexec_fn=bound,
                                    env={**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1',
                                         'MKL_NUM_THREADS': '1', 'NUMEXPR_NUM_THREADS': '1',
                                         'VECLIB_MAXIMUM_THREADS': '1', 'PYTHONHASHSEED': '0',
                                         'PYTHONPATH': '', 'PYTHONNOUSERSITE': '1'})
        operation['returncode'] = result.returncode
        if result.returncode != 0:
            raise ValueError(f'{label} failed with status {result.returncode}')
        operation['status'] = 'completed'
        return (store / (label + '.out')).read_text()
    except subprocess.TimeoutExpired:
        operation['timed_out'] = True
        operation['status'] = 'timeout'
        raise ValueError(label + ' timeout') from None
    except Exception:
        operation['status'] = 'failed'
        raise
    finally:
        for suffix in ('out', 'err'):
            path = store / (label + '.' + suffix)
            if path.exists():
                operation[suffix + '_bytes'] = path.stat().st_size
                operation[suffix + '_sha256'] = sha256(path)


def ingest_outputs(store, record):
    for label in ('native', 'reference', 'runtime'):
        path = store / (label + '.out')
        if path.is_file():
            try:
                record[label if label != 'runtime' else 'reference_runtime'] = parse(path.read_text())
            except Exception as error:
                record.setdefault('partial_output_errors', {})[label] = str(error)


def execute(args):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,95}', args.name):
        raise ValueError('attempt name must be one bounded path component')
    store = within(ROOT, 'results/lcdm-campaign/' + args.name)
    store.mkdir(parents=True, exist_ok=False)
    record = {'schema_version': 1, 'status': 'started', 'started_utc': datetime.now(timezone.utc).isoformat(),
              'gates': {'execution': 'unassessed', 'numerical': 'unassessed', 'inference': 'not performed',
                        'interpretation': 'conditional controls; full model blocked'},
              'controller_python': sys.version}
    started = time.monotonic()
    usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    before = None
    blobs = None
    packet_hashes = None
    source_before = None
    runtime_before = None
    executable_hash = None
    snapshot_hashes = None
    artifacts = getattr(args, 'engine_artifacts', None) or args.engine_source
    try:
        source_before = source_identity()
        record['source_before'] = source_before
        packet, _, _ = read_packet(FOLDER)
        validate_admission(packet)
        request = validate_request_bytes((FOLDER / 'request.json').read_bytes())
        record['request'] = request
        record['request_sha256'] = REQUEST_SHA256
        packet_hashes = {path.name: sha256(path) for path in FOLDER.iterdir() if path.is_file()}
        record['packet_hashes'] = packet_hashes
        snapshot = snapshot_source(store, source_before, record)
        for name in packet_hashes:
            (store / name).write_bytes((snapshot / FOLDER.relative_to(ROOT) / name).read_bytes())
            if sha256(store / name) != packet_hashes[name]:
                raise ValueError('copied immutable packet identity differs')
        snapshot_hashes = {str(path.relative_to(store)): sha256(path) for path in store.rglob('*')
                           if path.is_file() and (path.is_relative_to(snapshot) or path.name in packet_hashes or path.name == 'source.tar')}
        record['snapshot_hashes'] = snapshot_hashes
        blobs, rows, covariance = inputs(args.input_root, packet)
        record['input_sources'] = input_sources(packet)
        record['input_hashes'] = [hashlib.sha256(blob).hexdigest() for blob in blobs]
        for index, blob in enumerate(blobs):
            (store / f'original-input-{index}.txt').write_bytes(blob)
            snapshot_hashes[f'original-input-{index}.txt'] = hashlib.sha256(blob).hexdigest()
        before = fingerprint(args.engine_source, request['engine'], artifacts)
        record['engine_before'] = before
        limits = {**request['resources'], 'discovery_seconds': 30, 'runtime_seconds': 30}
        discovery = child([str(artifacts / 'target/release/irred'), 'describe', '--json'],
                          store, 'discovery', limits, record=record)
        if len(discovery.encode()) > 1048576:
            raise ValueError('discovery byte quota')
        (store / 'discovery.json').write_text(discovery)
        record['discovery_sha256'] = sha256(store / 'discovery.json')
        description = parse(discovery)
        if (description.get('product') != 'Irreducible' or
                description['build']['build_id'] != request['engine']['build_id'] or
                description['build']['git_head'] != request['engine']['revision'] or
                description['build']['git_status'] != ''):
            raise ValueError('compiled discovery identity differs')
        record['local_operation'] = {'id': request['operation'], 'interface': 'experiment-specific native consumer'}
        command = ['/usr/bin/c++', '-std=c++20', '-O2', '-Wall', '-Wextra', '-Wpedantic',
                   '-fno-fast-math', '-ffp-contract=off', str(store / 'consumer.cpp'), '-I',
                   str(args.engine_source / 'cpp/include'), str(artifacts / 'build/native-release/libirred_core.a'),
                   '-o', str(store / 'consumer')]
        record['compile_command'] = command
        child(command, store, 'compile', limits, record=record)
        executable_hash = sha256(store / 'consumer')
        record['consumer_sha256'] = executable_hash
        payload = transport(request, rows, covariance)
        (store / 'transport.txt').write_text(payload)
        snapshot_hashes['transport.txt'] = sha256(store / 'transport.txt')
        native = parse(child([str(store / 'consumer')], store, 'native', limits, payload, record))
        record['native'] = native
        check_native(native, request, rows)
        record['gates']['execution'] = 'passed'
        reference_input = {'schema_version': 1, 'request': request, 'request_sha256': REQUEST_SHA256,
                           'rows': rows, 'covariance': covariance, 'sources': input_sources(packet)}
        (store / 'reference-input.json').write_text(json.dumps(reference_input, allow_nan=False) + '\n')
        snapshot_hashes['reference-input.json'] = sha256(store / 'reference-input.json')
        record['reference_python'] = str(args.reference_python)
        record['reference_python_sha256'] = sha256(args.reference_python)
        runtime_before = parse(child([str(args.reference_python), str(store / 'reference.py'), '--runtime-fingerprint'],
                                     store, 'runtime', limits, record=record))
        verify_runtime(runtime_before, args.reference_python)
        record['reference_runtime'] = runtime_before
        reference_identity = {'runtime': runtime_before,
                              'reference_input_sha256': sha256(store / 'reference-input.json'),
                              'reference_script_sha256': sha256(store / 'reference.py')}
        reference = parse(child([str(args.reference_python), str(store / 'reference.py'),
                                 str(store / 'reference-input.json')], store, 'reference', limits, record=record))
        record['reference'] = reference
        record['comparisons'] = compare(native, reference, request, rows, reference_identity)
        record['gates']['numerical'] = 'passed named thermal-DESI direct-momentum/quadrature/LDLT comparisons'
        record['status'] = 'completed'
    except Exception as error:
        record['status'] = 'failed'
        record['error'] = str(error)
        if isinstance(error, ComparisonFailure):
            record['comparisons'] = error.checks
        if 'native' in record and record['gates']['execution'] != 'passed':
            record['gates']['execution'] = 'failed native admission'
    finally:
        ingest_outputs(store, record)
        errors = []
        if source_before is not None:
            try:
                record['source_after'] = source_identity()
                if record['source_after'] != source_before:
                    errors.append('committed Reproducible source changed')
            except Exception as error:
                errors.append('source final verification: ' + str(error))
        if before is not None:
            try:
                record['engine_after'] = fingerprint(args.engine_source, record['request']['engine'], artifacts)
                if record['engine_after'] != before:
                    errors.append('engine changed')
            except Exception as error:
                errors.append('engine final verification: ' + str(error))
        if runtime_before is not None:
            try:
                verify_runtime(runtime_before, args.reference_python)
            except Exception as error:
                errors.append('reference runtime final verification: ' + str(error))
        if blobs is not None:
            try:
                after, _, _ = inputs(args.input_root, packet)
                if after != blobs:
                    errors.append('inputs changed')
            except Exception as error:
                errors.append('input final verification: ' + str(error))
        if packet_hashes is not None:
            try:
                if {path.name: sha256(path) for path in FOLDER.iterdir() if path.is_file()} != packet_hashes:
                    errors.append('packet changed')
            except Exception as error:
                errors.append('packet final verification: ' + str(error))
        if snapshot_hashes is not None:
            try:
                if any(sha256(within(store, name)) != digest for name, digest in snapshot_hashes.items()):
                    errors.append('immutable committed snapshot changed')
            except Exception as error:
                errors.append('snapshot final verification: ' + str(error))
        if executable_hash is not None:
            try:
                if sha256(store / 'consumer') != executable_hash:
                    errors.append('compiled consumer changed')
            except Exception as error:
                errors.append('compiled consumer final verification: ' + str(error))
        if errors:
            record['status'] = 'failed'
            record.setdefault('error', 'post-run identity failure')
            record['gates']['numerical'] = 'not accepted: identity failure'
        record['integrity_errors'] = errors
        record['elapsed_seconds'] = time.monotonic() - started
        record['finished_utc'] = datetime.now(timezone.utc).isoformat()
        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        record['children_resources'] = {'cpu_user_seconds': usage.ru_utime - usage_before.ru_utime,
                                        'cpu_system_seconds': usage.ru_stime - usage_before.ru_stime,
                                        'process_lifetime_maximum_rss_kib': usage.ru_maxrss, 'jobs': 1, 'threads': 1}
        record['output_sha256'] = {str(path.relative_to(store)): sha256(path) for path in store.rglob('*') if path.is_file()}
        (store / 'record.json').write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')
        for path in store.rglob('*'):
            if path.is_file():
                path.chmod(0o555 if path.name == 'consumer' else 0o444)
    print(json.dumps({'status': record['status'], 'record': str(store / 'record.json'),
                      'gates': record['gates'], 'error': record.get('error')}))
    return 0 if record['status'] == 'completed' else 1


def verify_runtime(value, executable):
    keys(value, ('python', 'mpmath'), 'reference runtime')
    python = value['python']
    keys(python, ('version', 'implementation', 'executable'), 'reference Python')
    if not isinstance(python['version'], str) or not python['version'] or python['implementation'] != 'CPython':
        raise ValueError('reference Python version/implementation differs')
    path = Path(executable).absolute()
    expected = {'path': str(path), 'resolved_path': str(path.resolve()),
                'bytes': path.stat().st_size, 'sha256': sha256(path)}
    exact(python['executable'], expected, 'reference Python executable')
    module = value['mpmath']
    keys(module, ('version', 'backend', 'module_root', 'source_inventory', 'source_inventory_sha256'), 'reference mpmath')
    if module['version'] != '1.3.0' or module['backend'] != 'python':
        raise ValueError('reviewed pure-Python mpmath1.3.0 required')
    root = Path(module['module_root'])
    if not root.is_absolute() or str(root.resolve()) != str(root) or not root.is_dir():
        raise ValueError('reference module root differs')
    inventory = module['source_inventory']
    if not isinstance(inventory, list) or not 1 <= len(inventory) <= 256:
        raise ValueError('reference module inventory quota differs')
    paths = sorted(root.rglob('*.py'))
    actual = []
    total = 0
    for source in paths:
        if not source.is_file() or source.is_symlink():
            raise ValueError('reference module source type differs')
        size = source.stat().st_size
        total += size
        if total > 8388608:
            raise ValueError('reference module byte quota differs')
        actual.append({'path': str(source), 'relative_path': str(source.relative_to(root)),
                       'bytes': size, 'sha256': sha256(source)})
    exact(inventory, actual, 'reference complete module inventory')
    digest = hashlib.sha256(json.dumps(actual, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    exact(module['source_inventory_sha256'], digest, 'reference module inventory digest')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('engine-source', 'input-root', 'reference-python'):
        parser.add_argument('--' + key, type=Path, required=True)
    parser.add_argument('--engine-artifacts', type=Path)
    parser.add_argument('--name', required=True)
    args = parser.parse_args()
    args.engine_source = args.engine_source.resolve()
    args.engine_artifacts = args.engine_artifacts.resolve() if args.engine_artifacts else None
    args.input_root = args.input_root.resolve()
    args.reference_python = args.reference_python.absolute()
    sys.exit(execute(args))
