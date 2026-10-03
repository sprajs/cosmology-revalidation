"""Separate retained-CLASS DESI probe; new configured successor to attempt C.

Source preparation only until reviewed source/runtime bindings and a job grant.
No CLASS rerun, Gaussian implementation, fitting or combined-probe target.
"""
import hashlib
import json
import locale
import os
from pathlib import Path
import resource
import re
import signal
import stat
import sys
import time
import types

ROOT = Path(__file__).absolute().parents[2]
PROPOSAL = Path(__file__).absolute().parent
SDK = ATTEMPT = None
REQUEST_SHA = 'db50c6679a50be51bf2a80234c9cc893b1dd792307a4c8bff586f8fade29a95f'
BUILD = '6d495efb166006c6ce651359a366af5d87686ce516ecd7f0f49eed4e8d8be07e'
REV = 'f3b19b6539c72fda5a8b489a0546ca4d0266ace3'
MODULES = [
    ('theory', ROOT / 'experiments/lcdm-reference/theory.py', 'bb737cd123d72997ede5f62bc461d5c8f56fbc43aaf8f706ec8b83f6fd8ab273'),
    ('run', ROOT / 'experiments/lcdm-reference/run.py', 'e764e828e495c66eda296a3404759350fdcb0860e0ea2f817818c4cf7d35ec44'),
    ('bao', PROPOSAL / 'bao.py', '5a09d83648ca547a6d785ae234fb2e6a7edc747550a588d6395d4956af3d00bf'),
    ('transport', PROPOSAL / 'transport.py', 'd254c758d4c18248ba3802ce33916f5f2bd8450c39e3102a7116f55fd1437b05'),
    ('prepare_probe', PROPOSAL / 'prepare_probe.py', 'd084193d7568cb465e42573e39430acd09941e3fbc2511541baede727e3c7a31')]
TOOLS = [
    ('/usr/bin/bash', 'c14c8f498f3be96f076ceb07d65dba584a86dca3c86b5dce65ed0f6cae73358d'),
    ('/usr/bin/c++', 'f04191f6a7b2cd7d9a62e1745872b8a6088791e5af6955488c69c9b2c4668bc9'),
    ('/usr/lib/gcc/x86_64-pc-linux-gnu/16/libstdc++.so', 'f5fc7380f2ae46fa4053a64be04e7b98109f1066a4bbfff3c37042488aa0be0e'),
    ('/usr/lib/libm.so.6', '965106704753eefe8c31ae7da6daba3ead3d18c52b19b5fcc41178bcc4aa999d'),
    ('/usr/lib/libc.so.6', '02c8c3d06180beee9c7d0392aa85736c7d19f274a81e7f8f0b3d9e28e0641524'),
    ('/usr/lib/ld-linux-x86-64.so.2', 'f5e11cc62f8f2c24dff532982559f4303dba3b17e7450551372a88c8e7ea8757')]
LIMITS = {'address_bytes': 2147483648, 'case_cpu_seconds': 120,
          'case_wall_seconds': 120, 'file_bytes': 16777216, 'max_files': 32,
          'attempt_bytes': 33554432, 'combined_logs_bytes': 1048576}


def require(value, message):
    if not value:
        raise ValueError(message)


def read(path, expected=None, limit=33554432):
    path = Path(path).absolute()
    require('..' not in path.parts, 'path traversal')
    for parent in path.parents:
        require(parent.is_dir() and not parent.is_symlink(), 'source ancestor')
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    try:
        before = os.fstat(fd)
        require(stat.S_ISREG(before.st_mode) and before.st_size <= limit, 'source size/type')
        chunks, count = [], 0
        while True:
            block = os.read(fd, min(65536, limit + 1 - count))
            if not block:
                break
            chunks.append(block)
            count += len(block)
            require(count <= limit, 'consumed source size')
        after, linked = os.fstat(fd), path.lstat()
        facts = lambda x: (x.st_dev, x.st_ino, x.st_size, x.st_mode, x.st_mtime_ns, x.st_ctime_ns)
        require(facts(before) == facts(after) == facts(linked), 'source changed while consumed')
        raw = b''.join(chunks)
        pin = {'path': str(path), 'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest(),
               'device': after.st_dev, 'inode': after.st_ino, 'mode': oct(stat.S_IMODE(after.st_mode)),
               'mtime_ns': after.st_mtime_ns, 'ctime_ns': after.st_ctime_ns}
        if expected is not None:
            require(all(pin[k] == expected[k] for k in ('path', 'bytes', 'sha256')), 'source pin')
        return raw, pin
    finally:
        os.close(fd)


def write(path, value):
    raw = (json.dumps(value, sort_keys=True, allow_nan=False, indent=2) + '\n').encode()
    require(len(raw) <= 1048576, 'receipt/product byte cap')
    with Path(path).open('xb') as f:
        require(f.write(raw) == len(raw), 'write count')
        f.flush()
        os.fsync(f.fileno())
    Path(path).chmod(0o444)


def positive_pin(value):
    require(type(value) is dict and set(value) == {'path', 'bytes', 'sha256'}, 'closed file authority')
    require(type(value['path']) is str and len(value['path'].encode()) <= 4096 and
            Path(value['path']).is_absolute() and '..' not in Path(value['path']).parts, 'authority path')
    require(type(value['bytes']) is int and 0 <= value['bytes'] <= 33554432 and
            type(value['sha256']) is str and re.fullmatch(r'[0-9a-f]{64}', value['sha256']),
            'positive typed byte authority')
    return value


def pairs(items):
    value = {}
    for key, item in items:
        require(key not in value, 'duplicate configuration key')
        value[key] = item
    return value


def decode(raw):
    return json.loads(raw, object_pairs_hook=pairs,
                      parse_constant=lambda _: (_ for _ in ()).throw(ValueError('nonfinite JSON')))


def module(name, path, digest):
    raw, pin = read(path, limit=262144)
    require(pin['sha256'] == digest, 'RAM source hash ' + name)
    require(name not in sys.modules, 'unexpected project module ' + name)
    value = types.ModuleType(name)
    value.__file__ = str(path)
    sys.modules[name] = value
    exec(compile(raw, str(path), 'exec'), value.__dict__)
    return value


def runtime():
    modules, files = [], set()
    for name, value in sorted(sys.modules.items()):
        origin = getattr(value, '__file__', None)
        cached = getattr(value, '__cached__', None)
        modules.append({'name': name, 'origin': origin, 'cached': cached})
        for path in (origin, cached):
            if path and Path(path).is_file():
                files.add(str(Path(path).resolve()))
    with open('/proc/self/maps', 'rb') as f:
        raw = f.read(1048577)
    require(len(raw) <= 1048576, 'runtime map bound')
    mapped = []
    for line in raw.decode('ascii').splitlines():
        fields = line.split(None, 5)
        if len(fields) == 6 and fields[5].startswith('/'):
            require(not fields[5].endswith(' (deleted)'), 'deleted mapped library')
            path = fields[5]
            if '.so' in Path(path).name:
                mapped.append(path)
                files.add(str(Path(path).resolve()))
    require(len(files) <= 256 and len(modules) <= 256, 'runtime file/module cap')
    inventory = [read(path)[1] for path in sorted(files)]
    require(sum(pin['bytes'] for pin in inventory) <= 67108864, 'runtime file-byte cap')
    with open('/proc/self/status', 'rb') as f:
        status = f.read(65537)
    require(len(status) <= 65536 and b'\nThreads:\t1\n' in status, 'actual runtime thread count')
    require(not any(name.split('.')[0] in ('numpy', 'scipy', 'mpmath') for name in sys.modules), 'unexpected numerical runtime')
    return {'version': sys.version, 'executable': sys.executable, 'modules': modules,
            'mapped_libraries': sorted(set(mapped)), 'inventory': inventory, 'threads': 1}


def main():
    global SDK, ATTEMPT
    start = time.monotonic()
    whole_end, phase_end = start + 300, start + 120
    def expired(_signum, _frame):
        signal.signal(signal.SIGALRM, signal.SIG_DFL)
        signal.setitimer(signal.ITIMER_REAL, max(0.001, phase_end - time.monotonic()))
        raise TimeoutError('kernel timer: phase allowance exhausted; terminal reserve remains')
    def arm(end):
        nonlocal phase_end
        phase_end = end
        signal.signal(signal.SIGALRM, expired)
        signal.setitimer(signal.ITIMER_REAL, max(0.001, end - 10 - time.monotonic()))
    def deadline():
        require(time.monotonic() < phase_end - 10, 'phase preparation/native deadline')
    resource.setrlimit(resource.RLIMIT_AS, (2147483648,) * 2)
    resource.setrlimit(resource.RLIMIT_CPU, (300,) * 2)
    resource.setrlimit(resource.RLIMIT_FSIZE, (16777216,) * 2)
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    def interrupted(_signum, _frame):
        raise TimeoutError('outer interrupted; reviewed child capture performs cleanup')
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    arm(phase_end)
    require(len(sys.argv) == 5 and sys.argv[1] == '--config' and sys.argv[3] == '--attempt',
            'usage: controller.py --config ABSOLUTE_CONFIG --attempt FRESH_LABEL')
    require(len(sys.argv[2].encode()) <= 4096 and Path(sys.argv[2]).is_absolute(), 'config path')
    require(len(sys.argv[4]) <= 80 and re.fullmatch(r'[a-z0-9]+(?:-[a-z0-9]+)*', sys.argv[4]),
            'fresh attempt label')
    require(all(p.is_dir() and not p.is_symlink() for p in (ROOT, *ROOT.parents)), 'repository path')
    for folder in (ROOT / 'results', ROOT / 'results/lcdm-bao-reference'):
        if not folder.exists():
            folder.mkdir(mode=0o700)
        require(folder.is_dir() and not folder.is_symlink(), 'result directory')
    ATTEMPT = ROOT / 'results/lcdm-bao-reference' / sys.argv[4]
    ATTEMPT.mkdir(mode=0o700)
    record = {'schema': 'retained-class-separate-DESI-native-probe/v1',
              'source_route': 'lcdm-bao-reference-configured-successor/v1',
              'historical_C_record_sha256': '97dd8c95ddd5021ef387e34628b24e8b2580d6fae788f8cee91da967fef31b33', 'status': 'started',
              'sdk_revision': REV, 'sdk_build_id': BUILD, 'whole_cap_seconds': 300,
              'children': [{'label': n, 'status': 'not_started'} for n in ('describe', 'compile', 'native')],
              'before': [], 'after': [], 'prepared': None, 'native': None, 'failures': [],
              'runtime': {'version': sys.version, 'executable': sys.executable,
                          'stdlib_only': True, 'native_dynamic_maps_observed': False},
              'gates': {'execution': 'failed', 'native_numeric': 'unearned',
                        'likelihood_accuracy': None, 'joint_inference': 'blocked', 'SN_score': None},
              'descendant_cleanup_proven': False, 'aggregate_descendant_CPU_enforced': False}
    write(ATTEMPT / 'started.json', record)
    runner = None
    expected = []
    record['before'] = expected
    try:
        require(sys.dont_write_bytecode and os.environ.get('PYTHONDONTWRITEBYTECODE') == '1' and
                all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                    'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')) and
                not any(k in os.environ for k in ('PYTHONPATH', 'PYTHONHOME', 'LD_PRELOAD', 'LD_LIBRARY_PATH')),
                'reviewed bytecode/thread/no-loader-alias environment')
        raw, config_pin = read(sys.argv[2], limit=65536)
        expected.append(config_pin)
        record['configuration'] = config_pin
        config = decode(raw)
        require(type(config) is dict and set(config) == {'schema', 'controller', 'source_pins',
                'sdk_directory', 'native_archive'} and
                config['schema'] == 'lcdm-bao-reference-local-runtime/v1', 'closed local runtime configuration')
        controller_pin = positive_pin(config['controller'])
        require(controller_pin['path'] == str(Path(__file__).absolute()), 'actual controller source path')
        _, actual_controller = read(controller_pin['path'], controller_pin, limit=262144)
        expected.append(actual_controller)
        require(type(config['sdk_directory']) is str and Path(config['sdk_directory']).is_absolute() and
                '..' not in Path(config['sdk_directory']).parts and len(config['sdk_directory'].encode()) <= 4096,
                'read-only SDK directory path')
        SDK = Path(config['sdk_directory'])
        raw, request_pin = read(PROPOSAL / 'request.json', limit=65536)
        require(request_pin['sha256'] == REQUEST_SHA, 'frozen BAO design request')
        expected.append(request_pin)
        record['request'] = request_pin
        request = decode(raw)
        python_path = Path(sys.executable).resolve()
        raw, py = read(python_path)
        require(py['sha256'] == '8f9781a98200d9ecda7e00464e4c64b1327abae788ae8e6979d5c859311410c7'
                and sys.version_info[:2] == (3, 12), 'actual reviewed Python')
        expected.append(py)
        for name, path, digest in MODULES:
            _, pin = read(path, limit=262144)
            require(pin['sha256'] == digest, 'all project source admission before import')
            expected.append(pin)
        for path, digest in TOOLS:
            actual = Path(path).resolve()
            _, pin = read(actual)
            require(pin['sha256'] == digest, 'actual tool/dependency pin')
            expected.append(pin)
        source_authority = positive_pin(config['source_pins'])
        raw, source_pins = read(source_authority['path'], source_authority, limit=65536)
        expected.append(source_pins)
        input_sources = decode(raw)
        require(type(input_sources) is dict and set(input_sources) == {'schema_version', 'cases', 'scope', 'source_pins'}
                and type(input_sources['schema_version']) is int and input_sources['schema_version'] == 1,
                'closed retained source-pins object')
        require(type(input_sources['source_pins']) is list and len(input_sources['source_pins']) == 26,
                'fixed retained source-role count')
        require([{k: p[k] for k in ('role', 'bytes', 'sha256')} for p in input_sources['source_pins']] ==
                request['retained_source_roles'], 'unchanged retained source roles/byte identities')
        for pin in input_sources['source_pins']:
            require(type(pin) is dict and set(pin) == {'role', 'path', 'bytes', 'sha256'}, 'source-role authority')
            positive_pin({k: pin[k] for k in ('path', 'bytes', 'sha256')})
            _, consumed = read(pin['path'], pin)
            expected.append(consumed)
        raw, admission = read(SDK / 'sdk-admission.json', limit=1048576)
        require(admission['sha256'] == '406675b96e741b8bd1aee6d1b9008c581d47f413442e945357334c38d9dc3ba1', 'complete SDK source admission')
        expected.append(admission)
        sdk = decode(raw)
        require(len(sdk['sdk_header_copies']) == 46 and
                all(Path(p['path']).is_relative_to(SDK / 'sdk/include') for p in sdk['sdk_header_copies']),
                'actual consumed include path/header inventory')
        record['documentary_log_ancestry_unconsumed'] = sdk['artifacts_and_tools'][5:]
        for pin in sdk['sdk_header_copies'] + sdk['artifacts_and_tools'][:5]:
            actual = Path(pin['path']).resolve()
            _, consumed = read(actual, {**pin, 'path': str(actual)})
            expected.append(consumed)
        archive_authority = positive_pin(config['native_archive'])
        raw, archive = read(archive_authority['path'], archive_authority)
        require(archive['sha256'] == 'be08158b14c77f7692fe2a55caf9aaac4abe97eef1edd75e4966e57f27516d77', 'copied native archive')
        expected.append(archive)
        raw, consumer = read(PROPOSAL / 'consumer.cpp', limit=65536)
        require(consumer['sha256'] == 'b8047f44a9de2b23e02d98fb4ea1d622ed511cc9d5f24f819a134fe60dcbabc4', 'native caller source')
        expected.append(consumer)
        raw, manifest_pin = read(SDK / 'build-manifest-release.json', limit=1048576)
        require(manifest_pin['sha256'] == '77317e5a47fa9fa20c29f6cf7f4b8dfd5ce7262cf39ebcb1d1637af03d70cd85', 'retained manifest')
        expected.append(manifest_pin)
        manifest = decode(raw)
        record['before'] = expected
        loaded = {name: module(name, path, sha) for name, path, sha in MODULES}
        runner, preparer = loaded['run'], loaded['prepare_probe']
        record['runtime_before'] = runtime()
        def child(label, argv, seconds, cap):
            deadline()
            folder = ATTEMPT / label
            folder.mkdir(mode=0o700)
            index = ('describe', 'compile', 'native').index(label)
            record['children'][index] = {'label': label, 'status': 'entering-capture', 'argv': argv}
            write(ATTEMPT / (label + '.entering.json'), record['children'][index])
            result = runner.run_class(argv, str(ROOT), folder, ATTEMPT,
                                      {**LIMITS, 'case_cpu_seconds': seconds, 'case_wall_seconds': seconds},
                                      phase_end - 7, cap)
            record['children'][index] = {'label': label, **result}
            for field in ('stdout_identity', 'stderr_identity'):
                if field in result:
                    expected.append(result[field])
            require(result['status'] == 'completed', label + ' failed; retained raw prefix')
            return folder
        folder = child('describe', ['/home/szymon/Projects/irreducible/target/release/irred', 'describe', '--json'], 30, 1048576)
        raw, _ = read(folder / 'stdout.log', record['children'][0]['stdout_identity'], limit=1048576)
        describe = json.loads(raw, object_pairs_hook=runner.pairs)
        require(describe['build'] == manifest
                and describe['build']['build_id'] == BUILD and describe['build']['git_head'] == REV,
                'compiled discovery/retained full manifest identity')
        record['discovery'] = {'build_id': describe['build']['build_id'], 'git_head': describe['build']['git_head'],
                               'raw': record['children'][0]['stdout_identity']}
        argv = ['/usr/bin/c++', '-std=c++20', '-O2', '-Wall', '-Wextra', '-Wpedantic', '-fno-fast-math',
                '-ffp-contract=off', '-DIRRED_CLASS_BAO_BUILD_ID="' + BUILD + '"', str(PROPOSAL / 'consumer.cpp'),
                '-I', str(SDK / 'sdk/include'), archive['path'], '-o', str(ATTEMPT / 'consumer')]
        child('compile', argv, 120, 1048576)
        _, binary = read(ATTEMPT / 'consumer', limit=16777216)
        expected.append(binary)
        record['native_binary'] = binary
        arm(min(whole_end, time.monotonic() + 180))
        prepared = preparer.prepare(ATTEMPT, source_pins, deadline)
        record['prepared'] = {k: v for k, v in prepared.items() if k != 'dataset'}
        _, stdin = read(ATTEMPT / 'native.input', prepared['native_input'], limit=16384)
        expected.append(stdin)
        argv = ['/usr/bin/bash', '--noprofile', '--norc', '-c', 'exec "$1" <"$2"', '--',
                str(ATTEMPT / 'consumer'), str(ATTEMPT / 'native.input')]
        folder = child('native', argv, 180, 32768)
        raw, _ = read(folder / 'stdout.log', record['children'][2]['stdout_identity'], limit=32768)
        record['native'] = preparer.terminal_admission(raw, prepared)
        require(record['native']['status'] == 'accepted', 'native batch numeric refusal')
        record['gates'].update(execution='passed', native_numeric='finite-existing-empirical-solve-screen')
        record['status'] = 'completed'
    except BaseException as exc:
        record['failures'].append({'kind': type(exc).__name__, 'message': str(exc)[:2048],
                                   'error': getattr(exc, 'error', None)})
        if hasattr(exc, 'earned_prefix'):
            record['earned_failure_prefix'] = exc.earned_prefix
            if 'native_output' in exc.earned_prefix:
                record['native'] = exc.earned_prefix['native_output']
    finally:
        if runner is not None:
            try:
                record['runtime_after'] = runtime()
                require(record['runtime_after'] == record['runtime_before'], 'runtime before/after identity')
            except BaseException as exc:
                record['failures'].append({'stage': 'terminal-runtime', 'kind': type(exc).__name__, 'message': str(exc)[:512]})
        # Independent after checks preserve accepted objects and revoke disposition on drift.
        for pin in expected:
            try:
                _, after = read(pin['path'], pin)
                require(all(after[k] == pin[k] for k in ('device', 'inode', 'mode', 'mtime_ns', 'ctime_ns') if k in pin), 'terminal file drift')
                record['after'].append(after)
            except BaseException as exc:
                record['failures'].append({'stage': 'terminal-identity', 'path': pin['path'], 'kind': type(exc).__name__, 'message': str(exc)[:512]})
        if record['failures']:
            record['status'] = 'failed'
            record['gates'].update(execution='failed', native_numeric='withheld-on-failure')
        record['wall_seconds'] = time.monotonic() - start
        if runner is not None:
            record['sealed_outputs'] = runner.seal_inventory(ATTEMPT, LIMITS)
            if record['sealed_outputs']['errors']:
                record['status'] = 'failed'
                record['gates'].update(execution='failed', native_numeric='withheld-on-seal-failure')
        require(time.monotonic() < phase_end <= whole_end, 'whole/phase terminal cap')
        write(ATTEMPT / 'record.json', record)
        signal.setitimer(signal.ITIMER_REAL, 0)
        print(json.dumps({'status': record['status'], 'record': str(ATTEMPT / 'record.json'),
                          'wall_seconds': record['wall_seconds']}), flush=True)
    return 0 if record['status'] == 'completed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
