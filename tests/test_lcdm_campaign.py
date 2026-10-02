"""Adversarial admission and receipt checks; these are not scientific validation."""
import copy
from decimal import Decimal, localcontext
import importlib.util
import io
import json
import math
from pathlib import Path
import sys
import tarfile
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

FOLDER = Path(__file__).resolve().parents[1] / 'experiments/lcdm-campaign'
spec = importlib.util.spec_from_file_location('thermal_campaign_controller', FOLDER / 'controller.py')
controller = importlib.util.module_from_spec(spec)
spec.loader.exec_module(controller)


class CampaignTests(unittest.TestCase):
    def setUp(self):
        self.request = controller.validate_request_bytes((FOLDER / 'request.json').read_bytes())
        self.rows = [{'z': float(i + 1) / 10, 'observed': 2., 'kind': i % 3} for i in range(13)]
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.module_root = Path(temporary.name)
        source = self.module_root / '__init__.py'
        source.write_text('# synthetic runtime inventory for admission testing\n')
        inventory = [{'path': str(source), 'relative_path': '__init__.py',
                      'bytes': source.stat().st_size, 'sha256': controller.sha256(source)}]
        executable = Path(sys.executable).absolute()
        self.runtime = {'python': {'version': sys.version, 'implementation': 'CPython',
                                   'executable': {'path': str(executable), 'resolved_path': str(executable.resolve()),
                                                  'bytes': executable.stat().st_size, 'sha256': controller.sha256(executable)}},
                        'mpmath': {'version': '1.3.0', 'backend': 'python', 'module_root': str(self.module_root),
                                   'source_inventory': inventory,
                                   'source_inventory_sha256': controller.hashlib.sha256(json.dumps(inventory, sort_keys=True,
                                                                                                 separators=(',', ':')).encode()).hexdigest()}}

    def native(self):
        normalization = 13 * math.log(2 * math.pi)
        density = {'quadratic': 0., 'log_determinant': 0., 'normalization': normalization,
                   'log_density': -normalization / 2}
        state = {'availability': 1, 'status': 0, 'numerical_status': 0}
        slots = []
        for index in range(2):
            slots.append({'source_index': index, 'model_source': copy.deepcopy(self.request['models'][index]),
                          'numerical_status': 0, 'preparation_status': 0, 'density_status': 0,
                          'predictions': [1.] * 13, 'residuals': [1.] * 13,
                          'predictions_state': dict(state), 'residuals_state': dict(state),
                          'density_state': dict(state), 'projection_estimate': 1e-12,
                          'density': dict(density), 'callbacks': 10 + index,
                          'preparation_callbacks': index, 'outer_callbacks': 10,
                          'momentum_callbacks': index})
        return {'schema_version': 1, 'method': controller.METHOD, 'mapping': controller.MAPPING,
                'requested': 7, 'batch_status': 0, 'numerical_status': 0,
                'arithmetic': dict(controller.ARITHMETIC),
                'producer_policy': dict(self.request['producer_policy']),
                'policy_metadata': dict(controller.POLICY_METADATA),
                'model_order': list(controller.MODEL_ORDER), 'model_sources': copy.deepcopy(self.request['models']),
                'queries': controller.queries(self.rows), 'slots': slots, 'accepted': True,
                'callbacks': 21, 'preparation_callbacks': 1, 'outer_callbacks': 20, 'momentum_callbacks': 1}

    def reference(self):
        # Decimal normalization satisfies its own density identity exactly here.
        normalization = Decimal.from_float(13 * math.log(2 * math.pi))
        with localcontext() as context:
            context.prec = 180
            density = {'quadratic': '0', 'log_determinant': '0', 'normalization': str(normalization),
                       'log_density': str(-normalization / 2)}
        out = {'schema_version': 1, 'mpmath_version': '1.3.0',
               'method': controller.REFERENCE_METHOD, 'mapping': controller.MAPPING,
               'constants': controller.REFERENCE_CONSTANTS,
               'model_order': list(controller.MODEL_ORDER), 'models': copy.deepcopy(self.request['models']),
               'query_order': [row['id'] for row in controller.queries(self.rows)],
               'queries': copy.deepcopy(self.rows), 'sources': copy.deepcopy(controller.INPUT_IDENTITIES),
               'request_sha256': controller.REQUEST_SHA256, 'reference_input_sha256': 'a' * 64,
               'request_canonical_sha256': controller.hashlib.sha256(json.dumps(self.request, sort_keys=True,
                                                                              separators=(',', ':')).encode()).hexdigest(),
               'reference_script_sha256': 'b' * 64, 'runtime_before': copy.deepcopy(self.runtime),
               'runtime_after': copy.deepcopy(self.runtime), 'settings': {}}
        for name, policy in controller.REFERENCE_SETTINGS.items():
            out['settings'][name] = dict(policy)
            out[name] = [{'model': model, 'ruler_mpc': '140', 'predictions': ['1'] * 13,
                          'density': dict(density), 'tail_bounds': {'scaled_expansion': '0',
                                                                    'relative_scaled_expansion': '0',
                                                                    'predictions': ['0'] * 13,
                                                                    'quadratic': '0', 'log_density': '0'},
                          **policy} for model in controller.MODEL_ORDER]
        return out

    def test_request_hash_rejects_actual_changed_bytes_and_ignored_fields(self):
        blob = (FOLDER / 'request.json').read_bytes()
        with self.assertRaisesRegex(ValueError, 'request SHA'):
            controller.validate_request_bytes(blob + b' ')
        for mutate in (lambda q: q['producer_policy'].update(relative=1e-8),
                       lambda q: q['models'][0].update(Tcmb=2.8),
                       lambda q: q['models'].reverse(),
                       lambda q: q.update(ignored_extra_field=True)):
            request = copy.deepcopy(self.request)
            mutate(request)
            with self.assertRaisesRegex(ValueError, 'request SHA'):
                controller.validate_request_bytes(json.dumps(request).encode())

    def test_input_declaration_order_roles_and_source_identities(self):
        packet = controller.load(FOLDER / 'experiment.json')
        controller.validate_admission(packet)
        controller.input_sources(packet)
        for mutate in (lambda p: p['inputs'].reverse(),
                       lambda p: p['inputs'][0].update(role='synthetic_control'),
                       lambda p: p['inputs'][0].update(sha256='0' * 64),
                       lambda p: p['inputs'][0].update(source='unknown'),
                       lambda p: p['inputs'][0].update(bytes=True),
                       lambda p: p['inputs'].append(p['inputs'][0]),
                       lambda p: p['inputs'][1].update(path=p['inputs'][0]['path'])):
            bad = copy.deepcopy(packet)
            mutate(bad)
            with self.assertRaises(ValueError):
                controller.input_sources(bad)
        for mutate in (lambda p: p['origin'].update(sha256='0' * 64),
                       lambda p: p.update(status='completed'),
                       lambda p: p.update(execution={'operation': 'lcdm-campaign.thermal-control'})):
            bad = copy.deepcopy(packet)
            mutate(bad)
            with self.assertRaises(ValueError):
                controller.validate_admission(bad)

    def test_consumed_input_bytes_are_checked(self):
        packet = controller.load(FOLDER / 'experiment.json')
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / packet['inputs'][0]['path']
            path.parent.mkdir(parents=True)
            path.write_bytes(b'changed scientifically meaningful bytes')
            with self.assertRaisesRegex(ValueError, 'changed exact input'):
                controller.inputs(root, packet)

    def test_native_closed_contract_status_states_axes_and_work(self):
        controller.check_native(self.native(), self.request, self.rows)
        mutations = [
            lambda v: v.update(method='unrelated'), lambda v: v.update(mapping='unrelated'),
            lambda v: v.update(accepted=False), lambda v: v.update(accepted=1),
            lambda v: v.update(batch_status=1), lambda v: v.update(numerical_status=True),
            lambda v: v.update(requested=3), lambda v: v['arithmetic'].update(round_to_nearest=False),
            lambda v: v['producer_policy'].update(relative=1e-8),
            lambda v: v['policy_metadata'].update(momentum_maximum_depth=31),
            lambda v: v['model_order'].reverse(), lambda v: v['model_sources'].reverse(),
            lambda v: v['queries'].reverse(), lambda v: v['queries'][0].update(kind=2),
            lambda v: v['slots'][0].update(source_index=1),
            lambda v: v['slots'][0]['model_source'].update(H0=69.),
            lambda v: v['slots'][0].update(preparation_status=1),
            lambda v: v['slots'][0]['predictions_state'].update(availability=0),
            lambda v: v['slots'][0]['density_state'].update(numerical_status=4),
            lambda v: v['slots'][0].update(density_status=None),
            lambda v: v['slots'][0]['predictions'].pop(),
            lambda v: v['slots'][0]['residuals'].__setitem__(0, 1.1),
            lambda v: v['slots'][0].update(projection_estimate=-1.),
            lambda v: v['slots'][0].update(projection_estimate=None),
            lambda v: v['slots'][0].update(callbacks=True),
            lambda v: v['slots'][0].update(callbacks=500000001),
            lambda v: v.update(callbacks=22),
            lambda v: v['slots'][1].update(preparation_callbacks=2),
            lambda v: v['slots'][0]['density'].update(log_density=0.),
            lambda v: v['slots'][0]['density'].update(normalization=0.),
        ]
        for mutate in mutations:
            native = self.native()
            mutate(native)
            with self.subTest(mutate=mutate), self.assertRaises(ValueError):
                controller.check_native(native, self.request, self.rows)

    def test_native_nonfinite_and_boolean_numbers(self):
        for value in (float('nan'), float('inf'), float('-inf'), True, '1'):
            for field in ('projection_estimate', 'density', 'predictions'):
                native = self.native()
                if field == 'density':
                    native['slots'][0]['density']['quadratic'] = value
                elif field == 'predictions':
                    native['slots'][0]['predictions'][0] = value
                else:
                    native['slots'][0][field] = value
                with self.subTest(value=value, field=field), self.assertRaises(ValueError):
                    controller.check_native(native, self.request, self.rows)

    def test_tiny_reference_refinement_and_exact_native_difference_retained(self):
        reference = self.reference()
        reference['coarse'][0]['predictions'][0] = '1.000000000000000000000000000000000000000000000000000000000001'
        checks = controller.compare(self.native(), reference, self.request, self.rows)
        self.assertEqual(Decimal(checks[0]['reference_refinement']), Decimal('1e-60'))
        self.assertNotEqual(checks[0]['reference_refinement'], '0')
        reference['fine'][0]['predictions'][0] = '1.000000000000000000000000000000000000000000000000000000000002'
        checks = controller.compare(self.native(), reference, self.request, self.rows)
        self.assertEqual(Decimal(checks[0]['difference']), Decimal('2e-60'))
        # Native decimal presentation 0.1 is admitted as its actual binary64 value.
        self.assertEqual(controller.finite(.1, 'native'), Decimal.from_float(.1))

    def test_reference_model_order_settings_lengths_and_decimal_domains(self):
        for mutate in (lambda v: v['fine'].reverse(), lambda v: v['queries'].reverse(),
                       lambda v: v['settings']['coarse'].update(digits=90),
                       lambda v: v['fine'][0].update(momentum_tail=128),
                       lambda v: v['fine'][0]['predictions'].pop(),
                       lambda v: v['runtime_after'].update(changed=True),
                       lambda v: v['fine'][0]['density'].update(log_density='0')):
            reference = self.reference()
            mutate(reference)
            with self.assertRaises(ValueError):
                controller.compare(self.native(), reference, self.request, self.rows)
        for value in ('NaN', 'Infinity', '-Infinity', '1e1000000000', ' 1 ', '0x1', True, 1.):
            with self.subTest(value=value), self.assertRaises(ValueError):
                controller.decimal_string(value, 'reference')
        for text in ('{"a":1,"a":2}', '{"a":NaN}', '{"a":1e999}'):
            with self.assertRaises(ValueError):
                controller.parse(text)

    def test_reference_refinement_and_native_budget_fail_separately(self):
        reference = self.reference()
        reference['coarse'][0]['predictions'][0] = '1.000000001'
        with self.assertRaisesRegex(controller.ComparisonFailure, 'refinement') as failure:
            controller.compare(self.native(), reference, self.request, self.rows)
        self.assertEqual(Decimal(failure.exception.checks[0]['reference_refinement']), Decimal('1e-9'))
        native = self.native()
        native['slots'][0]['predictions'][0] = 1.00000001
        native['slots'][0]['residuals'][0] = 2. - native['slots'][0]['predictions'][0]
        with self.assertRaisesRegex(ValueError, 'native comparison'):
            controller.compare(native, self.reference(), self.request, self.rows)

    def test_reference_cross_axes_do_not_accept_cancellation_or_missing_tail_bounds(self):
        for mutate in (lambda v: v['momentum_refined'][0]['predictions'].__setitem__(0, '1.000000001'),
                       lambda v: v['outer_refined'][0]['predictions'].__setitem__(0, '0.999999999'),
                       lambda v: v['fine'][1]['tail_bounds']['predictions'].__setitem__(0, '1e-9'),
                       lambda v: v['fine'][1]['tail_bounds']['predictions'].pop(),
                       lambda v: v['fine'][0]['tail_bounds'].update(quadratic='1e-50')):
            reference = self.reference()
            mutate(reference)
            with self.assertRaises(ValueError):
                controller.compare(self.native(), reference, self.request, self.rows)
        reference = self.reference()
        reference['outer_refined'][1]['predictions'][0] = '1.00000000002'
        reference['momentum_refined'][1]['predictions'][0] = '0.99999999998'
        with self.assertRaisesRegex(ValueError, 'axis/tail'):
            controller.compare(self.native(), reference, self.request, self.rows)

    def test_runtime_complete_inventory_and_changed_module_bytes(self):
        controller.verify_runtime(self.runtime, Path(sys.executable))
        for mutate in (lambda r: r['mpmath'].update(backend='gmpy'),
                       lambda r: r['python']['executable'].update(sha256='0' * 64),
                       lambda r: r['mpmath']['source_inventory'].clear(),
                       lambda r: r['mpmath']['source_inventory'][0].update(bytes=True)):
            runtime = copy.deepcopy(self.runtime)
            mutate(runtime)
            with self.assertRaises(ValueError):
                controller.verify_runtime(runtime, Path(sys.executable))
        source = self.module_root / '__init__.py'
        source.write_text('# changed source\n')
        with self.assertRaisesRegex(ValueError, 'inventory'):
            controller.verify_runtime(self.runtime, Path(sys.executable))

    def test_fresh_failure_is_immutable_and_path_names_are_bounded(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            args = SimpleNamespace(name='failed', engine_source=root, input_root=root,
                                   reference_python=Path(sys.executable), engine_artifacts=None)
            identity = {'head': 'test', 'status': '', 'files': {}}
            with patch.object(controller, 'ROOT', root), patch.object(controller, 'source_identity', return_value=identity), \
                    patch.object(controller, 'read_packet', side_effect=ValueError('failed input admission')):
                self.assertEqual(controller.execute(args), 1)
                record = root / 'results/lcdm-campaign/failed/record.json'
                result = json.loads(record.read_text())
                self.assertEqual(result['status'], 'failed')
                self.assertEqual(result['gates']['execution'], 'unassessed')
                self.assertEqual(result['gates']['numerical'], 'unassessed')
                self.assertEqual(result['source_before'], result['source_after'])
                self.assertEqual(record.stat().st_mode & 0o222, 0)
                before = record.read_bytes()
                with self.assertRaises(FileExistsError):
                    controller.execute(args)
                self.assertEqual(record.read_bytes(), before)
                for name in ('../escape', 'sub/path', '/absolute', 'a' * 97, ''):
                    args.name = name
                    with self.assertRaises(ValueError):
                        controller.execute(args)
            outside = root / 'outside'
            outside.mkdir()
            (root / 'results/lcdm-campaign/escape').symlink_to(outside)
            args.name = 'escape'
            with patch.object(controller, 'ROOT', root), self.assertRaises(FileExistsError):
                controller.execute(args)
            external = tempfile.TemporaryDirectory()
            self.addCleanup(external.cleanup)
            (root / 'results/lcdm-campaign/external').symlink_to(external.name)
            args.name = 'external'
            with patch.object(controller, 'ROOT', root), self.assertRaises(ValueError):
                controller.execute(args)

    def test_request_changed_during_admission_is_retained_before_any_job(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            folder = root / 'experiments/lcdm-campaign'
            folder.mkdir(parents=True)
            for name in ('request.json', 'experiment.json', 'candidate.json'):
                (folder / name).write_bytes((FOLDER / name).read_bytes())
            packet = controller.load(folder / 'experiment.json')
            identity = {'head': 'reviewed-source', 'status': '', 'files': {}}
            def changed(_folder):
                request = controller.load(folder / 'request.json')
                request['producer_policy']['relative'] = 1e-8
                (folder / 'request.json').write_text(json.dumps(request))
                return packet, None, {}
            args = SimpleNamespace(name='changed-request', engine_source=root, input_root=root,
                                   reference_python=Path(sys.executable), engine_artifacts=None)
            with patch.object(controller, 'ROOT', root), patch.object(controller, 'FOLDER', folder), \
                    patch.object(controller, 'source_identity', return_value=identity), \
                    patch.object(controller, 'read_packet', side_effect=changed), patch.object(controller, 'child') as child:
                self.assertEqual(controller.execute(args), 1)
                child.assert_not_called()
            record = root / 'results/lcdm-campaign/changed-request/record.json'
            result = json.loads(record.read_text())
            self.assertIn('request SHA', result['error'])
            self.assertEqual(result['status'], 'failed')
            self.assertEqual(record.stat().st_mode & 0o222, 0)

    def test_committed_snapshot_rejects_mutation_inventory_and_symlink(self):
        for mode in ('good', 'changed', 'missing', 'symlink'):
            with tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                store = root / 'results'
                store.mkdir()
                original = b'reviewed source'
                (root / 'consumer.cpp').write_bytes(b'changed' if mode == 'changed' else original)
                identity = {'head': 'reviewed', 'status': '', 'files': {'consumer.cpp': controller.hashlib.sha256(original).hexdigest()}}
                def archive(*args, **kwargs):
                    with tarfile.open(store / 'source.tar', 'w') as output:
                        if mode == 'missing':
                            return
                        member = tarfile.TarInfo('consumer.cpp')
                        member.size = len(original)
                        if mode == 'symlink':
                            member.type = tarfile.SYMTYPE
                            member.linkname = '/etc/passwd'
                            output.addfile(member)
                        else:
                            output.addfile(member, io.BytesIO(original))
                with patch.object(controller, 'ROOT', root), patch.object(controller, 'child', side_effect=archive):
                    if mode == 'good':
                        destination = controller.snapshot_source(store, identity, {})
                        self.assertEqual((destination / 'consumer.cpp').read_bytes(), original)
                    else:
                        with self.assertRaises(ValueError):
                            controller.snapshot_source(store, identity, {})
        with self.assertRaisesRegex(ValueError, 'clean committed'):
            controller.snapshot_source(Path('.'), {'status': '?? consumer.cpp'}, {})

    def test_child_failure_timeout_and_partial_json_are_preserved(self):
        with tempfile.TemporaryDirectory() as temporary:
            store = Path(temporary)
            limits = {'native_seconds': 2, 'memory_bytes': 1073741824, 'output_bytes': 65536}
            record = {}
            command = [sys.executable, '-c', "import json,sys; print(json.dumps({'accepted':False,'slots':[{'numerical_status':4}]})); print('native work failure',file=sys.stderr); sys.exit(3)"]
            with self.assertRaisesRegex(ValueError, 'status 3'):
                controller.child(command, store, 'native', limits, record=record)
            controller.ingest_outputs(store, record)
            self.assertIs(record['native']['accepted'], False)
            self.assertEqual(record['native']['slots'][0]['numerical_status'], 4)
            self.assertEqual(record['subprocesses']['native']['returncode'], 3)
            self.assertIn('native work failure', (store / 'native.err').read_text())
            limits['timeout_seconds'] = .05
            with self.assertRaisesRegex(ValueError, 'timeout'):
                controller.child([sys.executable, '-c', "import time; print('partial',flush=True); time.sleep(3)"],
                                 store, 'timeout', limits, record=record)
            self.assertTrue(record['subprocesses']['timeout']['timed_out'])
            self.assertIn('partial', (store / 'timeout.out').read_text())
            self.assertTrue((store / 'timeout.err').is_file())
            limits['overflow_seconds'] = 2
            limits['output_bytes'] = 64
            with self.assertRaisesRegex(ValueError, 'overflow failed'):
                controller.child([sys.executable, '-c', "print('x'*100000)"], store, 'overflow', limits, record=record)
            self.assertLessEqual((store / 'overflow.out').stat().st_size, 64)
            self.assertNotEqual(record['subprocesses']['overflow']['returncode'], 0)

    def test_engine_artifacts_source_inventory_and_real_changed_bytes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / 'source'
            artifacts = root / 'artifacts'
            source.mkdir()
            names = ('Cargo.toml', 'Cargo.lock', 'build.rs', 'cpp/CMakeLists.txt', 'cpp/include/irred/unused.hpp')
            for name in names:
                path = source / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text('reviewed ' + name)
            library = root / 'standard-library'
            library.write_bytes(b'pinned standard library bytes')
            paths = {'archive': artifacts / 'build/native-release/libirred_core.a',
                     'cli': artifacts / 'target/release/irred',
                     'manifest': artifacts / 'build/build-manifest-release.json'}
            for name in ('archive', 'cli'):
                paths[name].parent.mkdir(parents=True, exist_ok=True)
                paths[name].write_bytes(('pinned ' + name).encode())
            build = {'profile': 'release', 'backend': 'portable_cpu',
                     'sources': {name: controller.sha256(source / name) for name in names},
                     'compiler_executable_digest': controller.sha256(Path('/usr/bin/c++').resolve()),
                     'standard_library': str(library), 'standard_library_digest': controller.sha256(library),
                     'tool_executable_digests': {}}
            build_id = controller.hashlib.sha256(json.dumps(build, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
            build.update(git_head='reviewed-engine', git_status='', build_id=build_id)
            paths['manifest'].write_text(json.dumps(build))
            identity = {'revision': 'reviewed-engine', 'build_id': build_id,
                        **{name + '_sha256': controller.sha256(path) for name, path in paths.items()}}
            with patch.object(controller, 'git', side_effect=lambda _root, *args: 'reviewed-engine' if args[0] == 'rev-parse' else ''):
                admitted = controller.fingerprint(source, identity, artifacts)
                self.assertEqual(admitted['artifact_path'], str(artifacts))
                self.assertEqual(admitted['sources'], build['sources'])
                for name in ('archive', 'cli', 'manifest'):
                    original = paths[name].read_bytes()
                    paths[name].write_bytes(original + b'changed')
                    with self.subTest(name=name), self.assertRaisesRegex(ValueError, 'changed engine'):
                        controller.fingerprint(source, identity, artifacts)
                    paths[name].write_bytes(original)
                header = source / 'cpp/include/irred/unused.hpp'
                original = header.read_bytes()
                header.write_bytes(b'changed unused header')
                with self.assertRaisesRegex(ValueError, 'changed build source'):
                    controller.fingerprint(source, identity, artifacts)
                header.unlink()
                with self.assertRaisesRegex(ValueError, 'inventory'):
                    controller.fingerprint(source, identity, artifacts)
                header.write_bytes(original)
                extra = source / 'cpp/include/irred/extra.hpp'
                extra.write_text('not inventoried')
                with self.assertRaisesRegex(ValueError, 'inventory'):
                    controller.fingerprint(source, identity, artifacts)
                extra.unlink()
                library.write_bytes(b'changed standard library')
                with self.assertRaisesRegex(ValueError, 'compiler/standard library'):
                    controller.fingerprint(source, identity, artifacts)


if __name__ == '__main__':
    unittest.main()
