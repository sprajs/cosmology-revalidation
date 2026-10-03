"""Synthetic transport fixtures only; never run the original reference fit.

Parent may copy this into tests/ and select the source using
TEMPORAL_REFERENCE_SOURCE. No archived reference/native numerical result is used.
"""
import contextlib
import copy
import decimal
import hashlib
import importlib.util
import io
import json
import os
import pathlib
import re
import sys
import tempfile
import types
import unittest
from unittest import mock


REFERENCE_FILE = pathlib.Path(os.environ.get(
    "TEMPORAL_REFERENCE_SOURCE",
    str(pathlib.Path(__file__).resolve().parents[1] / "experiments" / "temporal-shared-optical-control" / "reference.py")))
SPEC = importlib.util.spec_from_file_location("temporal_reference_transport_fixture", REFERENCE_FILE)
REFERENCE = importlib.util.module_from_spec(SPEC)


def forbidden_science(*args, **kwargs):
    raise AssertionError("Synthetic transport fixture forbids scientific initialization/GL/Fourier/reference calculation")


def metadata_module(name):
    module = types.ModuleType(name)
    module.__getattr__ = forbidden_science
    return module


# The ordinary repository test suite requires no science package. These scoped
# metadata fixtures supply only attributes used during module definition, never
# numerical functions. The actual getter is a separate pinned-runtime check.
FIXTURE_MPMATH = metadata_module("mpmath")
FIXTURE_MPMATH.__path__ = []
FIXTURE_MPMATH.__file__ = "/synthetic-fixture/mpmath/__init__.py"
FIXTURE_MPMATH.__version__ = "1.3.0"
FIXTURE_MPMATH.mp = types.SimpleNamespace(dps=60, prec=203)
FIXTURE_LIBMP = metadata_module("mpmath.libmp")
FIXTURE_LIBMP.__path__ = []
FIXTURE_BACKEND = metadata_module("mpmath.libmp.backend")
FIXTURE_BACKEND.BACKEND = "python"
FIXTURE_LIBMP.backend = FIXTURE_BACKEND
FIXTURE_MPMATH.libmp = FIXTURE_LIBMP
with mock.patch.dict(sys.modules, {"mpmath": FIXTURE_MPMATH, "mpmath.libmp": FIXTURE_LIBMP,
                                  "mpmath.libmp.backend": FIXTURE_BACKEND}):
    SPEC.loader.exec_module(REFERENCE)
for forbidden_name in ("_initialize_scientific_constants", "rule", "fourier", "run"):
    setattr(REFERENCE, forbidden_name, forbidden_science)


def reviewed_request():
    value = REFERENCE.strict_json(REFERENCE._EXPECTED_REQUEST_TEXT.encode("utf8"))
    value.update(request_id="temporal-shared-optical-reproducible-request/v1",
                 interface_id="temporal-shared-optical-bounded-controller/v1",
                 status="reviewed-native-synthetic-control", execution=None)
    for name in ("reference", "controller", "sdk_verifier"):
        value["source_ports"][name]["sha256"] = "a" * 64
    return value


class SourcePreservationTests(unittest.TestCase):
    def test_fixture_loader_cannot_call_scientific_routes(self):
        self.assertIs(REFERENCE.mp, FIXTURE_MPMATH)
        self.assertIs(REFERENCE.mp_backend, FIXTURE_BACKEND)
        for name in ("_initialize_scientific_constants", "rule", "fourier", "run"):
            with self.assertRaises(AssertionError):
                getattr(REFERENCE, name)()
        with self.assertRaises(AssertionError):
            REFERENCE.mp.gauss_quadrature

    def test_original_scientific_functions_and_initializer_sequence(self):
        source = REFERENCE_FILE.read_text(encoding="utf8")
        helpers = source[source.index("def reserve(x):\n"):source.index("def run():\n")]
        self.assertEqual(hashlib.sha256(helpers.encode("utf8")).hexdigest(),
                         "77fbb65f144782178e3f14d82a9a75f20be106911ab0284995813ee794274739")
        run = source[source.index("def run():\n"):source.index("# All helpers below implement")]
        run = re.sub(r"^[ \t]*_observe_partial\([^\n]*\)\n", "", run, flags=re.MULTILINE).rstrip() + "\n"
        self.assertEqual(hashlib.sha256(run.encode("utf8")).hexdigest(),
                         "21d5a5aeef3c0d587c13374fb3cd2dec99f93fb2ca61b0d857daf805419a9016")
        initialized = source[source.index("    mpf = mp.mpf\n"):source.index("\n\n", source.index("    mpf = mp.mpf\n"))]
        initialized = "\n".join(line[4:] for line in initialized.splitlines()) + "\n"
        self.assertEqual(hashlib.sha256(initialized.encode("utf8")).hexdigest(),
                         "5958939f4e1385ccf08f21cf365f6e5776865742a5d21ec152c7df2844f12c6d")

    def test_runtime_branch_never_initializes_or_calls_science(self):
        with mock.patch.object(REFERENCE, "_bounded_limits"), \
             mock.patch.object(REFERENCE, "runtime_fingerprint", return_value={"fixture": "metadata"}), \
             mock.patch.object(REFERENCE, "_initialize_scientific_constants", side_effect=AssertionError("science forbidden")) as initialize, \
             mock.patch.object(REFERENCE, "run", side_effect=AssertionError("science forbidden")) as route, \
             contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(REFERENCE.main(["--runtime-fingerprint"]), 0)
        self.assertEqual((initialize.call_count, route.call_count), (0, 0))


class ClosedJSONTests(unittest.TestCase):
    def test_decimal_transport_retains_emitted_digits(self):
        value = REFERENCE.strict_json(b'{"native":1.01000000000000000888,"flag":true,"index":1}')
        self.assertEqual(value["native"], decimal.Decimal("1.01000000000000000888"))
        self.assertIs(type(value["native"]), decimal.Decimal)
        self.assertIs(type(value["flag"]), bool)
        self.assertIs(type(value["index"]), int)
        self.assertEqual(REFERENCE.strict_json(REFERENCE._json_bytes(value)), value)

    def test_duplicate_nonfinite_excessive_and_deep_input_refuse(self):
        for raw in (b'{"x":0,"x":1}', b'{"x":NaN}', b'{"x":Infinity}', b'{"x":1e1001}',
                    b"[" * 33 + b"0" + b"]" * 33, b"0" * 257, b"{}x", b"\xff"):
            with self.subTest(raw=raw[:40]), self.assertRaises((ValueError, UnicodeError)):
                REFERENCE.strict_json(raw)

    def test_unknown_keys_bool_integer_and_order_remain_distinct(self):
        for actual in ({"x": True}, {"x": decimal.Decimal(1)}, {"x": 1, "extra": 0}):
            with self.assertRaises(REFERENCE.TransportRefusal):
                REFERENCE._same_constants({"x": 1}, actual)
        with self.assertRaises(REFERENCE.TransportRefusal):
            REFERENCE._same_constants(["S0", "S1"], ["S1", "S0"])

    def test_frozen_request_allows_only_actual_script_hash_slots(self):
        value = reviewed_request()
        REFERENCE._request_expected(value)
        for mutate in (
            lambda x: x.update(execution="run"),
            lambda x: x["resources"].update(reference_cpu_seconds=181),
            lambda x: x["shape_contract"]["unchanged_scientific_settings"].update(decimal_digits=90),
            lambda x: x["source_ports"]["reference"].update(sha256=None),
            lambda x: x.update(extra=True),
        ):
            changed = copy.deepcopy(value)
            mutate(changed)
            with self.assertRaises(REFERENCE.TransportRefusal):
                REFERENCE._request_expected(changed)

    def test_output_bound_and_decimal_grammar(self):
        with self.assertRaises(REFERENCE.TransportRefusal):
            REFERENCE._json_bytes(["x" * 8192] * 129)
        for token in ("+1", "01", "NaN", "1_0", " 1", "1e1001"):
            with self.assertRaises(REFERENCE.TransportRefusal):
                REFERENCE._finite_decimal(token)


class PrefixAndErrorTests(unittest.TestCase):
    def setUp(self):
        self.previous = REFERENCE._PARTIAL
        REFERENCE._PARTIAL = REFERENCE._empty_partial()

    def tearDown(self):
        REFERENCE._PARTIAL = self.previous

    def test_sigma_and_record_prefix_order_preserves_unearned_group(self):
        REFERENCE._observe_partial("sigma-start", {"noise": 0})
        row = {"count_probability": "1", "count_probability_error": "0", "below": "0", "below_error": "0"}
        REFERENCE._observe_partial("sigma-row", {"sigma_rows": [row]})
        self.assertEqual(REFERENCE._PARTIAL["active_sigma_rows"], [row])
        self.assertEqual(REFERENCE._PARTIAL["detector_prefix"], [])
        group = {"sigma_electrons": 0, "rows": [row] * 8}
        REFERENCE._observe_partial("sigma-complete", {"noise": 0, "detector_rows": [group]})
        REFERENCE._observe_partial("record", {"results": [{"fixture_record": i} for i in range(3)]})
        self.assertEqual(REFERENCE._PARTIAL["completed_groups"], ["detector-sigma0"])
        self.assertIsNone(REFERENCE._PARTIAL["active_sigma_index"])
        self.assertEqual(len(REFERENCE._PARTIAL["record_prefix"]), 3)
        self.assertNotIn("records", REFERENCE._PARTIAL["completed_groups"])

    def test_existing_arithmetic_failure_details_are_retained(self):
        detail = {"stage": "positive-mixture-error-gate", "value": "0", "error": "1",
                  "conditional_values": ["0", "0"], "conditional_errors": ["1", "1"]}
        result = REFERENCE._error(ArithmeticError(json.dumps(detail)))
        self.assertEqual(result["details"], detail)
        self.assertEqual(result["stage"], detail["stage"])
        self.assertEqual(result["kind"], "reference")
        self.assertFalse(result["message_truncated"])

    def test_long_error_keeps_digest_and_bounded_message(self):
        message = "error" * 2000
        result = REFERENCE._error(ValueError(message))
        self.assertTrue(result["message_truncated"])
        self.assertLessEqual(len(result["message"].encode("utf8")), 4096)
        self.assertEqual(result["message_sha256"], hashlib.sha256(message.encode("utf8")).hexdigest())

    def test_conditional_positive_failure_binds_sigma_row_and_exact_indices(self):
        for sigma in (0, 1):
            row = ({"count_probability": "0", "count_probability_error": "1", "below": "0", "below_error": "1"}
                   if sigma == 0 else {"density32_per_adu": "-1", "density64_per_adu": "0",
                                      "density_error_per_adu": "1", "below32": "0", "below64": "0",
                                      "below_error": "1", "density_tail_per_adu": "0", "cdf_tail": "0"})
            detail = {"stage": "conditional-positive-reference-gate", "sigma_index": sigma, "channel_index": 7,
                      "lambda": "1", "density": "0", "density_error": "1", "nondetection": "0",
                      "nondetection_error": "1", "coarse_refined_attempt": row}
            self.assertTrue(REFERENCE._failure_details(detail))
            self.assertEqual(REFERENCE._error(ArithmeticError(json.dumps(detail)))["details"], detail)
            for field, wrong in (("sigma_index", True), ("sigma_index", 2), ("channel_index", True),
                                 ("channel_index", 8), ("density_error", "-1"), ("coarse_refined_attempt", {})):
                changed = dict(detail, **{field: wrong})
                self.assertFalse(REFERENCE._failure_details(changed))


class RuntimeCollectorTests(unittest.TestCase):
    @staticmethod
    def builtin():
        return types.SimpleNamespace(__spec__=types.SimpleNamespace(origin="built-in"))

    def test_module_set_cap_is_distinct_from_file_union(self):
        for count, refused in ((4096, False), (4097, True)):
            modules = {"fixture_%04d" % i: self.builtin() for i in range(count)}
            with mock.patch.object(REFERENCE.sys, "modules", modules):
                if refused:
                    with self.assertRaises(REFERENCE.TransportRefusal):
                        REFERENCE._collect_imported_modules(set())
                else:
                    expanded = set()
                    result, _ = REFERENCE._collect_imported_modules(expanded)
                    self.assertEqual(len(result), count)
                    self.assertEqual(expanded, set())

    def test_base_alias_must_resolve_to_pinned_interpreter_prefix(self):
        pinned = pathlib.Path(REFERENCE.RESOLVED_PYTHON).parent.parent
        for destination, refused in ((pinned, False), (pinned.parent / "future-version", True)):
            def resolve(path, strict=False):
                self.assertTrue(strict)
                if str(path) == REFERENCE.PYTHON_PATH:
                    return pathlib.Path(REFERENCE.RESOLVED_PYTHON)
                self.assertEqual(str(path), REFERENCE.BASE_PREFIX)
                return destination
            with mock.patch.object(REFERENCE.sys, "executable", REFERENCE.PYTHON_PATH), \
                 mock.patch.object(REFERENCE.sys, "version", REFERENCE.PYTHON_VERSION), \
                 mock.patch.object(REFERENCE.sys, "prefix", REFERENCE.VENV_PREFIX), \
                 mock.patch.object(REFERENCE.sys, "base_prefix", REFERENCE.BASE_PREFIX), \
                 mock.patch.object(REFERENCE.platform, "python_implementation", return_value="CPython"), \
                 mock.patch.object(REFERENCE.pathlib.Path, "resolve", resolve):
                if refused:
                    with self.assertRaises(REFERENCE.TransportRefusal):
                        REFERENCE._runtime_python_identity()
                else:
                    REFERENCE._runtime_python_identity()

    @staticmethod
    def map_bytes(paths):
        return "".join("1000-2000 r-xp 00000000 00:00 1 " + str(path) + "\n" for path in paths).encode("utf8")

    def test_mapped_cap_deleted_and_unreadable_paths_refuse(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            paths = []
            for i in range(257):
                path = root / ("lib%03d.so" % i)
                path.write_bytes(b"library fixture")
                paths.append(path)
            for selected, refused in ((paths[:256], False), (paths, True)):
                with mock.patch.object(REFERENCE, "_proc_bytes", return_value=self.map_bytes(selected)):
                    if refused:
                        with self.assertRaises(REFERENCE.TransportRefusal):
                            REFERENCE._collect_mapped_libraries()
                    else:
                        self.assertEqual(len(REFERENCE._collect_mapped_libraries()), 256)
            for malformed in ((str(paths[0]) + " (deleted)",), (str(root / "missing.so"),)):
                with mock.patch.object(REFERENCE, "_proc_bytes", return_value=self.map_bytes(malformed)), \
                     self.assertRaises((REFERENCE.TransportRefusal, OSError)):
                    REFERENCE._collect_mapped_libraries()

    def fingerprint_fixture(self, root, modules, mappings, import_during_hash=False):
        executable, config, init = root / "python", root / "pyvenv.cfg", root / "__init__.py"
        for path in (executable, config, init):
            path.write_bytes(b"runtime fixture")
        artifact_sha = "a" * 64
        def discovery():
            # Both libraries already belong to the file union, even when one map disappears.
            return {init}, {init, root / "a-" / "m.so", root / "a" / "m.so"}
        def identity(path):
            if import_during_hash:
                modules["late_helper"] = self.builtin()
            return {"path": str(path), "bytes": 1, "sha256": artifact_sha}
        def proc(path):
            return b"Threads:\t1\n" if path.endswith("status") else self.map_bytes(mappings)
        with mock.patch.object(REFERENCE, "_metadata_discovery", side_effect=discovery), \
             mock.patch.object(REFERENCE, "_runtime_python_identity"), \
             mock.patch.object(REFERENCE, "_proc_bytes", side_effect=proc), \
             mock.patch.object(REFERENCE, "_runtime_file", side_effect=identity), \
             mock.patch.object(REFERENCE.sys, "modules", modules), \
             mock.patch.object(REFERENCE.sys, "executable", str(executable)), \
             mock.patch.object(REFERENCE.sys, "prefix", str(root)), \
             mock.patch.object(REFERENCE.sys, "dont_write_bytecode", True), \
             mock.patch.object(REFERENCE.mp, "__file__", str(init)), \
             mock.patch.object(REFERENCE, "LEGACY_SHA", REFERENCE._sha(REFERENCE._canonical({str(init): artifact_sha}))), \
             mock.patch.object(REFERENCE, "CFG_SHA", artifact_sha), \
             mock.patch.object(REFERENCE, "PYTHON_SHA", artifact_sha), \
             mock.patch.object(REFERENCE, "MPMATH_INIT_SHA", artifact_sha), \
             mock.patch.dict(REFERENCE.os.environ, {key: "1" for key in REFERENCE.ENV_KEYS}):
            return REFERENCE.runtime_fingerprint()

    def test_same_union_cannot_hide_module_or_map_changes_and_sort_is_lexical(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            maps = [root / "a" / "m.so", root / "a-" / "m.so"]
            for path in maps:
                path.parent.mkdir()
                path.write_bytes(b"library fixture")
            modules = {"first": self.builtin()}
            before_names = tuple(modules)
            before = self.fingerprint_fixture(root, modules, maps)
            self.assertEqual(tuple(modules), before_names)  # Helpers imported no module during collection.
            self.assertEqual([record["path"] for record in before["mapped_libraries"]], sorted(map(str, maps)))
            modules["already_in_package_union"] = self.builtin()
            imported = self.fingerprint_fixture(root, modules, maps)
            self.assertEqual(before["inventory"], imported["inventory"])
            self.assertNotEqual(before["imported_modules"], imported["imported_modules"])
            disappeared = self.fingerprint_fixture(root, modules, maps[:1])
            self.assertEqual(imported["inventory"], disappeared["inventory"])
            self.assertNotEqual(imported["mapped_libraries"], disappeared["mapped_libraries"])
            with self.assertRaises(REFERENCE.TransportRefusal):
                self.fingerprint_fixture(root, modules, maps, import_during_hash=True)


class AdmissionBoundaryTests(unittest.TestCase):
    @contextlib.contextmanager
    def fixture(self):
        with tempfile.TemporaryDirectory() as temporary:
            attempt = pathlib.Path(temporary) / "results" / "temporal-shared-optical-control" / "synthetic"
            packet = attempt / "source" / "experiments" / "temporal-shared-optical-control"
            packet.mkdir(parents=True)
            request = reviewed_request()
            for name in ("consumer", "reference", "controller", "sdk_verifier"):
                leaf = request["source_ports"][name]["path"]
                raw = ("synthetic source: " + leaf).encode("ascii")
                (packet / leaf).write_bytes(raw)
                request["source_ports"][name]["sha256"] = hashlib.sha256(raw).hexdigest()
            contract = b"{}" + b" " * 26348
            (attempt / "contract.snapshot.json").write_bytes(contract)
            (attempt / "request.json").write_bytes(REFERENCE._json_bytes(request))
            # Synthetic immutable pins only in this fixture; actual admission logic remains active.
            with mock.patch.object(REFERENCE, "__file__", str(packet / "reference.py")), \
                 mock.patch.object(REFERENCE, "CONTRACT_SHA", hashlib.sha256(contract).hexdigest()), \
                 mock.patch.object(REFERENCE, "_EXPECTED_REQUEST_TEXT", REFERENCE._json_bytes(request).decode("utf8")):
                yield attempt, packet

    def test_genuine_admission_binds_consumed_buffers_and_creates_no_output(self):
        with self.fixture() as (attempt, _):
            _, identities, seal, attempt_fd, packet_fd = REFERENCE._admit(str(attempt / "request.json"), str(attempt / "reference.json"))
            try:
                self.assertEqual(identities["original_contract_sha256"], seal["current_source_sha256"]["contract.snapshot.json"])
                self.assertEqual(identities["request_sha256"], hashlib.sha256((attempt / "request.json").read_bytes()).hexdigest())
                self.assertFalse((attempt / "reference.json").exists())
            finally:
                os.close(attempt_fd)
                os.close(packet_fd)

    def test_contract_reread_drift_refuses_before_initialization_or_science(self):
        with self.fixture() as (attempt, _):
            original_read = REFERENCE._read_leaf
            contract_reads = 0
            def read(directory_fd, name, maximum=REFERENCE.JSON_LIMIT):
                nonlocal contract_reads
                raw = original_read(directory_fd, name, maximum)
                if name == "contract.snapshot.json":
                    contract_reads += 1
                    if contract_reads == 1:
                        (attempt / name).write_bytes(b"[]" + b" " * 26348)
                return raw
            with mock.patch.object(REFERENCE, "_read_leaf", side_effect=read), \
                 mock.patch.object(REFERENCE, "runtime_fingerprint", return_value={"fixture": "runtime"}), \
                 mock.patch.object(REFERENCE, "_initialize_scientific_constants") as initialize, \
                 mock.patch.object(REFERENCE, "run") as route, \
                 contextlib.redirect_stdout(io.StringIO()) as summary, contextlib.redirect_stderr(io.StringIO()) as error:
                code = REFERENCE._normal(str(attempt / "request.json"), str(attempt / "reference.json"))
            self.assertEqual((code, initialize.call_count, route.call_count, contract_reads), (1, 0, 0, 2))
            self.assertIn("contract-drift", error.getvalue())
            self.assertIsNone(REFERENCE.strict_json(summary.getvalue().encode("utf8"))["reference_json_sha256"])
            self.assertFalse((attempt / "reference.json").exists())

    def test_consumed_request_drift_and_wrong_source_bytes_refuse(self):
        with self.fixture() as (attempt, packet):
            (packet / "controller.py").write_bytes(b"changed source")
            with self.assertRaises(REFERENCE.TransportRefusal):
                REFERENCE._admit(str(attempt / "request.json"), str(attempt / "reference.json"))
        with self.fixture() as (attempt, _):
            original_read = REFERENCE._read_leaf
            request_reads = 0
            def read(directory_fd, name, maximum=REFERENCE.JSON_LIMIT):
                nonlocal request_reads
                raw = original_read(directory_fd, name, maximum)
                if name == "request.json":
                    request_reads += 1
                    if request_reads == 1:
                        (attempt / name).write_bytes(raw + b" ")
                return raw
            with mock.patch.object(REFERENCE, "_read_leaf", side_effect=read), self.assertRaises(REFERENCE.TransportRefusal):
                REFERENCE._admit(str(attempt / "request.json"), str(attempt / "reference.json"))

    def test_layout_rejects_unapproved_output_and_path_before_read(self):
        with self.fixture() as (attempt, _):
            for request, output in ((attempt / "request.json", attempt / "wrong.json"),
                                    (attempt / "wrong.json", attempt / "reference.json"),
                                    (attempt / "request.json", attempt.parent / "reference.json")):
                with self.assertRaises(REFERENCE.TransportRefusal):
                    REFERENCE._admit(str(request), str(output))


class FileAndFailureTests(unittest.TestCase):
    def test_fifo_leaf_and_runtime_file_use_nonblocking_open_and_refuse(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            os.mkfifo(root / "fifo")
            fd = REFERENCE._directory(str(root))
            original_open = os.open
            def opened(name, flags, *args, **kwargs):
                if name == "fifo":
                    self.assertTrue(flags & os.O_NONBLOCK)
                return original_open(name, flags, *args, **kwargs)
            try:
                with mock.patch.object(REFERENCE.os, "open", side_effect=opened):
                    with self.assertRaises(REFERENCE.TransportRefusal):
                        REFERENCE._read_leaf(fd, "fifo")
                    with self.assertRaises(REFERENCE.TransportRefusal):
                        REFERENCE._runtime_file(root / "fifo")
            finally:
                os.close(fd)

    def test_symlink_leaf_and_parent_are_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            (root / "real").mkdir()
            (root / "real" / "source").write_bytes(b"fixture")
            (root / "link").symlink_to(root / "real", target_is_directory=True)
            with self.assertRaises(OSError):
                REFERENCE._directory(str(root / "link"))
            (root / "alias").symlink_to(root / "real" / "source")
            fd = REFERENCE._directory(str(root))
            try:
                with self.assertRaises(OSError):
                    REFERENCE._read_leaf(fd, "alias")
            finally:
                os.close(fd)

    def exercise_normal(self, root, maps, earned_payload, existing=False, fail_science=False, source_result=None, prefix_failure=False):
        request = reviewed_request()
        identities = {"request_sha256": "b" * 64}
        baseline = {"fixture": "source-seal"}
        packet_fd, attempt_fd = REFERENCE._directory(str(root)), REFERENCE._directory(str(root))
        if existing:
            (root / "reference.json").write_bytes(b"old immutable output")
        def science():
            if prefix_failure:
                REFERENCE._observe_partial("sigma-complete", {"noise": 0, "detector_rows": [{"fixture": "earned sigma0"}]})
                REFERENCE._observe_partial("record", {"results": [{"fixture": "earned record", "index": i} for i in range(3)]})
                REFERENCE._observe_partial("sigma-start", {"noise": 1})
                REFERENCE._observe_partial("sigma-row", {"sigma_rows": [{"fixture": "passed sigma1 row"}]})
                raise ArithmeticError('{"stage":"conditional-positive-reference-gate","sigma_index":1,"channel_index":1,"lambda":"1","density":"0","density_error":"1","nondetection":"0","nondetection_error":"1","coarse_refined_attempt":{"density32_per_adu":"-1","density64_per_adu":"0","density_error_per_adu":"1","below32":"0","below64":"0","below_error":"1","density_tail_per_adu":"0","cdf_tail":"0"}}')
            if fail_science:
                REFERENCE._PARTIAL["photon_prefix"] = ["1"]
                REFERENCE._PARTIAL["photon_error_prefix"] = ["0"]
                REFERENCE._PARTIAL["frequency_prefix"] = [{"fixture": "earned prefix"}]
                raise ArithmeticError('{"stage":"positive-mixture-error-gate","value":"0","error":"1","conditional_values":["0","0"],"conditional_errors":["1","1"]}')
            return earned_payload
        output = io.StringIO()
        with mock.patch.object(REFERENCE, "_admit", return_value=(request, identities, baseline, attempt_fd, packet_fd)), \
             mock.patch.object(REFERENCE, "runtime_fingerprint", side_effect=maps), \
             mock.patch.object(REFERENCE, "_source_seal", return_value=baseline if source_result is None else source_result) as source_seal, \
             mock.patch.object(REFERENCE, "_initialize_scientific_constants") as initialize, \
             mock.patch.object(REFERENCE, "run", side_effect=science) as route, \
             contextlib.redirect_stdout(output), contextlib.redirect_stderr(io.StringIO()):
            code = REFERENCE._normal(str(root / "request.json"), str(root / "reference.json"))
        return code, REFERENCE.strict_json(output.getvalue().encode("utf8")), route.call_count, initialize.call_count, source_seal.call_count

    def test_complete_payload_survives_terminal_runtime_drift(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            code, summary, calls, initialized, source_checks = self.exercise_normal(root, [{"fixture": "before"}, {"fixture": "changed"}], {"fixture": "earned complete payload"})
            envelope = REFERENCE.strict_json((root / "reference.json").read_bytes())
            self.assertEqual(code, 1)
            self.assertEqual(envelope["status"], "refused")
            self.assertEqual(envelope["payload"], {"fixture": "earned complete payload"})
            self.assertEqual(envelope["error"]["kind"], "runtime")
            self.assertEqual((calls, initialized, source_checks), (1, 1, 1))
            self.assertEqual(summary["reference_json_sha256"], hashlib.sha256((root / "reference.json").read_bytes()).hexdigest())

    def test_original_positive_failure_keeps_prefix_and_never_fabricates_payload(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            code, _, _, _, _ = self.exercise_normal(root, [{"fixture": "same"}] * 2, None, fail_science=True)
            envelope = REFERENCE.strict_json((root / "reference.json").read_bytes())
            self.assertEqual(code, 1)
            self.assertIsNone(envelope["payload"])
            self.assertEqual(envelope["partial"]["photon_prefix"], ["1"])
            self.assertEqual(envelope["partial"]["control_failure"]["stage"], "positive-mixture-error-gate")

    def test_existing_output_never_changes_or_runs_science(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            code, summary, calls, initialized, source_checks = self.exercise_normal(root, [{"fixture": "same"}], None, existing=True)
            self.assertEqual(code, 1)
            self.assertEqual((root / "reference.json").read_bytes(), b"old immutable output")
            self.assertEqual((calls, initialized, source_checks), (0, 0, 1))
            self.assertIsNone(summary["reference_json_sha256"])

    def test_first_three_records_survive_next_sigma_failure_without_completion_marker(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            code, _, _, _, _ = self.exercise_normal(root, [{"fixture": "same"}] * 2, None, prefix_failure=True)
            envelope = REFERENCE.strict_json((root / "reference.json").read_bytes())
            partial = envelope["partial"]
            self.assertEqual(code, 1)
            self.assertIsNone(envelope["payload"])
            self.assertEqual(partial["completed_groups"], ["detector-sigma0"])
            self.assertEqual((len(partial["record_prefix"]), partial["active_sigma_index"], len(partial["active_sigma_rows"])), (3, 1, 1))
            self.assertEqual(partial["control_failure"]["channel_index"], 1)
            self.assertEqual(partial["control_failure"]["coarse_refined_attempt"]["density32_per_adu"], "-1")

    def test_terminal_source_refusal_retains_complete_payload_and_both_causes(self):
        for runtime_drift in (False, True):
            with tempfile.TemporaryDirectory() as temporary:
                root = pathlib.Path(temporary)
                maps = [{"fixture": "before"}, {"fixture": "after" if runtime_drift else "before"}]
                code, _, _, _, checks = self.exercise_normal(root, maps, {"fixture": "earned"}, source_result={"fixture": "source drift"})
                envelope = REFERENCE.strict_json((root / "reference.json").read_bytes())
                self.assertEqual((code, checks), (1, 1))
                self.assertEqual(envelope["payload"], {"fixture": "earned"})
                if runtime_drift:
                    self.assertEqual(envelope["error"]["kind"], "runtime")
                    self.assertTrue(any("Additional source terminal refusal" in x for x in envelope["qualification"]["limits"]))
                else:
                    self.assertEqual(envelope["error"]["stage"], "source-terminal-drift")

    def test_normal_oversize_serialization_preserves_empty_file_and_null_identity(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            code, summary, calls, initialized, checks = self.exercise_normal(root, [{"fixture": "same"}] * 2, {"oversized": ["x" * 8192] * 129})
            self.assertEqual((code, calls, initialized, checks), (1, 1, 1, 1))
            self.assertEqual((root / "reference.json").read_bytes(), b"")
            self.assertIsNone(summary["reference_json_sha256"])
            self.assertIsNone(summary["reference_json_bytes"])

    def test_partial_write_preserves_earned_bytes_without_retry_or_output_identity(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            original_write = os.write
            writes = []
            def write(fd, raw):
                writes.append(raw)
                if len(writes) == 1:
                    return original_write(fd, raw[:37])
                raise OSError("synthetic partial write failure")
            with mock.patch.object(REFERENCE.os, "write", side_effect=write):
                code, summary, calls, initialized, checks = self.exercise_normal(root, [{"fixture": "same"}] * 2, {"fixture": "earned"})
            self.assertEqual((code, calls, initialized, checks, len(writes)), (1, 1, 1, 1, 2))
            self.assertEqual((root / "reference.json").read_bytes(), writes[0][:37])
            self.assertEqual(writes[1], writes[0][37:])
            self.assertIsNone(summary["reference_json_sha256"])
            self.assertIsNone(summary["reference_json_bytes"])


if __name__ == "__main__":
    unittest.main()
