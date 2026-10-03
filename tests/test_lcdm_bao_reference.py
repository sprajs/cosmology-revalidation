"""Synthetic BAO axes and metadata admission; no CLASS/native Gaussian calls."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import stat
import sys
import tempfile
import types
import unittest
from unittest import mock

PACKET = Path(__file__).resolve().parents[1] / "experiments/lcdm-bao-reference"


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, PACKET / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


with mock.patch.dict(sys.modules, {"theory": types.ModuleType("theory")}):
    bao = load("bao_structure_controls", "bao.py")
controller = load("bao_authority_controls", "controller.py")

# Independent literal axes; observation and covariance fixtures are synthetic.
AXES = ((.295, "DV_over_rs"), (.510, "DM_over_rs"), (.510, "DH_over_rs"),
        (.706, "DM_over_rs"), (.706, "DH_over_rs"), (.934, "DM_over_rs"),
        (.934, "DH_over_rs"), (1.321, "DM_over_rs"), (1.321, "DH_over_rs"),
        (1.484, "DM_over_rs"), (1.484, "DH_over_rs"), (2.33, "DH_over_rs"),
        (2.33, "DM_over_rs"))
HEADER = b"# [z] [value at z] [quantity]\n"


class ReleaseHeaderControls(unittest.TestCase):
    def setUp(self):
        self.mean = "".join(f"{z} {i + 1} {kind}\n" for i, (z, kind) in enumerate(AXES)).encode()
        self.C = [[float(10 + i) if i == j else 0.0 for j in range(13)] for i in range(13)]
        self.C[0][12] = self.C[12][0] = -.125
        self.C[2][4] = self.C[4][2] = .375
        self.cov = "".join(" ".join(str(x) for x in row) + "\n" for row in self.C).encode()

    def test_exact_optional_header_preserves_order_and_all_covariance_entries(self):
        plain = bao.parse_data(self.mean, self.cov)
        labeled = bao.parse_data(HEADER + self.mean, self.cov)
        self.assertEqual(plain, labeled)
        self.assertEqual([(r["z"], r["kind"]) for r in labeled["rows"]], list(AXES))
        self.assertEqual([r["id"] for r in labeled["rows"]],
                         [f"desi-dr2-{i:02d}-{kind}" for i, (_, kind) in enumerate(AXES)])
        self.assertEqual([r["observed"] for r in labeled["rows"]], list(range(1, 14)))
        self.assertEqual(labeled["covariance"], self.C)
        self.assertEqual(sum(len(row) for row in labeled["covariance"]), 169)
        self.assertEqual(labeled["covariance"][0][12], -.125)

    def test_spoof_arbitrary_and_nonliteral_headers_refuse(self):
        for header in (b"# [z] [value] [quantity]\n", b"# arbitrary\n",
                       b" # [z] [value at z] [quantity]\n",
                       b"# [z] [value at z] [quantity]\r\n"):
            with self.subTest(header=header), self.assertRaises(ValueError):
                bao.parse_data(header + self.mean, self.cov)

    def test_repeated_header_refuses(self):
        with self.assertRaises(ValueError):
            bao.parse_data(HEADER * 2 + self.mean, self.cov)

    def test_nonleading_and_midstream_header_refuse(self):
        first, rest = self.mean.split(b"\n", 1)
        for raw in (b"\n" + HEADER + self.mean, first + b"\n" + HEADER + rest):
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                bao.parse_data(raw, self.cov)

    def test_header_does_not_relax_axes_or_full_covariance(self):
        swapped = self.mean.splitlines(keepends=True)
        swapped[1], swapped[2] = swapped[2], swapped[1]
        for raw in (b"".join(swapped), b"".join(swapped[:-1])):
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                bao.parse_data(HEADER + raw, self.cov)
        with self.assertRaises(ValueError):
            bao.parse_data(HEADER + self.mean, b"\n".join(self.cov.splitlines()[:-1]))
        asymmetric = [row[:] for row in self.C]
        asymmetric[0][12] = .125
        cov = "".join(" ".join(str(x) for x in row) + "\n" for row in asymmetric).encode()
        with self.assertRaises(ValueError):
            bao.parse_data(HEADER + self.mean, cov)


class AuthorityControls(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name).resolve()
        self.path = self.root / "opaque-source.txt"
        self.path.write_bytes(b"source")
        self.pin = {"path": str(self.path), "bytes": 6,
                    "sha256": hashlib.sha256(b"source").hexdigest()}

    def test_positive_closed_authority_and_matching_consumed_bytes(self):
        self.assertEqual(controller.positive_pin(self.pin), self.pin)
        raw, actual = controller.read(self.path, self.pin, limit=6)
        self.assertEqual(raw, b"source")
        self.assertEqual({key: actual[key] for key in self.pin}, self.pin)
        for extra in ("command", "model", "tolerance"):
            with self.subTest(extra=extra), self.assertRaises(ValueError):
                controller.positive_pin({**self.pin, extra: "unreviewed"})
        with self.assertRaises(ValueError):
            controller.positive_pin({key: value for key, value in self.pin.items() if key != "sha256"})

    def test_authority_exact_types_size_and_digest_domains(self):
        for count in (True, False, 6.0, "6", None, -1, 33554433, float("inf")):
            with self.subTest(bytes=count), self.assertRaises(ValueError):
                controller.positive_pin({**self.pin, "bytes": count})
        for digest in (None, 1, "", "x" * 64, "a" * 63, "a" * 65, "A" * 64):
            with self.subTest(digest=digest), self.assertRaises(ValueError):
                controller.positive_pin({**self.pin, "sha256": digest})

    def test_authority_relative_traversal_and_long_paths_refuse(self):
        for path in ("relative.txt", str(self.root / ".." / "elsewhere"), "/" + "a" * 4096, None):
            with self.subTest(path=path), self.assertRaises(ValueError):
                controller.positive_pin({**self.pin, "path": path})

    def test_duplicate_json_at_top_and_nested_levels_refuses(self):
        for raw in (b'{"schema":1,"schema":2}',
                    b'{"controller":{"bytes":6,"bytes":7}}'):
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                controller.decode(raw)

    def test_nonfinite_json_constants_and_overflowed_authority_refuse(self):
        for token in ("NaN", "Infinity", "-Infinity"):
            with self.subTest(token=token), self.assertRaises(ValueError):
                controller.decode(('{"bytes":' + token + '}').encode())
        # JSON exponent overflow reaches the exact integer type/domain gate.
        raw = json.dumps(self.pin).replace('"bytes": 6', '"bytes": 1e309').encode()
        with self.assertRaises(ValueError):
            controller.positive_pin(controller.decode(raw))

    def test_changed_bytes_and_byte_limit_refuse(self):
        self.path.write_bytes(b"drift!")
        with self.assertRaisesRegex(controller.SourceIdentityError, "source pin") as caught:
            controller.read(self.path, self.pin, limit=6)
        observed = caught.exception.observed_identity
        self.assertEqual(caught.exception.expected_identity, self.pin)
        self.assertEqual(observed["consumed_bytes"], 6)
        self.assertEqual(observed["consumed_sha256"], hashlib.sha256(b"drift!").hexdigest())
        self.assertTrue(observed["reached_eof"])
        self.assertEqual(observed["capture_errors"], [])
        for key in ("fd_before", "fd_after", "link_stat"):
            self.assertEqual(observed[key]["bytes"], 6)
            self.assertEqual(observed[key]["inode"], self.path.stat().st_ino)
        with self.assertRaisesRegex(controller.SourceIdentityError, "source size/type") as capped:
            controller.read(self.path, limit=5)
        prefix = capped.exception.observed_identity
        self.assertIsNone(capped.exception.expected_identity)
        self.assertEqual(prefix["consumed_bytes"], 0)
        self.assertEqual(prefix["consumed_sha256"], hashlib.sha256(b"").hexdigest())
        self.assertFalse(prefix["reached_eof"])
        self.assertEqual(prefix["fd_before"]["bytes"], 6)
        self.assertEqual(prefix["fd_after"]["bytes"], 6)
        self.assertEqual(prefix["link_stat"]["bytes"], 6)

    def test_actual_consumed_read_drift_refuses(self):
        original_read = os.read
        changed = False
        before_ns = self.path.stat().st_mtime_ns

        def changed_read(fd, count):
            nonlocal changed
            raw = original_read(fd, count)
            if raw and not changed:
                changed = True
                self.path.write_bytes(b"drift!")
                os.utime(self.path, ns=(before_ns + 1, before_ns + 1))
            return raw

        with mock.patch.object(controller.os, "read", side_effect=changed_read), \
                mock.patch.object(controller.os, "open", wraps=os.open) as opened, \
                self.assertRaisesRegex(controller.SourceIdentityError,
                                       "source changed while consumed") as caught:
            controller.read(self.path, self.pin, limit=6)
        self.assertTrue(changed)
        self.assertEqual(opened.call_count, 1)
        observed = caught.exception.observed_identity
        self.assertEqual(caught.exception.expected_identity, self.pin)
        self.assertEqual(observed["consumed_bytes"], 6)
        self.assertEqual(observed["consumed_sha256"], self.pin["sha256"])
        self.assertNotEqual(observed["consumed_sha256"], hashlib.sha256(self.path.read_bytes()).hexdigest())
        self.assertTrue(observed["reached_eof"])
        self.assertEqual(observed["fd_before"]["mtime_ns"], before_ns)
        self.assertEqual(observed["fd_after"]["mtime_ns"], before_ns + 1)
        self.assertEqual(observed["link_stat"]["mtime_ns"], before_ns + 1)

    def test_missing_link_retains_consumed_stream_and_descriptor_observations(self):
        original_read = os.read
        removed = False

        def unlinked_read(fd, count):
            nonlocal removed
            raw = original_read(fd, count)
            if raw and not removed:
                removed = True
                self.path.unlink()
            return raw

        with mock.patch.object(controller.os, "read", side_effect=unlinked_read), \
                mock.patch.object(controller.os, "open", wraps=os.open) as opened, \
                self.assertRaisesRegex(controller.SourceIdentityError,
                                       "source identity stat capture failed") as caught:
            controller.read(self.path, self.pin, limit=6)
        self.assertEqual(opened.call_count, 1)
        observed = caught.exception.observed_identity
        self.assertEqual(caught.exception.expected_identity, self.pin)
        self.assertEqual(observed["consumed_bytes"], 6)
        self.assertEqual(observed["consumed_sha256"], self.pin["sha256"])
        self.assertTrue(observed["reached_eof"])
        self.assertEqual(observed["fd_before"]["bytes"], 6)
        self.assertEqual(observed["fd_after"]["bytes"], 6)
        self.assertIsNone(observed["link_stat"])
        self.assertEqual(len(observed["capture_errors"]), 1)
        self.assertEqual(observed["capture_errors"][0]["operation"], "link_stat")
        self.assertEqual(observed["capture_errors"][0]["kind"], "FileNotFoundError")

    def test_terminal_stat_drift_retains_consumed_identity_and_original_authority(self):
        _, original = controller.read(self.path, self.pin, limit=6)
        self.path.chmod(0o444)
        with self.assertRaisesRegex(controller.SourceIdentityError, "terminal file drift") as caught:
            controller.read(self.path, original, limit=6, check_stat=True)
        observed = caught.exception.observed_identity
        self.assertEqual(caught.exception.expected_identity, self.pin)
        self.assertEqual(observed["consumed_sha256"], self.pin["sha256"])
        self.assertEqual(observed["consumed_bytes"], 6)
        self.assertTrue(observed["reached_eof"])
        self.assertEqual(observed["fd_after"]["mode"], "0o444")
        self.assertNotEqual(observed["fd_after"]["mode"], original["mode"])

    def test_initial_and_terminal_failure_formatter_preserves_identity_in_json(self):
        self.path.write_bytes(b"drift!")
        with self.assertRaises(controller.SourceIdentityError) as caught:
            controller.read(self.path, self.pin, limit=6)
        exc = caught.exception
        # These are the two receipt contexts, using their shared failure formatter.
        initial = {**controller.failure(exc), "error": getattr(exc, "error", None)}
        terminal = {"stage": "terminal-identity", "path": str(self.path),
                    **controller.failure(exc, 512)}
        record = json.loads(json.dumps({"failure": initial, "failures": [terminal]}, allow_nan=False))
        for entry in (record["failure"], record["failures"][0]):
            self.assertEqual(entry["kind"], "SourceIdentityError")
            self.assertEqual(entry["expected_identity"], self.pin)
            self.assertEqual(entry["observed_identity"], exc.observed_identity)
            self.assertEqual(entry["observed_identity"]["consumed_sha256"],
                             hashlib.sha256(b"drift!").hexdigest())
            self.assertTrue(entry["observed_identity"]["reached_eof"])
            self.assertIsNone(entry["source_cause"])
        self.assertEqual(record["failures"][0]["stage"], "terminal-identity")

    def test_symlink_leaf_and_ancestor_and_traversal_refuse(self):
        leaf = self.root / "alias"
        leaf.symlink_to(self.path)
        with self.assertRaises(OSError):
            controller.read(leaf)
        directory = self.root / "alias-dir"
        directory.symlink_to(self.root, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, "source ancestor"):
            controller.read(directory / self.path.name)
        with self.assertRaisesRegex(ValueError, "path traversal"):
            controller.read(self.root / ".." / "unrelated")

    def test_exclusive_receipt_write_retains_original_bytes_and_mode(self):
        path = self.root / "record.json"
        controller.write(path, {"status": "refused", "cause": "synthetic-control"})
        before = path.read_bytes()
        digest = hashlib.sha256(before).hexdigest()
        mode = stat.S_IMODE(path.stat().st_mode)
        self.assertEqual(mode, 0o444)
        with self.assertRaises(FileExistsError):
            controller.write(path, {"status": "replacement"})
        self.assertEqual(path.read_bytes(), before)
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), digest)
        self.assertEqual(stat.S_IMODE(path.stat().st_mode), mode)


if __name__ == "__main__":
    unittest.main()
