"""Challenge transport/provenance behavior without duplicating engine physics."""
import contextlib
import io
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import packet
import run
import restore_inputs

REPOSITORY = packet.ROOT
REVISION = "a" * 40
FAKE_ENGINE = '''#!/usr/bin/env python3
import hashlib,json,pathlib,sys,time
self=pathlib.Path(__file__)
mode=(self.parent/'mode').read_text()
if sys.argv[1]=='describe':
 if mode=='request_changed':
  (self.parent/'experiments/example/request.json').write_text(json.dumps({'schema_version':2,'operation':'background.evaluate','changed':True}))
 if mode=='packet_changed':
  p=self.parent/'experiments/example/experiment.json'; config=json.loads(p.read_text());config['question']='Changed during admission';p.write_text(json.dumps(config))
 print(json.dumps({'product':'Irreducible','build':{'git_head':'a'*40,'git_status':' M bad' if mode=='dirty' else '', 'build_id':'test-build'},'capabilities':[{'id':'background.evaluate','implementation':'implemented'}]}));sys.exit()
if mode=='timeout':time.sleep(5)
if mode=='transport':sys.exit(3)
request=pathlib.Path(sys.argv[2])
receipt={'accepted':mode not in ('failure','qualified'),'accepted_scope':'numerical_contract',
'executable_digest':hashlib.sha256(self.read_bytes()).hexdigest(),'build_id':'test-build','source_revision':'a'*40,'source_status':'',
'input_digest':'wrong' if mode=='mismatch' else hashlib.sha256(request.read_bytes()).hexdigest(),
'numerical':'checks_failed' if mode=='failure' else 'checks_passed','inference':'not_applicable','interpretation':'unqualified'}
print(json.dumps({'receipt':receipt,'result':{}}))
sys.exit(2 if mode=='failure' else 6 if mode=='qualified' else 0)
'''


class ConsumerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        (self.root / "schemas").mkdir()
        shutil.copyfile(REPOSITORY / "schemas/experiment.schema.json", self.root / "schemas/experiment.schema.json")
        self.folder = self.root / "experiments/example"
        self.folder.mkdir(parents=True)
        (self.folder / "README.md").write_text("A test, not scientific evidence.\n")
        self.config = packet.load(REPOSITORY / "experiments/expansion-background/experiment.json")
        self.config["id"] = "example"
        self.config["execution"]["engine_revision"] = REVISION
        self.config["execution"]["timeout_seconds"] = 1
        (self.folder / "request.json").write_text(json.dumps({"schema_version": 2, "operation": "background.evaluate"}))
        self.save()
        self.binary = self.root / "irred-test"
        self.binary.write_text(FAKE_ENGINE)
        self.binary.chmod(0o755)
        self.mode("success")
        self.stack = contextlib.ExitStack()
        self.stack.enter_context(patch.object(packet, "ROOT", self.root))
        self.stack.enter_context(patch.object(run, "ROOT", self.root))
        self.stack.enter_context(patch.object(run, "source_identity", return_value={"revision": REVISION, "status": ""}))
        self.addCleanup(self.temporary.cleanup)
        self.addCleanup(self.stack.close)

    def mode(self, value):
        (self.root / "mode").write_text(value)

    def save(self):
        (self.folder / "experiment.json").write_text(json.dumps(self.config))

    def execute(self, name="attempt"):
        with contextlib.redirect_stdout(io.StringIO()):
            code = run.execute("example", self.binary, name)
        self.record = packet.load(self.root / "results/example" / name / "run.json")
        return code

    def test_success_preserves_identity_and_unqualified_interpretation(self):
        self.assertEqual(self.execute(), 0)
        self.assertEqual(self.record["execution"], "completed")
        self.assertEqual(self.record["qualification"]["interpretation"], "unqualified")
        self.assertEqual(self.record["request_sha256"], packet.sha256(self.folder / "request.json"))

    def test_engine_failures_and_qualification_are_retained(self):
        for mode, code in [("failure", 2), ("qualified", 6), ("transport", 3)]:
            self.mode(mode)
            self.assertEqual(self.execute(mode), code)
            self.assertEqual(self.record["execution"], "failed")
            self.assertEqual(self.record["returncode"], code)

    def test_receipt_identity_mismatch_fails(self):
        self.mode("mismatch")
        self.assertEqual(self.execute(), 1)
        self.assertIn("identity differs", self.record["error"])

    def test_admitted_sources_cannot_change_during_discovery(self):
        for mode in ("request_changed", "packet_changed"):
            self.mode(mode)
            self.assertEqual(self.execute(mode), 1)
            self.assertNotIn("command", self.record)
            self.assertIn("changed after validation", self.record["error"])

    def test_dirty_build_refused_and_recorded(self):
        self.mode("dirty")
        self.assertEqual(self.execute(), 1)
        self.assertNotIn("command", self.record)

    def test_timeout_retains_attempt(self):
        self.mode("timeout")
        self.assertEqual(self.execute(), 124)
        self.assertEqual(self.record["execution"], "timed_out")

    def test_existing_destination_is_not_overwritten(self):
        self.execute()
        before = (self.root / "results/example/attempt/run.json").read_bytes()
        with self.assertRaises(FileExistsError):
            self.execute()
        self.assertEqual(before, (self.root / "results/example/attempt/run.json").read_bytes())

    def test_paths_cannot_escape_by_traversal_or_symlink(self):
        with self.assertRaises(ValueError):
            run.execute("example", self.binary, "../../outside")
        (self.folder / "external").symlink_to(self.root.parent, target_is_directory=True)
        with self.assertRaises(ValueError):
            packet.within(self.folder, "external/anything")

    def test_blocked_experiment_cannot_execute(self):
        self.config.update(status="blocked", execution=None, blockers=["Missing physics"])
        self.save()
        self.assertIsNone(packet.read_packet(self.folder)[1])
        with self.assertRaises(ValueError):
            run.execute("example", self.binary)

    def test_dedicated_blocked_request_pin_checks_exact_bytes_without_execution(self):
        self.config.update(status="blocked", execution=None, blockers=["Missing physics"],
                           request_sha256=packet.sha256(self.folder / "request.json"))
        self.save()
        _, request, identities = packet.read_packet(self.folder)
        self.assertIsNone(request)
        self.assertEqual(identities["controller_request"], self.config["request_sha256"])
        with self.assertRaises(ValueError):
            run.execute("example", self.binary)
        changed = b'{"schema_version":2,"operation":"different.control"}\n'
        (self.folder / "request.json").write_bytes(changed)
        with self.assertRaisesRegex(ValueError, "Dedicated controller request hash differs"):
            packet.read_packet(self.folder)
        self.assertEqual((self.folder / "request.json").read_bytes(), changed)

    def test_dedicated_request_pin_requires_hash_blocked_status_and_confined_file(self):
        self.config["request_sha256"] = packet.sha256(self.folder / "request.json")
        self.save()
        with self.assertRaises(Exception):
            packet.read_packet(self.folder)  # A runnable packet uses its execution binding.
        self.config.update(status="blocked", execution=None, blockers=["Missing physics"])
        for invalid in (True, 12, "A" * 64, "a" * 63):
            self.config["request_sha256"] = invalid
            self.save()
            with self.subTest(invalid=invalid), self.assertRaises(Exception):
                packet.read_packet(self.folder)
        self.config["request_sha256"] = "a" * 64
        self.save()
        (self.folder / "request.json").unlink()
        with self.assertRaises(FileNotFoundError):
            packet.read_packet(self.folder)
        (self.folder / "request.json").symlink_to(self.root.parent / "external-request.json")
        with self.assertRaises(ValueError):
            packet.read_packet(self.folder)

    def test_unknown_fields_duplicates_and_nonfinite_values_rejected(self):
        self.config["shell"] = "anything"
        self.save()
        with self.assertRaises(Exception):
            packet.read_packet(self.folder)
        for text in ['{"x":1,"x":2}', '{"x":NaN}', '{"x":1e999}']:
            path = self.root / "invalid.json"
            path.write_text(text)
            with self.assertRaises(ValueError):
                packet.load(path)

    def test_changed_input_is_preserved_and_refused(self):
        path = self.root / "data/observation.txt"
        path.parent.mkdir()
        path.write_text("original")
        self.config["inputs"] = [{"path": "data/observation.txt", "bytes": 8,
            "sha256": packet.sha256(path), "role": "observation",
            "semantics": "Test bytes, no scientific claim", "source": "synthetic test"}]
        self.save()
        path.write_text("changed!")
        with self.assertRaises(ValueError):
            run.execute("example", self.binary)
        self.assertEqual(path.read_text(), "changed!")

    def test_candidate_identity_and_readiness(self):
        design = {"schema_version": 1, "kind": "candidate_design", "readiness": "ready_for_consumer_review",
                  "blockers": [], "unknowns": [], "review_id": "review",
                  "consumer": {"repository": "sprajs/irreducible", "revision": REVISION, "dirty": False},
                  "executable_request": {"sha256": packet.sha256(self.folder / "request.json"),
                                         "consumer_validation": "passed_at_inspected_revision"}}
        source = self.folder / "candidate.json"
        source.write_text(json.dumps(design))
        self.config["origin"] = {"kind": "prospector_candidate", "repository": "sprajs/prospector",
            "revision": REVISION, "path": "designs/test.json", "snapshot": "candidate.json", "sha256": packet.sha256(source)}
        self.save()
        packet.read_packet(self.folder)
        design["consumer"]["revision"] = "b" * 40
        source.write_text(json.dumps(design))
        self.config["origin"]["sha256"] = packet.sha256(source)
        self.save()
        with self.assertRaises(ValueError):
            packet.read_packet(self.folder)

    def test_restore_preserves_existing_changed_input(self):
        path = self.root / "data/test.txt"
        path.parent.mkdir()
        path.write_text("original")
        row = {"path": "data/test.txt", "bytes": 8, "sha256": packet.sha256(path), "restore": "download", "urls": []}
        self.assertEqual(restore_inputs.restore(row, self.root, None, None), "already-verified")
        path.write_text("changed!")
        with self.assertRaises(ValueError):
            restore_inputs.restore(row, self.root, None, None)
        self.assertEqual(path.read_text(), "changed!")


if __name__ == "__main__":
    unittest.main()
