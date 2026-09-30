"""Independent restoration checks: corrupt data, archive safety and publication."""
import hashlib
import importlib.util
import io
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("restore_data", ROOT / "tools/restore_data.py")
restore_data = importlib.util.module_from_spec(spec)
spec.loader.exec_module(restore_data)


class RestorationChecks(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.data = b"frozen observation\n"
        self.row = {"path": "data/test.dat", "bytes": len(self.data),
                    "sha256": hashlib.sha256(self.data).hexdigest(),
                    "restore": "download", "urls": ["https://example.test/input"]}

    def test_exact_download_and_existing_identity(self):
        with patch.object(restore_data.urllib.request, "urlopen", return_value=io.BytesIO(self.data)):
            self.assertEqual(restore_data.restore(self.row, self.root, None, None), "restored")
        self.assertEqual((self.root / self.row["path"]).read_bytes(), self.data)
        self.assertEqual(restore_data.restore(self.row, self.root, None, None), "already-verified")

    def test_corrupt_same_length_and_oversized_responses_are_rejected(self):
        for content in [b"x" * len(self.data), self.data + b"extra"]:
            with patch.object(restore_data.urllib.request, "urlopen", return_value=io.BytesIO(content)):
                with self.assertRaises(ValueError):
                    restore_data.restore(self.row, self.root, None, None)
            self.assertFalse((self.root / self.row["path"]).exists())
            self.assertEqual(list(self.root.rglob("*.partial")), [])

    def test_existing_changed_input_is_preserved(self):
        dest = self.root / self.row["path"]
        dest.parent.mkdir()
        dest.write_bytes(b"changed")
        with self.assertRaises(ValueError):
            restore_data.restore(self.row, self.root, None, None)
        self.assertEqual(dest.read_bytes(), b"changed")

    def test_archive_restores_only_selected_regular_member(self):
        archive = self.root / "inputs.tar.gz"
        with tarfile.open(archive, "w:gz") as stream:
            item = tarfile.TarInfo(self.row["path"])
            item.size = len(self.data)
            stream.addfile(item, io.BytesIO(self.data))
        row = dict(self.row, restore="frozen-bundle")
        record = {"sha256": restore_data.sha(archive)}
        self.assertEqual(restore_data.restore(row, self.root, archive, record), "restored")
        archive.write_bytes(b"corrupt")
        (self.root / row["path"]).unlink()
        with self.assertRaises(ValueError):
            restore_data.restore(row, self.root, archive, record)

    def test_destination_escape_and_symlink_escape_are_rejected(self):
        for path in ["../outside", "/tmp/outside"]:
            with self.assertRaises(ValueError):
                restore_data.destination(self.root, path)
        (self.root / "data").symlink_to(self.root.parent, target_is_directory=True)
        with self.assertRaises(ValueError):
            restore_data.destination(self.root, "data/outside")


if __name__ == "__main__":
    unittest.main()
