"""Keep public storage limits separate from scientific and local evidence."""
import contextlib
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import check_repository


class StorageTests(unittest.TestCase):
    def check_tree(self, files):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name, size in files.items():
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                with path.open("wb") as stream:
                    stream.truncate(size)
            inventory = ("\0".join(files) + "\0").encode()
            with patch.object(check_repository, "ROOT", root), \
                    patch.object(check_repository, "validate_all", return_value=0), \
                    patch.object(check_repository.subprocess, "check_output", return_value=inventory), \
                    contextlib.redirect_stdout(io.StringIO()):
                check_repository.check()

    def test_two_mib_boundary_and_excess(self):
        self.check_tree({"source.txt": 2 * 1024 * 1024})
        with self.assertRaisesRegex(ValueError, "exceeds 2 MiB"):
            self.check_tree({"source.txt": 2 * 1024 * 1024 + 1})

    def test_local_and_generated_files_remain_excluded(self):
        for name in ("results/receipt.json", "data/input.txt", ".work/source.py", "figure.png"):
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, "Generated/local artifact"):
                self.check_tree({name: 1})

    def test_eight_source_files_per_packet_remain_the_limit(self):
        files = {f"experiments/example/source{i}.txt": 1 for i in range(8)}
        self.check_tree(files)
        files["experiments/example/ninth.txt"] = 1
        with self.assertRaisesRegex(ValueError, "exceeds eight"):
            self.check_tree(files)


if __name__ == "__main__":
    unittest.main()
