"""Release labels must not change migration guidance or historical dates."""
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.update_version import TARGETS, proposed_changes


class TestVersionUpdate(unittest.TestCase):
    def test_updates_only_declared_release_labels(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ["VERSION", *TARGETS, "CHANGELOG.md", "docs/installation.md"]:
                target = root / name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(ROOT / name, target)
            before = {name: (root / name).read_bytes() for name in ("CHANGELOG.md", "docs/installation.md")}
            cmd = [sys.executable, str(ROOT / "scripts/update_version.py"), "v1.2.3", "--root", str(root)]
            subprocess.run(cmd, check=True, capture_output=True)
            self.assertEqual((root / "VERSION").read_text(), "v1.2.3\n")
            self.assertIn("# Shiba (v1.2.3)", (root / "README.md").read_text())
            self.assertIn("v0.8.2 package/image", (root / "README.md").read_text())
            self.assertIn("shiba:v1.2.3", (root / "example/Shiba_SnakeShiba_beta_binomial_config_example.yaml").read_text())
            for name, original in before.items():
                self.assertEqual((root / name).read_bytes(), original)
            self.assertEqual(proposed_changes(root, "v1.2.3"), {})
            subprocess.run(cmd + ["--check"], check=True, capture_output=True)

    def test_invalid_versions_are_rejected_before_writing(self):
        for version in ("1.0.0", "v01.0.0", "v1.0.0\n", "v1.0.0; echo BAD", "v1.0.0-rc1"):
            with self.subTest(version=version), self.assertRaises(ValueError):
                proposed_changes(ROOT, version)

    def test_missing_label_fails_before_writing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "VERSION").write_text("v0.8.2\n")
            (root / "README.md").write_text("Unexpected heading\n")
            with self.assertRaisesRegex(ValueError, "exactly one"):
                proposed_changes(root, "v1.0.0")
            self.assertEqual((root / "VERSION").read_text(), "v0.8.2\n")
