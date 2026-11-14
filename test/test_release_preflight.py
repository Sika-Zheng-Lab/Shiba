from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.prepare_release import release_notes, validate_commit


class TestReleasePreflight(unittest.TestCase):
    def test_notes_require_version_and_real_date(self):
        text = "# Change log\n\n## [v1.0.0] - 2026-11-13\n\nBreaking changes.\n\n## [v0.8.2] - 2026-02-23\nOld.\n"
        self.assertIn("Breaking changes.", release_notes(text, "v1.0.0"))
        self.assertNotIn("Old.", release_notes(text, "v1.0.0"))
        for candidate, version in [(text, "v0.8.2"), (text.replace("2026-11-13", "2026-XX-XX"), "v1.0.0")]:
            with self.assertRaises(ValueError):
                release_notes(candidate, version)

    def test_same_tag_retry_allowed_but_mismatched_tag_and_stale_main_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            def git(*args):
                return subprocess.check_output(["git", "-C", tmp, *args], text=True, stderr=subprocess.DEVNULL).strip()
            git("init")
            git("-c", "user.name=test", "-c", "user.email=test@example.com", "commit", "--allow-empty", "-m", "first")
            first = git("rev-parse", "HEAD")
            git("update-ref", "refs/remotes/origin/main", first)
            validate_commit(Path(tmp), first, "v1.0.0")
            git("tag", "v1.0.0")
            validate_commit(Path(tmp), first, "v1.0.0")
            git("-c", "user.name=test", "-c", "user.email=test@example.com", "commit", "--allow-empty", "-m", "second")
            second = git("rev-parse", "HEAD")
            with self.assertRaisesRegex(ValueError, "tip"):
                validate_commit(Path(tmp), second, "v1.1.0")
            git("update-ref", "refs/remotes/origin/main", second)
            with self.assertRaisesRegex(ValueError, "different commit"):
                validate_commit(Path(tmp), second, "v1.0.0")
            with self.assertRaisesRegex(ValueError, "Checked-out"):
                validate_commit(Path(tmp), first, "v1.0.0")
