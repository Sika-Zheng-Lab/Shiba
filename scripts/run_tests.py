#!/usr/bin/env python3
"""Run the complete test suite and reject skipped tests in release environments."""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
suite = unittest.defaultTestLoader.discover(str(ROOT / "test"), top_level_dir=str(ROOT))
result = unittest.TextTestRunner(verbosity=2).run(suite)
if result.skipped:
    print("Release validation requires zero skips:", result.skipped, file=sys.stderr)
raise SystemExit(0 if result.wasSuccessful() and not result.skipped else 1)
