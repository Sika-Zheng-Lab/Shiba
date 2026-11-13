#!/usr/bin/env python3
"""Update release labels without rewriting historical version references."""
import argparse
import difflib
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
VERSION_PATTERN = r"v(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)"
IMAGE = r"(docker://naotokubota/shiba:)" + VERSION_PATTERN
TARGETS = {
    "README.md": r"(^# Shiba \()" + VERSION_PATTERN + r"(?=\))",
    "docs/quickstart/diff_splicing_bulk.md": IMAGE,
    "docs/quickstart/diff_splicing_sc.md": IMAGE,
    "example/Shiba_SnakeShiba_config_example.yaml": IMAGE,
    "example/Shiba_SnakeShiba_beta_binomial_config_example.yaml": IMAGE,
    "example/scShiba_SnakeScShiba_config_example.yaml": IMAGE,
    **{f".github/ISSUE_TEMPLATE/{name}.yaml": r"(placeholder: ['\"]e\.g\., )" + VERSION_PATTERN
       for name in ("bug_report", "feature_request", "question")},
}


def proposed_changes(root, version):
    if not re.fullmatch(VERSION_PATTERN, version):
        raise ValueError("Version must be vMAJOR.MINOR.PATCH without leading zeroes")
    changes = {"VERSION": version + "\n"}
    for name, pattern in TARGETS.items():
        original = (root / name).read_text()
        updated, count = re.subn(pattern, lambda match: match[1] + version,
                                 original, flags=re.MULTILINE)
        if count != 1:
            raise ValueError(f"{name}: expected exactly one release label, found {count}")
        changes[name] = updated
    return {name: value for name, value in changes.items()
            if (root / name).read_text() != value}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("version", nargs="?", help="Defaults to the VERSION file")
    parser.add_argument("--root", type=Path, default=ROOT)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", help="Fail if release labels are out of sync")
    mode.add_argument("--dry-run", action="store_true", help="Print the patch without writing")
    args = parser.parse_args()
    try:
        version = args.version or (args.root / "VERSION").read_text().strip()
        changes = proposed_changes(args.root, version)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    for name, updated in changes.items():
        if args.check:
            print(f"Out of sync: {name}")
        elif args.dry_run:
            print("".join(difflib.unified_diff((args.root / name).read_text().splitlines(True),
                                             updated.splitlines(True), fromfile=name, tofile=name)), end="")
        else:
            (args.root / name).write_text(updated)
            print(f"Updated {name}")
    if args.check and changes:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
