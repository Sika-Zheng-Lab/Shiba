#!/usr/bin/env python3
"""Validate a release commit and prepare deterministic source/notes artifacts.

This command never creates tags, publishes releases, or updates containers.
"""
import argparse
from datetime import date
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from update_version import ROOT, VERSION_PATTERN, proposed_changes


def release_notes(changelog, version):
    entries = list(re.finditer(r"^## \[(v[^\]]+)\] - (.+)$", changelog, re.MULTILINE))
    if not entries or entries[0][1] != version:
        raise ValueError("VERSION must match the first CHANGELOG release entry")
    date.fromisoformat(entries[0][2])
    end = entries[1].start() if len(entries) > 1 else len(changelog)
    return changelog[entries[0].start():end].strip() + "\n"


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


def validate_commit(root, expected_sha, version):
    if not re.fullmatch(r"[0-9a-f]{40}", expected_sha):
        raise ValueError("Expected SHA must be a full commit hash")
    if git(root, "rev-parse", "HEAD") != expected_sha:
        raise ValueError("Checked-out commit differs from the requested release SHA")
    if git(root, "rev-parse", "origin/main") != expected_sha:
        raise ValueError("Candidate is no longer the tip of origin/main; do not move latest backwards")
    exists = subprocess.run(["git", "-C", str(root), "show-ref", "--verify", "--quiet", f"refs/tags/{version}"])
    if exists.returncode == 0:
        if git(root, "rev-parse", f"refs/tags/{version}^{{commit}}") != expected_sha:
            raise ValueError(f"{version} already points to a different commit")
    elif exists.returncode != 1:
        raise RuntimeError("Unable to inspect release tag")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--expected-sha", required=True)
    parser.add_argument("--output", type=Path, default=Path("release-artifacts"))
    args = parser.parse_args()
    version = (args.root / "VERSION").read_text().strip()
    if not re.fullmatch(VERSION_PATTERN, version):
        parser.error("Invalid VERSION")
    try:
        notes = release_notes((args.root / "CHANGELOG.md").read_text(), version)
        if proposed_changes(args.root, version):
            raise ValueError("Release labels are out of sync; run scripts/update_version.py")
        validate_commit(args.root, args.expected_sha, version)
    except ValueError as error:
        parser.error(str(error))
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "release-notes.md").write_text(notes)
    archive = args.output / f"shiba-{version[1:]}.tar.gz"
    with archive.open("wb") as handle, gzip.GzipFile(filename="", mode="wb", fileobj=handle, mtime=0) as compressed:
        process = subprocess.Popen(["git", "-C", str(args.root), "archive", "--format=tar",
                                    f"--prefix=Shiba-{version[1:]}/", args.expected_sha], stdout=subprocess.PIPE)
        with process.stdout:
            shutil.copyfileobj(process.stdout, compressed)
        if process.wait() != 0:
            raise RuntimeError("git archive failed")
    manifest = dict(version=version, commit=args.expected_sha, platform="linux/amd64",
                    source=archive.name, sha256=hashlib.sha256(archive.read_bytes()).hexdigest())
    (args.output / "source-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest))


if __name__ == "__main__":
    main()
