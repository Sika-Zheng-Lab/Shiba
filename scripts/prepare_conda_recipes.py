#!/usr/bin/env python3
"""Generate Shiba/MameShiba recipes from a verified release source manifest.

Without --manifest, generate local-only candidate recipes from this checkout.
This does not change VERSION, submit recipes, or publish packages.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil

import yaml

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--output", type=Path, default=ROOT / "release-artifacts/recipes")
    args = parser.parse_args()
    version = (ROOT / "VERSION").read_text().strip()
    if args.manifest:
        manifest = json.loads(args.manifest.read_text())
        if manifest["version"] != version:
            parser.error("Manifest version must match VERSION")
        archive = args.manifest.parent / manifest["source"]
        if hashlib.sha256(archive.read_bytes()).hexdigest() != manifest["sha256"]:
            parser.error("Source archive checksum does not match the manifest")
        source = {"url": f"https://github.com/Sika-Zheng-Lab/Shiba/releases/download/{version}/{archive.name}",
                  "sha256": manifest["sha256"]}
    else:
        # Stage only distributable files: source.path must not copy .git,
        # virtual environments, or large validation data into conda-build.
        staging = args.output.resolve() / "source"
        staging.mkdir(parents=True, exist_ok=True)
        for name in ("src", "scripts", "example"):
            shutil.copytree(ROOT / name, staging / name, dirs_exist_ok=True,
                            ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        for name in ("shiba.py", "scshiba.py", "snakeshiba.smk", "snakescshiba.smk", "VERSION", "LICENSE", "environment.yml"):
            shutil.copy2(ROOT / name, staging / name)
        source = {"path": str(staging)}
    if not re.fullmatch(r"v\d+\.\d+\.\d+", version):
        parser.error("Invalid VERSION")
    dependencies = yaml.safe_load((ROOT / "environment.yml").read_text())["dependencies"]
    # Conda recipes use spaces rather than environment.yml's equals syntax.
    dependencies = [d.replace("=", " ", 1) for d in dependencies if isinstance(d, str) and d not in ("git", "pip")]
    for name in ("shiba", "mameshiba"):
        runtime = dependencies.copy()
        if name == "mameshiba":
            runtime = [d for d in runtime if not d.startswith(("bioconductor-deseq2", "r-data.table", "r-locfit", "subread", "matplotlib-base", "seaborn", "scikit-learn", "plotly"))]
        recipe = {
            "package": {"name": name, "version": version[1:]},
            "source": source,
            "build": {"number": 0, "noarch": "generic", "script": [
                'mkdir -p "$PREFIX/bin" "$PREFIX/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM"',
                'cp -R src scripts example shiba.py scshiba.py snakeshiba.smk snakescshiba.smk VERSION LICENSE environment.yml "$PREFIX/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM/"',
                'chmod 0755 "$PREFIX/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM/"{shiba,scshiba}.py',
                'ln -s "$PREFIX/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM/shiba.py" "$PREFIX/bin/shiba.py"',
                'ln -s "$PREFIX/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM/scshiba.py" "$PREFIX/bin/scshiba.py"',
            ]},
            "requirements": {"host": ["python 3.12"], "run": runtime},
            "test": {"commands": [
                f'shiba.py --help | grep "{version}"',
                f'scshiba.py --help | grep "{version}"',
                'test -s "$PREFIX/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM/src/beta_binomial_glm.R"',
                'test -d "$PREFIX/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM/src/templates"',
                'tosa --version | grep "1.0.0"',
                "Rscript --vanilla -e 'options(warn=2); stopifnot(requireNamespace(\"glmmTMB\"), requireNamespace(\"jsonlite\"))'",
            ]},
            "about": {"home": "https://github.com/Sika-Zheng-Lab/Shiba", "license": "MIT",
                      "license_file": "LICENSE", "summary": "Differential RNA splicing analysis with Tosa and beta-binomial regression",
                      "doc_url": "https://sika-zheng-lab.github.io/Shiba"},
            "extra": {"recipe-maintainers": ["NaotoKubota"]},
        }
        directory = args.output / name
        directory.mkdir(parents=True, exist_ok=True)
        header = "# Local candidate only; regenerate with --manifest before submission.\n" if not args.manifest else "# Generated from the verified release source archive.\n"
        (directory / "meta.yaml").write_text(header + yaml.safe_dump(recipe, sort_keys=False, width=120))
        print(directory / "meta.yaml")


if __name__ == "__main__":
    main()
