#!/usr/bin/env python3
"""Fail before integration tests when a required runtime/backend is missing."""
import importlib.metadata
import json
from pathlib import Path
import shutil
import subprocess
import sys


def main():
    root = Path(__file__).resolve().parents[1]
    packages = ["numpy", "pandas", "scipy", "statsmodels", "plotly", "matplotlib",
                "scikit-learn", "pyyaml", "pysam", "openpyxl", "snakemake"]
    metadata = {"python": sys.version, "packages": {p: importlib.metadata.version(p) for p in packages}}
    for executable in ("tosa", "Rscript", "stringtie", "featureCounts", "samtools", "snakemake"):
        if not shutil.which(executable):
            raise RuntimeError(f"Missing required executable: {executable}")
    tosa = subprocess.check_output(["tosa", "--version"], text=True).strip()
    if tosa.split()[-1] != "1.0.0":
        raise RuntimeError(f"Expected Tosa 1.0.0, found {tosa}")
    metadata["tosa"] = tosa
    metadata["snakemake"] = subprocess.check_output(["snakemake", "--version"], text=True).strip()
    r_code = '''options(warn=2)
packages <- c("glmmTMB", "TMB", "Matrix", "jsonlite", "DESeq2")
for (p in packages) stopifnot(requireNamespace(p, quietly=TRUE))
cat(jsonlite::toJSON(list(R=R.version.string, packages=setNames(lapply(packages, function(p) as.character(packageVersion(p))), packages)), auto_unbox=TRUE))'''
    metadata["R"] = json.loads(subprocess.check_output(["Rscript", "--vanilla", "-e", r_code], text=True))
    version = (root / "VERSION").read_text().strip()
    for name in ("shiba.py", "scshiba.py"):
        command = shutil.which(name) or str(root / name)
        output = subprocess.check_output([sys.executable, command, "--help"], text=True)
        if version not in output:
            raise RuntimeError(f"{name} does not resolve the checked-out VERSION")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
