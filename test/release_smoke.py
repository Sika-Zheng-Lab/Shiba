"""Run all four workflows on reproducible BAMs, including assembly and expression.

Synthetic data exercise execution and output contracts, not biological accuracy.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pysam
import yaml
from openpyxl import load_workbook

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "test"))
from make_sc_tosa_fixture import build_fixture, make_read

EVENTS = ("SE", "FIVE", "THREE", "MXE", "RI", "MSE", "AFE", "ALE")


def bulk_fixture(directory, genes=80):
    directory.mkdir(parents=True, exist_ok=True)
    gtf = directory / "annotation.gtf"
    with gtf.open("w") as handle:
        for gene in range(genes):
            base = 3000 * gene + 1000
            strand = "+" if gene % 2 == 0 else "-"
            for isoform, exons in enumerate([
                [(1, 200), (500, 700), (1000, 1200)],
                [(1, 200), (1000, 1200)],
                [(1, 700), (1000, 1200)],
            ]):
                for start, end in exons:
                    attributes = f'gene_id "G{gene}"; gene_name "Gene{gene}"; transcript_id "G{gene}T{isoform}";'
                    handle.write(f"chr1\tsynthetic\texon\t{base+start}\t{base+end}\t.\t{strand}\t.\t{attributes}\n")
    rng = np.random.default_rng(20261113)
    expression = rng.uniform(25, 100, genes)
    metadata, experiment = [], []
    for genotype in ("WT", "KO"):
        for treatment in ("Control", "Drug"):
            for replicate in range(2):
                sample = f"sample_{genotype}_{treatment}_{replicate}"
                group = f"{genotype}_{treatment}"
                metadata.append([sample, genotype, treatment, f"B{replicate}"])
                bam = directory / f"{sample}.bam"
                header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": genes*3000+5000}]}
                reads = []
                for gene in range(genes):
                    base = gene * 3000 + 1000
                    interaction = genotype == "KO" and treatment == "Drug"
                    mean = expression[gene] * (2 if interaction and gene < 10 else 1)
                    total = max(16, int(rng.negative_binomial(8, 8 / (8 + mean))))
                    probability = .7 if interaction and gene < 20 else .35
                    inc = int(rng.binomial(total, rng.beta(probability*30, (1-probability)*30)))
                    retained = max(2, int(rng.binomial(total, .1)))
                    patterns = [(base+150, "50M299N50M", inc), (base+650, "50M299N50M", inc),
                                (base+150, "50M799N50M", total-inc),
                                (base+150, "100M", retained), (base+449, "100M", retained)]
                    for offset in (0, 499, 999):
                        for _ in range(total):
                            patterns.append((base+offset+int(rng.integers(0, 100)), "100M", 1))
                    for start, cigar, count in patterns:
                        for _ in range(count):
                            read = make_read(f"{sample}_{gene}_{len(reads)}", start, cigar, "CELL", str(len(reads)))
                            read.set_tag("XS", "+" if gene % 2 == 0 else "-")
                            reads.append(read)
                with pysam.AlignmentFile(bam, "wb", header=header) as handle:
                    for read in sorted(reads, key=lambda read: read.reference_start):
                        handle.write(read)
                pysam.index(str(bam))
                experiment.append([sample, str(bam), group])
    pd.DataFrame(metadata, columns=["sample", "genotype", "treatment", "batch"]).to_csv(directory / "samples.tsv", sep="\t", index=False)
    pd.DataFrame(experiment, columns=["sample", "bam", "group"]).to_csv(directory / "experiment.tsv", sep="\t", index=False)
    return dict(gtf=str(gtf), experiment_table=str(directory / "experiment.tsv"),
                container="unused-local-runtime", unannotated=True, strand="unstranded",
                minimum_reads=1, reference_group="WT_Control", alternative_group="KO_Drug",
                ttest=False, excel=True, individual_psi=True,
                sample_metadata=str(directory / "samples.tsv"))


def execute(command, root, name):
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    with (root / f"{name}.log").open("w") as log:
        result = subprocess.run(list(map(str, command)), cwd=root, env=env, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"{name} failed:\n" + (root / f"{name}.log").read_text()[-7000:])


def check_outputs(path, bulk, regression):
    splicing = path / "results/splicing" if bulk else path / "results"
    for event in EVENTS:
        frame = pd.read_csv(splicing / f"PSI_{event}.txt", sep="\t")
        assert "event_id" in frame.columns, event
    assert (path / "report.json").is_file()
    assert (path / "plots/summary.html").stat().st_size > 1000
    if regression:
        table = pd.read_csv(splicing / "event_statistics.tsv", sep="\t")
        assert not table.empty
        assert set(table.event_type) <= set(EVENTS)
        assert (splicing / "design.json").exists()
        assert (splicing / "factorial_results.xlsx").exists()
    else:
        book = load_workbook(splicing / "results.xlsx", read_only=True)
        assert book.sheetnames == list(EVENTS)
        book.close()
    if bulk:
        for name in ("counts.txt", "TPM.txt", "CPM.txt", "DEG.txt"):
            frame = pd.read_csv(path / "results/expression" / name, sep="\t")
            assert not frame.empty and len(frame.columns) > 1, name
        for name in ("psi_pca.tsv", "tpm_pca.tsv"):
            assert (path / "results/pca" / name).stat().st_size > 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--pipelines", nargs="+", default=["shiba-fisher", "shiba-factorial", "snakeshiba-fisher", "snakeshiba-factorial", "scshiba", "snakescshiba"])
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    bulk = bulk_fixture(root / "bulk-input")
    sc_input = build_fixture(root / "sc-input")
    sc = yaml.safe_load((sc_input / "config.yaml").read_text())
    sc["excel"] = True
    completed = []
    for name in args.pipelines:
        is_bulk = "factorial" in name or "fisher" in name
        regression = "factorial" in name
        config = dict(bulk if is_bulk else sc)
        workdir = root / name
        config["workdir"] = str(workdir)
        if regression:
            config.update(stat_method="beta_binomial", formula="genotype * treatment + batch",
                          categorical=["genotype", "treatment", "batch"],
                          reference_levels={"genotype": "WT", "treatment": "Control"},
                          coef=["genotypeKO:treatmentDrug"],
                          expression_reference_group="WT_Control", expression_alternative_group="KO_Drug")
        path = root / f"{name}.yaml"
        path.write_text(yaml.safe_dump(config))
        if name.startswith("snake"):
            command = ["snakemake", "-s", ROOT / ("snakeshiba.smk" if is_bulk else "snakescshiba.smk"),
                       "--configfile", path, "--cores", "2"]
        else:
            command = [sys.executable, ROOT / ("shiba.py" if is_bulk else "scshiba.py"), path, "-p", "2"]
        print(f"Running {name}", flush=True)
        execute(command, root, name)
        check_outputs(workdir, is_bulk, regression)
        completed.append(name)
        (root / "validation.json").write_text(json.dumps({"completed": completed, "fixture": "synthetic", "version": (ROOT / "VERSION").read_text().strip()}, indent=2))
    if all(name in completed for name in ("scshiba", "snakescshiba")):
        first = pd.read_csv(root / "scshiba/junctions/junctions.bed", sep="\t").sort_values("ID").reset_index(drop=True)
        second = pd.read_csv(root / "snakescshiba/junctions/junctions.bed", sep="\t").sort_values("ID").reset_index(drop=True)
        pd.testing.assert_frame_equal(first, second)
    print(json.dumps({"completed": completed}))


if __name__ == "__main__":
    main()
