"""Create a documented regional validation subset from existing public BAMs.

Never changes the original alignments. Region coordinates are 0-based, half-open.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import pysam
import yaml


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--groups", nargs=2, required=True)
    parser.add_argument("--samples-per-group", type=int, default=2)
    parser.add_argument("--chromosome", default="1")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--end", type=int, default=10000000)
    parser.add_argument("--technology", choices=["short", "long"], default="short")
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    with args.experiment.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    selected = []
    for group in args.groups:
        selected += [row for row in rows if row["group"] == group][:args.samples_per_group]
    if len(selected) != 2 * args.samples_per_group:
        parser.error("Not enough samples for both groups")
    manifest = dict(source_experiment=str(args.experiment), source_gtf=str(args.gtf),
                    region=[args.chromosome, args.start, args.end], technology=args.technology, samples=[])
    output_rows = []
    for row in selected:
        source = Path(row["bam"])
        target = root / f"{row['sample']}.bam"
        with pysam.AlignmentFile(source, "rb") as original, pysam.AlignmentFile(target, "wb", template=original) as output:
            count = 0
            for read in original.fetch(args.chromosome, args.start, args.end):
                output.write(read)
                count += 1
        if not count:
            raise RuntimeError(f"Empty subset: {source}")
        pysam.index(str(target))
        manifest["samples"].append(dict(sample=row["sample"], source_bam=str(source), reads=count,
                                        subset_sha256=hashlib.sha256(target.read_bytes()).hexdigest()))
        output_rows.append([row["sample"], str(target), row["group"], args.technology])
    gtf = root / "annotation.gtf"
    with args.gtf.open() as original, gtf.open("w") as output:
        for line in original:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if fields[0] == args.chromosome and int(fields[3]) > args.start and int(fields[4]) <= args.end:
                output.write(line)
    with (root / "experiment.tsv").open("w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "bam", "group", "technology"])
        writer.writerows(output_rows)
    config = dict(workdir=str(root / "fisher"), gtf=str(gtf), experiment_table=str(root / "experiment.tsv"),
                  container="unused-local-runtime", unannotated=True, strand="unstranded", minimum_reads=10,
                  reference_group=args.groups[0], alternative_group=args.groups[1],
                  ttest=False, individual_psi=True, excel=True)
    (root / "fisher.yaml").write_text(yaml.safe_dump(config))
    with (root / "samples.tsv").open("w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "group"])
        writer.writerows([[r[0], r[2]] for r in output_rows])
    config.update(workdir=str(root / "regression"), stat_method="beta_binomial", formula="group",
                  sample_metadata=str(root / "samples.tsv"), categorical=["group"],
                  reference_levels={"group": args.groups[0]}, coef=["group" + args.groups[1]],
                  expression_reference_group=args.groups[0], expression_alternative_group=args.groups[1])
    (root / "regression.yaml").write_text(yaml.safe_dump(config))
    (root / "input-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
