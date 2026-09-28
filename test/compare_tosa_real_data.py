"""Build a Shiba-format count table and compare Tosa bulk counts with an old run.

This is a migration validation utility. It does not change the production pipeline.
"""

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path


def normalize_chrom(chrom):
    return "chr" + chrom if chrom.isdecimal() or len(chrom) <= 2 else chrom


def parse_coordinate(value):
    chrom, positions = value.split(":", 1)
    start, end = map(int, positions.split("-", 1))
    return normalize_chrom(chrom), start, end


def read_tosa(prefix, boundary_anchor):
    junctions = Counter()
    boundaries = Counter()
    with gzip.open(f"{prefix}_junction.tsv.gz", "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["Junction", "Strand", "Count"]:
            raise ValueError(f"Unexpected junction columns: {reader.fieldnames}")
        for row in reader:
            chrom, start, end = parse_coordinate(row["Junction"])
            junctions[f"{chrom}:{start - 1}-{end + 1}"] += int(row["Count"])

    with gzip.open(f"{prefix}_boundary.tsv.gz", "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["Boundary", "Type", "Strand", "Count"]:
            raise ValueError(f"Unexpected boundary columns: {reader.fieldnames}")
        for row in reader:
            chrom, start, end = parse_coordinate(row["Boundary"])
            if row["Type"] == "5p":
                site = start + boundary_anchor
            elif row["Type"] == "3p":
                site = end - boundary_anchor
            else:
                raise ValueError(f"Unexpected boundary type: {row['Type']}")
            boundaries[f"{chrom}:{site}-{site + 1}"] += int(row["Count"])
    return junctions, boundaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old-bed", required=True)
    parser.add_argument("--out-bed", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--boundary-anchor", type=int, default=1)
    parser.add_argument("--sample", action="append", required=True,
                        help="SAMPLE=PATH_PREFIX; repeat for each sample")
    args = parser.parse_args()
    samples = dict(spec.split("=", 1) for spec in args.sample)
    if len(samples) != len(args.sample):
        parser.error("Sample names must be unique")

    with open(args.old_bed) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames[:4] != ["chr", "start", "end", "ID"] or not all(
            sample in reader.fieldnames[4:] for sample in samples
        ):
            raise ValueError(f"Unexpected old BED columns: {reader.fieldnames}")
        old_rows = {row["ID"]: row for row in reader}

    boundary_ids = {
        key for key, row in old_rows.items()
        if int(row["end"]) - int(row["start"]) == 1
    }
    old_junction_ids = set(old_rows) - boundary_ids
    counts = {}
    report = {"old_rows": len(old_rows), "old_boundary_rows": len(boundary_ids), "samples": {}}
    novel_junction_ids = set()

    for sample, prefix in samples.items():
        junctions, boundaries = read_tosa(prefix, args.boundary_anchor)
        novel_junction_ids.update(junctions.keys() - old_junction_ids)
        counts[sample] = junctions + Counter({key: value for key, value in boundaries.items()
                                               if key in boundary_ids})
        report["samples"][sample] = {}
        for kind, ids in (("junction", old_junction_ids), ("boundary", boundary_ids)):
            exact = sum(int(old_rows[key][sample]) == counts[sample].get(key, 0) for key in ids)
            old_positive = {key for key in ids if int(old_rows[key][sample]) > 0}
            exact_positive = sum(int(old_rows[key][sample]) == counts[sample].get(key, 0)
                                 for key in old_positive)
            report["samples"][sample][kind] = {
                "old_rows": len(ids), "exact_rows": exact,
                "old_positive_rows": len(old_positive), "exact_positive_rows": exact_positive,
                "old_total": sum(int(old_rows[key][sample]) for key in ids),
                "tosa_total_on_old_ids": sum(counts[sample].get(key, 0) for key in ids),
            }
        report["samples"][sample]["new_junction_ids"] = len(junctions.keys() - old_junction_ids)
        report["samples"][sample]["boundary_ids_not_in_old_run"] = len(boundaries.keys() - boundary_ids)

    Path(args.out_bed).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_bed, "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chr", "start", "end", "ID", *samples])
        for key in list(old_rows) + sorted(novel_junction_ids - set(old_rows)):
            chrom, start, end = parse_coordinate(key)
            writer.writerow([chrom, start, end, key, *(counts[sample].get(key, 0) for sample in samples)])
    report["novel_junction_ids_across_samples"] = len(novel_junction_ids)
    with open(args.report, "w") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
