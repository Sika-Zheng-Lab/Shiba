"""Count single-cell junctions and RI boundaries from BAM/CRAM with Tosa."""

import argparse
import csv
import gzip
import logging
import os
import re
import subprocess
import tempfile
from collections import Counter

from lib.general import normalize_tosa_strand


logger = logging.getLogger(__name__)
SAMPLE_NAME = re.compile(r"^[A-Za-z0-9._-]+$")


def parse_coordinate(value):
    chrom, positions = value.rsplit(":", 1)
    start, end = map(int, positions.split("-", 1))
    return chrom, start, end


def read_experiment_table(path):
    samples = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"sample", "alignment", "barcode"}.issubset(reader.fieldnames):
            raise ValueError("Single-cell experiment table needs sample, alignment, and barcode columns")
        for row in reader:
            sample = row["sample"]
            if not sample or not SAMPLE_NAME.fullmatch(sample) or sample in samples:
                raise ValueError(f"Invalid or duplicate sample name: {sample!r}")
            if not row["alignment"] or not row["barcode"]:
                raise ValueError(f"Missing alignment or barcode file for {sample}")
            samples[sample] = row
    if not samples:
        raise ValueError("Single-cell experiment table has no samples")
    return samples


def read_barcode_groups(path):
    groups = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"barcode", "group"}.issubset(reader.fieldnames):
            raise ValueError(f"Barcode table needs barcode and group columns: {path}")
        for row in reader:
            barcode, group = row["barcode"], row["group"]
            if not barcode or not group:
                raise ValueError(f"Empty barcode or group in {path}")
            if barcode in groups and groups[barcode] != group:
                raise ValueError(f"Barcode {barcode} has conflicting groups in {path}")
            groups[barcode] = group
    if not groups:
        raise ValueError(f"Barcode table has no cells: {path}")
    return groups


def ri_boundary_ids(path):
    ids = set()
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "intron_a" not in reader.fieldnames:
            raise ValueError("RI event file needs an intron_a column")
        for row in reader:
            chrom, start, end = parse_coordinate(row["intron_a"])
            ids.add(f"{chrom}:{start}-{start + 1}")
            ids.add(f"{chrom}:{end - 1}-{end}")
    return ids


def tosa_command(alignment, gtf, prefix, whitelist, threads, anchor,
                 boundary_anchor, min_intron, max_intron, strand):
    command = ["tosa", "single", "-c", whitelist, "-g", gtf,
               "-a", str(anchor), "-b", str(boundary_anchor),
               "-m", str(min_intron), "-M", str(max_intron), "-p", str(threads)]
    strand = normalize_tosa_strand(strand)
    if strand is not None:
        command.extend(["-s", strand])
    command.extend([alignment, prefix])
    return command


def run_tosa(alignment, barcode_table, gtf, prefix, threads=1, anchor=8,
             boundary_anchor=1, min_intron=20, max_intron=500000,
             strand="unstranded"):
    groups = read_barcode_groups(barcode_table)
    directory = os.path.dirname(os.path.abspath(prefix))
    os.makedirs(directory, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", prefix="barcodes_", suffix=".tsv",
                                     dir=directory, delete=False) as handle:
        whitelist = handle.name
        handle.writelines(f"{barcode}\n" for barcode in sorted(groups))
    try:
        command = tosa_command(alignment, gtf, prefix, whitelist, threads, anchor,
                               boundary_anchor, min_intron, max_intron, strand)
        logger.info("Running Tosa single on %s", alignment)
        subprocess.run(command, check=True)
    except FileNotFoundError as error:
        raise RuntimeError("Tosa executable was not found; install tosa and add it to PATH") from error
    finally:
        os.unlink(whitelist)


def add_tosa_counts(prefix, barcode_groups, boundary_anchor, wanted_boundaries, counts):
    junction_file = f"{prefix}_junction_barcodes.tsv.gz"
    with gzip.open(junction_file, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["Feature", "Strand", "Barcode", "Count"]:
            raise ValueError(f"Unexpected Tosa junction columns in {junction_file}: {reader.fieldnames}")
        for row in reader:
            group = barcode_groups.get(row["Barcode"])
            if group is None:
                continue
            chrom, start, end = parse_coordinate(row["Feature"])
            key = f"{chrom}:{start - 1}-{end + 1}"
            counts[group][key] += int(row["Count"])

    boundary_file = f"{prefix}_boundary_barcodes_detail.tsv.gz"
    with gzip.open(boundary_file, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["Boundary", "Type", "Strand", "Barcode", "Count"]:
            raise ValueError(f"Unexpected Tosa boundary columns in {boundary_file}: {reader.fieldnames}")
        for row in reader:
            group = barcode_groups.get(row["Barcode"])
            if group is None:
                continue
            chrom, start, end = parse_coordinate(row["Boundary"])
            if row["Type"] == "5p":
                site = start + boundary_anchor
            elif row["Type"] == "3p":
                site = end - boundary_anchor
            else:
                raise ValueError(f"Unexpected Tosa boundary type: {row['Type']}")
            key = f"{chrom}:{site}-{site + 1}"
            if key in wanted_boundaries:
                counts[group][key] += int(row["Count"])


def merge_tosa_samples(samples, prefixes, ri_event, output, boundary_anchor=1):
    wanted_boundaries = ri_boundary_ids(ri_event)
    counts = {}
    for sample, row in samples.items():
        barcode_groups = read_barcode_groups(row["barcode"])
        for group in barcode_groups.values():
            counts.setdefault(group, Counter())
        add_tosa_counts(prefixes[sample], barcode_groups, boundary_anchor,
                        wanted_boundaries, counts)
    groups = sorted(counts)
    ids = set(wanted_boundaries)
    for group_counts in counts.values():
        ids.update(group_counts)
    rows = []
    for key in ids:
        chrom, start, end = parse_coordinate(key)
        rows.append((chrom, start, end, key))
    rows.sort(key=lambda row: (row[0], row[1], row[2], row[3]))
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    with open(output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chr", "start", "end", "ID", *groups])
        for chrom, start, end, key in rows:
            writer.writerow([chrom, start, end, key,
                             *(counts[group].get(key, 0) for group in groups)])
    logger.info("Wrote %d junction and RI boundary rows to %s", len(rows), output)


def add_tosa_options(parser, require_gtf):
    parser.add_argument("-g", "--gtf", required=require_gtf)
    parser.add_argument("-p", "--processors", type=int, default=1)
    parser.add_argument("-a", "--anchor", type=int, default=8)
    parser.add_argument("-b", "--boundary-anchor", type=int, default=1)
    parser.add_argument("-m", "--min-intron", type=int, default=20)
    parser.add_argument("-M", "--max-intron", type=int, default=500000)
    parser.add_argument("-s", "--strand", default="unstranded")
    parser.add_argument("-v", "--verbose", action="store_true")


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", help="sample/alignment/barcode experiment table")
    parser.add_argument("-r", "--ri-event", help="EVENT_RI.txt")
    parser.add_argument("-o", "--output", help="Shiba junctions.bed")
    add_tosa_options(parser, False)
    commands = parser.add_subparsers(dest="subcommand")
    run = commands.add_parser("run", help="Run Tosa on one alignment")
    run.add_argument("--alignment", required=True)
    run.add_argument("--barcode", required=True)
    run.add_argument("--prefix", required=True)
    add_tosa_options(run, True)
    merge = commands.add_parser("merge", help="Combine per-sample Tosa results")
    merge.add_argument("-i", "--input", required=True)
    merge.add_argument("-d", "--directory", required=True)
    merge.add_argument("-r", "--ri-event", required=True)
    merge.add_argument("-o", "--output", required=True)
    merge.add_argument("-b", "--boundary-anchor", type=int, default=1)
    merge.add_argument("-v", "--verbose", action="store_true")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="[%(asctime)s] %(levelname)7s %(message)s")
    if args.subcommand == "run":
        run_tosa(args.alignment, args.barcode, args.gtf, args.prefix, args.processors,
                 args.anchor, args.boundary_anchor, args.min_intron, args.max_intron, args.strand)
    elif args.subcommand == "merge":
        samples = read_experiment_table(args.input)
        prefixes = {sample: os.path.join(args.directory, sample) for sample in samples}
        merge_tosa_samples(samples, prefixes, args.ri_event, args.output, args.boundary_anchor)
    else:
        if not all((args.input, args.ri_event, args.output, args.gtf)):
            parser.error("-i, -r, -o, and -g are required for the all-in-one run")
        samples = read_experiment_table(args.input)
        output_dir = os.path.dirname(os.path.abspath(args.output))
        os.makedirs(output_dir, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="tosa_sc_", dir=output_dir) as tmp:
            prefixes = {}
            for sample, row in samples.items():
                prefixes[sample] = os.path.join(tmp, sample)
                run_tosa(row["alignment"], row["barcode"], args.gtf, prefixes[sample],
                         args.processors, args.anchor, args.boundary_anchor,
                         args.min_intron, args.max_intron, args.strand)
            merge_tosa_samples(samples, prefixes, args.ri_event, args.output, args.boundary_anchor)


if __name__ == "__main__":
    main()
