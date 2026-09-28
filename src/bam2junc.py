"""Count bulk splice junctions and RI boundaries with Tosa for Shiba."""

import argparse
import csv
import gzip
import logging
import os
import subprocess
import tempfile
from collections import Counter

from lib.general import normalize_tosa_strand


logger = logging.getLogger(__name__)


def normalize_chrom(chrom):
    """Match the chromosome spelling used by gtf2event and the old merger."""
    return "chr" + chrom if chrom.isdecimal() or len(chrom) <= 2 else chrom


def parse_coordinate(value):
    chrom, positions = value.rsplit(":", 1)
    start, end = map(int, positions.split("-", 1))
    return normalize_chrom(chrom), start, end


def read_experiment_table(path):
    samples = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"sample", "bam"}.issubset(reader.fieldnames):
            raise ValueError("Experiment table must contain sample and bam columns")
        for row in reader:
            sample = row["sample"]
            if not sample or sample in samples:
                raise ValueError(f"Missing or duplicate sample name: {sample!r}")
            samples[sample] = row["bam"]
    if not samples:
        raise ValueError("Experiment table has no samples")
    return samples


def ri_boundary_ids(path):
    """Return Shiba's two boundary IDs for each RI event."""
    ids = set()
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "intron_a" not in reader.fieldnames:
            raise ValueError("RI event file must contain an intron_a column")
        for row in reader:
            chrom, start, end = parse_coordinate(row["intron_a"])
            ids.add(f"{chrom}:{start}-{start + 1}")
            ids.add(f"{chrom}:{end - 1}-{end}")
    return ids


def tosa_command(bam, gtf, prefix, threads, anchor, boundary_anchor,
                 min_intron, max_intron, strand):
    command = [
        "tosa", "bulk", "-a", str(anchor), "-b", str(boundary_anchor),
        "-m", str(min_intron), "-M", str(max_intron), "-p", str(threads),
        "-g", gtf,
    ]
    strand = normalize_tosa_strand(strand)
    if strand is not None:
        command.extend(["-s", strand])
    command.extend([bam, prefix])
    return command


def run_tosa(bam, gtf, prefix, threads=1, anchor=8, boundary_anchor=1,
             min_intron=70, max_intron=500000, strand="XS", log_file=None):
    os.makedirs(os.path.dirname(os.path.abspath(prefix)), exist_ok=True)
    command = tosa_command(bam, gtf, prefix, threads, anchor, boundary_anchor,
                           min_intron, max_intron, strand)
    logger.info("Running Tosa on %s", bam)
    logger.debug("Command: %s", command)
    try:
        if log_file:
            os.makedirs(os.path.dirname(os.path.abspath(log_file)), exist_ok=True)
            with open(log_file, "w") as handle:
                subprocess.run(command, check=True, stdout=handle, stderr=subprocess.STDOUT)
        else:
            subprocess.run(command, check=True)
    except FileNotFoundError as error:
        raise RuntimeError("Tosa executable was not found; install tosa and add it to PATH") from error
    except subprocess.CalledProcessError as error:
        raise RuntimeError(f"Tosa failed for {bam} (exit code {error.returncode}); see {log_file}") from error


def read_tosa_counts(prefix, boundary_anchor, wanted_boundaries):
    """Convert Tosa coordinates to Shiba IDs and aggregate strand rows."""
    counts = Counter()
    with gzip.open(f"{prefix}_junction.tsv.gz", "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["Junction", "Strand", "Count"]:
            raise ValueError(f"Unexpected Tosa junction columns: {reader.fieldnames}")
        for row in reader:
            chrom, start, end = parse_coordinate(row["Junction"])
            counts[f"{chrom}:{start - 1}-{end + 1}"] += int(row["Count"])

    with gzip.open(f"{prefix}_boundary.tsv.gz", "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["Boundary", "Type", "Strand", "Count"]:
            raise ValueError(f"Unexpected Tosa boundary columns: {reader.fieldnames}")
        for row in reader:
            chrom, start, end = parse_coordinate(row["Boundary"])
            if row["Type"] == "5p":
                site = start + boundary_anchor
            elif row["Type"] == "3p":
                site = end - boundary_anchor
            else:
                raise ValueError(f"Unexpected Tosa boundary type: {row['Type']}")
            key = f"{chrom}:{site}-{site + 1}"
            if key in wanted_boundaries:
                counts[key] += int(row["Count"])
    return counts


def merge_tosa_samples(prefixes, ri_event, output, boundary_anchor=1):
    wanted_boundaries = ri_boundary_ids(ri_event)
    samples = sorted(prefixes)
    counts = {sample: read_tosa_counts(prefixes[sample], boundary_anchor, wanted_boundaries)
              for sample in samples}
    ids = set(wanted_boundaries)
    for sample_counts in counts.values():
        ids.update(sample_counts)
    rows = []
    for key in ids:
        chrom, start, end = parse_coordinate(key)
        rows.append((chrom, start, end, key))
    rows.sort(key=lambda row: (row[0], row[1], row[2], row[3]))
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    with open(output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chr", "start", "end", "ID", *samples])
        for chrom, start, end, key in rows:
            writer.writerow([chrom, start, end, key,
                             *(counts[sample].get(key, 0) for sample in samples)])
    logger.info("Wrote %d junction and RI boundary rows to %s", len(rows), output)


def add_tosa_options(parser, require_gtf):
    parser.add_argument("-g", "--gtf", required=require_gtf, help="GTF used to generate events")
    parser.add_argument("-p", "--processors", type=int, default=1)
    parser.add_argument("-a", "--anchor", type=int, default=8)
    parser.add_argument("-b", "--boundary-anchor", type=int, default=1)
    parser.add_argument("-m", "--min-intron", type=int, default=70)
    parser.add_argument("-M", "--max-intron", type=int, default=500000)
    parser.add_argument("-s", "--strand", default="XS",
                        help="XS, RF, FR, unstranded, or legacy 0/1/2")
    parser.add_argument("-v", "--verbose", action="store_true")


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="subcommand")
    parser.add_argument("-i", "--input", help="Experiment table")
    parser.add_argument("-r", "--ri-event", help="EVENT_RI.txt")
    parser.add_argument("-o", "--output", help="Shiba junctions.bed")
    add_tosa_options(parser, False)

    run = subparsers.add_parser("run", help="Run Tosa for one BAM")
    run.add_argument("--bam", required=True)
    run.add_argument("--prefix", required=True)
    add_tosa_options(run, True)

    merge = subparsers.add_parser("merge", help="Merge per-sample Tosa outputs")
    merge.add_argument("-i", "--input", required=True, help="Experiment table")
    merge.add_argument("-d", "--directory", required=True, help="Directory of sample prefixes")
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
        run_tosa(args.bam, args.gtf, args.prefix, args.processors, args.anchor,
                 args.boundary_anchor, args.min_intron, args.max_intron, args.strand)
    elif args.subcommand == "merge":
        samples = read_experiment_table(args.input)
        prefixes = {sample: os.path.join(args.directory, sample) for sample in samples}
        merge_tosa_samples(prefixes, args.ri_event, args.output, args.boundary_anchor)
    else:
        if not all((args.input, args.ri_event, args.output, args.gtf)):
            parser.error("-i, -r, -o, and -g are required for the all-in-one run")
        samples = read_experiment_table(args.input)
        output_dir = os.path.dirname(os.path.abspath(args.output))
        os.makedirs(output_dir, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="tosa_", dir=output_dir) as tmp:
            prefixes = {}
            for sample, bam in samples.items():
                prefixes[sample] = os.path.join(tmp, sample)
                run_tosa(bam, args.gtf, prefixes[sample], args.processors, args.anchor,
                         args.boundary_anchor, args.min_intron, args.max_intron, args.strand,
                         os.path.join(output_dir, "logs", f"{sample}_tosa.log"))
            merge_tosa_samples(prefixes, args.ri_event, args.output, args.boundary_anchor)


if __name__ == "__main__":
    main()
