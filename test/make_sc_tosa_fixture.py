"""Generate two tiny tagged single-cell BAMs for the Tosa workflow test.

Requires pysam. The same CB appears in both libraries with different groups.
"""

import argparse
import csv
import json
from pathlib import Path

import pysam


def make_read(name, start, cigar, barcode, umi):
    read = pysam.AlignedSegment()
    read.query_name = name
    read.query_sequence = "A" * 100
    read.flag = 0
    read.reference_id = 0
    read.reference_start = start
    read.mapping_quality = 255
    read.cigarstring = cigar
    read.query_qualities = pysam.qualitystring_to_array("I" * 100)
    read.set_tag("CB", barcode)
    read.set_tag("UB", umi)
    read.set_tag("NH", 1)
    return read


def write_bam(path, specs):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr1", "LN": 10000}]}
    reads = [make_read(*spec) for spec in specs]
    with pysam.AlignmentFile(path, "wb", header=header) as handle:
        for read in sorted(reads, key=lambda item: item.reference_start):
            handle.write(read)
    pysam.index(str(path))


def write_table(path, header, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def build_fixture(output_dir):
    directory = Path(output_dir).resolve()
    directory.mkdir(parents=True, exist_ok=True)

    # chr1:1200-1500 is the Shiba ID for the first intron. TX2 retains it.
    gtf = directory / "annotation.gtf"
    with open(gtf, "w") as handle:
        for tx, exons in (("TX1", [(1000, 1200), (1500, 1700), (2000, 2300)]),
                          ("TX2", [(1000, 1700), (2000, 2300)])):
            for start, end in exons:
                attributes = f'gene_id "GENE1"; transcript_id "{tx}"; gene_name "ExampleGene";'
                handle.write(f"chr1\tsynthetic\texon\t{start}\t{end}\t.\t+\t.\t{attributes}\n")

    first = [
        ("j1_ref_1", 1150, "50M299N50M", "AAAA-1", "U1"),
        ("j1_ref_2", 1150, "50M299N50M", "AAAA-1", "U2"),
        ("j1_ref_dup", 1150, "50M299N50M", "AAAA-1", "U2"),
        ("j1_alt_1", 1150, "50M299N50M", "BBBB-1", "U1"),
        ("j1_alt_2", 1150, "50M299N50M", "BBBB-1", "U2"),
        ("b5_ref", 1150, "100M", "AAAA-1", "B1"),
        ("b5_alt", 1150, "100M", "BBBB-1", "B1"),
        ("b5_alt_dup", 1150, "100M", "BBBB-1", "B1"),
        ("b3_ref", 1450, "100M", "AAAA-1", "B2"),
        ("b3_alt", 1450, "100M", "BBBB-1", "B2"),
        ("excluded", 1150, "50M299N50M", "CCCC-1", "U1"),
    ]
    second = [
        ("j1_alt_library2", 1150, "50M299N50M", "AAAA-1", "U1"),
        ("j2_alt_1", 1650, "50M299N50M", "AAAA-1", "V1"),
        ("j2_alt_2", 1650, "50M299N50M", "AAAA-1", "V2"),
        ("b5_alt_library2", 1150, "100M", "AAAA-1", "B1"),
        ("b3_alt_library2", 1450, "100M", "AAAA-1", "B2"),
    ]
    write_bam(directory / "library1.bam", first)
    write_bam(directory / "library2.bam", second)
    write_table(directory / "library1_barcodes.tsv", ["barcode", "group"],
                [["AAAA-1", "Ref"], ["BBBB-1", "Alt"]])
    write_table(directory / "library2_barcodes.tsv", ["barcode", "group"],
                [["AAAA-1", "Alt"]])
    write_table(directory / "experiment.tsv", ["sample", "alignment", "barcode"],
                [["library1", directory / "library1.bam", directory / "library1_barcodes.tsv"],
                 ["library2", directory / "library2.bam", directory / "library2_barcodes.tsv"]])
    expected = {
        "chr1:1200-1500": {"Ref": 2, "Alt": 3},
        "chr1:1700-2000": {"Ref": 0, "Alt": 2},
        "chr1:1200-1201": {"Ref": 1, "Alt": 2},
        "chr1:1499-1500": {"Ref": 1, "Alt": 2},
    }
    with open(directory / "expected_counts.json", "w") as handle:
        json.dump(expected, handle, indent=2)
        handle.write("\n")
    with open(directory / "config.yaml", "w") as handle:
        handle.write(
            f"workdir: {directory / 'scshiba'}\n"
            "container: docker://naotokubota/shiba:v1.0.0\n"
            f"gtf: {gtf}\n"
            f"experiment_table: {directory / 'experiment.tsv'}\n"
            "minimum_anchor_length: 6\n"
            "boundary_anchor_length: 1\n"
            "minimum_intron_length: 70\n"
            "maximum_intron_length: 500000\n"
            "strand: unstranded\n"
            "only_psi: false\n"
            "fdr: 0.05\n"
            "delta_psi: 0.1\n"
            "reference_group: Ref\n"
            "alternative_group: Alt\n"
            "minimum_reads: 1\n"
            "beta_binomial: false\n"
            "excel: false\n"
        )
    return directory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir")
    args = parser.parse_args()
    build_fixture(args.output_dir)


if __name__ == "__main__":
    main()
