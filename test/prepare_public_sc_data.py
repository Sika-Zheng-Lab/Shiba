"""Prepare a regional 10x PBMC execution fixture (not a biological comparison).

Source: https://www.10xgenomics.com/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0
10x Genomics, 2018, CC BY 4.0. The two barcode groups are artificial splits.
"""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile
import urllib.request

import pysam
import yaml

BASE = "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--gtf", type=Path, required=True, help="GRCh38 GTF with chromosome 1")
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    bam = root / "pbmc.bam"
    if not bam.exists():
        subprocess.run(["samtools", "view", "-b", "-o", str(bam), BASE + "_possorted_genome_bam.bam", "1:1-10000000"], cwd=root, check=True)
    pysam.index(str(bam))
    request = urllib.request.Request(BASE + "_filtered_feature_bc_matrix.tar.gz", headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(request, timeout=120) as response:
        archive = response.read()
    with tarfile.open(fileobj=io.BytesIO(archive), mode="r:gz") as handle:
        member = next(m for m in handle.getmembers() if m.name.endswith("barcodes.tsv.gz"))
        barcodes = set(gzip.decompress(handle.extractfile(member).read()).decode().splitlines())
    with pysam.AlignmentFile(bam) as handle:
        if handle.get_reference_length("1") != 248956422:
            raise ValueError("Expected GRCh38 chromosome 1")
        tagged = 0
        observed = set()
        for read in handle:
            if read.has_tag("CB") and read.has_tag("UB") and read.get_tag("CB") in barcodes:
                tagged += 1
                observed.add(read.get_tag("CB"))
    with (root / "barcodes.tsv").open("w") as handle:
        handle.write("barcode\tgroup\n")
        for i, barcode in enumerate(sorted(observed)):
            handle.write(f"{barcode}\t{'A' if i % 2 else 'B'}\n")
    with args.gtf.open() as source, (root / "annotation.gtf").open("w") as target:
        for line in source:
            fields = line.rstrip().split("\t")
            if len(fields) == 9 and fields[0] == "1" and int(fields[3]) >= 1 and int(fields[4]) <= 10000000:
                target.write(line)
    (root / "experiment.tsv").write_text(f"sample\talignment\tbarcode\npbmc\t{bam}\t{root / 'barcodes.tsv'}\n")
    config = dict(workdir=str(root / "scshiba"), container="unused-local-runtime", gtf=str(root / "annotation.gtf"),
                  experiment_table=str(root / "experiment.tsv"), strand="unstranded", minimum_reads=1,
                  reference_group="A", alternative_group="B", excel=True)
    (root / "config.yaml").write_text(yaml.safe_dump(config))
    manifest = dict(source=BASE + "_possorted_genome_bam.bam", region="1:1-10000000", tagged_reads=tagged,
                    observed_cell_barcodes=len(observed), bam_sha256=hashlib.sha256(bam.read_bytes()).hexdigest(),
                    grouping="Alternating sorted cell barcodes; technical validation only", license="CC BY 4.0", attribution="10x Genomics, 1k PBMCs v3, 2018")
    (root / "input-manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps(manifest))


if __name__ == "__main__":
    main()
