
import os
import datetime
_version_path = os.path.join(os.path.dirname(workflow.snakefile), "VERSION")
with open(_version_path, "r") as _vf:
    VERSION = _vf.read().strip()

# Capture pipeline start time as ISO string for report generation
start_time_str = datetime.datetime.now().isoformat(timespec="microseconds")

'''
SnakeScShiba: A snakemake-based workflow of scShiba

Usage:
    snakemake -s snakescshiba.smk --configfile config.yaml --cores <int> --use-singularity --singularity-args "--bind $HOME:$HOME"
'''

import os
from pathlib import Path
import sys
import csv
import shlex
args = sys.argv

with open(config["experiment_table"], newline="") as _handle:
    _samples = {row["sample"]: row for row in csv.DictReader(_handle, delimiter="\t")}
EVENT_TYPES = ["SE", "FIVE", "THREE", "MXE", "RI", "MSE", "AFE", "ALE"]

workdir: config["workdir"]
container: config["container"]
base_dir = os.path.dirname(workflow.snakefile)

# Validate configuration before pipeline execution
sys.path.insert(0, os.path.join(base_dir, "src", "lib"))
from general import validate_config, apply_config_defaults
apply_config_defaults(config, mode="sc")
_validation_errors = validate_config(config, mode="sc")
if _validation_errors:
    for _err in _validation_errors:
        print(f"ERROR: Configuration error: {_err}", file=sys.stderr)
    print(f"ERROR: {len(_validation_errors)} configuration error(s) found. Aborting.", file=sys.stderr)
    sys.exit(1)

command = shlex.join([sys.executable, *args])

rule all:
    input:
        event_all = expand("events/EVENT_{sample}.txt", sample = EVENT_TYPES),
        PSI = expand("results/PSI_{sample}.txt", sample = EVENT_TYPES),
        report = "report.json",
        plots = "plots/summary.html"

rule generate_report:
    input:
        event_all = expand("events/EVENT_{sample}.txt", sample = EVENT_TYPES),
        PSI = expand("results/PSI_{sample}.txt", sample = EVENT_TYPES)
    output:
        report = "report.json"
    params:
        workdir = config["workdir"],
        version = VERSION,
        command = command,
        experiment_table = config["experiment_table"],
        start_time = start_time_str
    shell:
        """
        export PYTHONPATH={base_dir:q}/src/lib:$PYTHONPATH
        python -c 'import sys; from general import generate_report; generate_report("SnakeScShiba", *sys.argv[1:])' {params.workdir:q} {params.version:q} {params.command:q} {params.experiment_table:q} {params.start_time:q}
        """

rule gtf2event:
    input:
        gtf = config["gtf"]
    output:
        events = directory("events"),
        events_all = expand("events/EVENT_{sample}.txt", sample = EVENT_TYPES)
    threads:
        workflow.cores
    benchmark:
        "benchmark/gtf2event.txt"
    log:
        "log/gtf2event.log"
    params:
        base_dir = base_dir
    shell:
        """
        python {params.base_dir}/src/gtf2event.py \
        -i {input.gtf} \
        -o {output.events} \
        -p {threads} \
        -v \
        >& {log}
        """

rule tosa_single:
    wildcard_constraints:
        sample = "|".join(_samples)
    input:
        alignment = lambda wildcards: _samples[wildcards.sample]["alignment"],
        barcode = lambda wildcards: _samples[wildcards.sample]["barcode"],
        gtf = config["gtf"]
    output:
        matrix = temp("junctions/{sample}_matrix.mtx.gz"),
        barcodes = temp("junctions/{sample}_barcodes.tsv.gz"),
        features = temp("junctions/{sample}_features.tsv.gz"),
        junction_detail = temp("junctions/{sample}_junction_barcodes.tsv.gz"),
        boundary_matrix = temp("junctions/{sample}_boundary_matrix.mtx.gz"),
        boundary_barcodes = temp("junctions/{sample}_boundary_barcodes.tsv.gz"),
        boundary_features = temp("junctions/{sample}_boundary_features.tsv.gz"),
        boundary_detail = temp("junctions/{sample}_boundary_barcodes_detail.tsv.gz")
    threads:
        8
    benchmark:
        "benchmark/sc2junc/{sample}_tosa.txt"
    log:
        "log/sc2junc/{sample}_tosa.log"
    params:
        base_dir = base_dir,
        anchor = config.get("minimum_anchor_length", 8),
        boundary_anchor = config.get("boundary_anchor_length", 1),
        min_intron = config.get("minimum_intron_length", 20),
        max_intron = config.get("maximum_intron_length", 500000),
        strand = config.get("strand", "unstranded")
    shell:
        """
        python {params.base_dir}/src/sc2junc.py run \
        --alignment {input.alignment} \
        --barcode {input.barcode} \
        --prefix junctions/{wildcards.sample} \
        -g {input.gtf} -p {threads} \
        -a {params.anchor} -b {params.boundary_anchor} \
        -m {params.min_intron} -M {params.max_intron} \
        -s {params.strand} -v >& {log}
        """

rule sc2junc:
    input:
        experiment = config["experiment_table"],
        junctions = expand("junctions/{sample}_junction_barcodes.tsv.gz", sample = _samples),
        boundaries = expand("junctions/{sample}_boundary_barcodes_detail.tsv.gz", sample = _samples),
        RI = "events/EVENT_RI.txt"
    output:
        "junctions/junctions.bed"
    benchmark:
        "benchmark/sc2junc.txt"
    log:
        "log/sc2junc.log"
    params:
        base_dir = base_dir,
        boundary_anchor = config.get("boundary_anchor_length", 1)
    shell:
        """
        python {params.base_dir}/src/sc2junc.py merge \
        -i {input.experiment} -d junctions -r {input.RI} \
        -o {output} \
        -b {params.boundary_anchor} -v >& {log}
        """

rule scpsi:
    input:
        junc = "junctions/junctions.bed",
        events_all = expand("events/EVENT_{sample}.txt", sample = EVENT_TYPES)
    output:
        results = directory("results"),
        summary = [] if config["only_psi"] else ["results/summary.txt"],
        PSI = expand("results/PSI_{sample}.txt", sample = EVENT_TYPES)
    threads:
        1
    benchmark:
        "benchmark/scpsi.txt"
    log:
        "log/scpsi.log"
    params:
        base_dir = base_dir
    shell:
        """
        python {params.base_dir}/src/scpsi.py \
        -p {threads} \
        -f {config[fdr]} \
        -d {config[delta_psi]} \
        -m {config[minimum_reads]} \
        -r {config[reference_group]} \
        -a {config[alternative_group]} \
        --onlypsi {config[only_psi]} \
        --excel {config[excel]} \
        -v \
        {input.junc} \
        events \
        {output.results} >& {log}
        """

rule plots:
    input:
        PSI = expand("results/PSI_{sample}.txt", sample=EVENT_TYPES),
        summary = [] if config["only_psi"] else ["results/summary.txt"]
    output:
        "plots/summary.html"
    params:
        script = os.path.join(base_dir, "src", "plots.py"),
        command = command,
        only_psi = "--only-psi" if config["only_psi"] else ""
    log:
        "log/plots.log"
    shell:
        """
        python {params.script:q} --single-cell {params.only_psi} \
        -i results -o plots -s {params.command:q} > {log:q} 2>&1
        """
