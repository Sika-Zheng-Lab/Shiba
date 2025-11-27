
VERSION = "v0.8.0"

'''
SnakeScShiba: A snakemake-based workflow of scShiba

Usage:
    snakemake -s snakescshiba.smk --configfile config.yaml --cores <int> --use-singularity --singularity-args "--bind $HOME:$HOME"
'''

import os
from pathlib import Path
import sys
args = sys.argv

workdir: config["workdir"]
container: config["container"]
base_dir = os.path.dirname(workflow.snakefile)
command = " ".join(args)
# Replace snakefile path with the absolute path
command = command.replace(workflow.snakefile, os.path.join(base_dir, workflow.snakefile))
# Replace config yaml path with the absolute path
configfile_path= args[args.index("--configfile") + 1]
# Absolute path of the directory containing the config yaml file
configfile_dir_path = Path(configfile_path).resolve().parent
command = command.replace(configfile_path, os.path.join(str(configfile_dir_path), configfile_path))

rule all:
    input:
        event_all = expand("events/EVENT_{sample}.txt", sample = ["SE", "FIVE", "THREE", "MXE", "MSE", "AFE", "ALE"]),
        PSI = expand("results/PSI_{sample}.txt", sample = ["SE", "FIVE", "THREE", "MXE", "MSE", "AFE", "ALE"]),
        report = "report.json"

rule generate_report:
    input:
        event_all = expand("events/EVENT_{sample}.txt", sample = ["SE", "FIVE", "THREE", "MXE", "MSE", "AFE", "ALE"]),
        PSI = expand("results/PSI_{sample}.txt", sample = ["SE", "FIVE", "THREE", "MXE", "MSE", "AFE", "ALE"])
    output:
        report = "report.json"
    params:
        workdir = config["workdir"],
        version = VERSION,
        command = command,
        experiment_table = config["experiment_table"]
    shell:
        """
        export PYTHONPATH={base_dir}/src/lib:$PYTHONPATH
        python -c 'from general import generate_report; generate_report("SnakeScShiba", "{params.workdir}", "{params.version}", "{params.command}", "{params.experiment_table}")'
        """

rule gtf2event:
    input:
        gtf = config["gtf"]
    output:
        events = directory("events"),
        events_all = expand("events/EVENT_{sample}.txt", sample = ["SE", "FIVE", "THREE", "MXE", "MSE", "AFE", "ALE"])
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

rule sc2junc:
    input:
        config["experiment_table"]
    output:
        "junctions/junctions.bed"
    benchmark:
        "benchmark/sc2junc.txt"
    log:
        "log/sc2junc.log"
    params:
        base_dir = base_dir
    shell:
        """
        python {params.base_dir}/src/sc2junc.py \
        -i {input} \
        -o {output} \
        -v \
        >& {log}
        """

rule scpsi:
    input:
        junc = "junctions/junctions.bed",
        events_all = expand("events/EVENT_{sample}.txt", sample = ["SE", "FIVE", "THREE", "MXE", "MSE", "AFE", "ALE"])
    output:
        results = directory("results"),
        PSI = expand("results/PSI_{sample}.txt", sample = ["SE", "FIVE", "THREE", "MXE", "MSE", "AFE", "ALE"])
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
        python {params.base_dir}/src/scpsi_snakemake.py \
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
