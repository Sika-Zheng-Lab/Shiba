
VERSION = "v0.5.2"

'''
SnakeScShiba: A snakemake-based workflow of scShiba

Usage:
    snakemake -s snakescshiba.smk --configfile config.yaml --cores <int> --use-singularity --singularity-args "--bind $HOME:$HOME"
'''

workdir: config["workdir"]
container: config["container"]
base_dir = os.path.dirname(workflow.snakefile)

rule all:
    input:
        event_all = expand("events/EVENT_{sample}.txt", sample = ["SE", "FIVE", "THREE", "MXE", "MSE", "AFE", "ALE"]),
        PSI = expand("results/PSI_{sample}.txt", sample = ["SE", "FIVE", "THREE", "MXE", "MSE", "AFE", "ALE"])
    params:
        version = VERSION
    shell:
        """
        echo \
        '

            All analysis finished successfully!
            SnakeScShiba version: {params.version}
            https://github.com/Sika-Zheng-Lab/Shiba

        '
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
