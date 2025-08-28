
VERSION = "v0.7.1"

'''
SnakeShiba: A snakemake-based workflow of Shiba for differential RNA splicing analysis between two groups of samples

Usage:
    snakemake -s snakeshiba.smk --configfile config.yaml --cores 4 --use-singularity --singularity-args "--bind $HOME:$HOME"
'''

import os
from pathlib import Path
import sys
args = sys.argv

def load_experiment(file):
    experiment_dict = {}
    with open(file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("sample"):
                continue
            try:
                sample, bam, group, technology = line.split(maxsplit=3)
            except ValueError:
                sample, bam, group = line.split(maxsplit=2)
                technology = "short"
            experiment_dict[sample] = {"bam": bam, "group": group, "technology": technology}
    return experiment_dict
experiment_dict = load_experiment(config["experiment_table"])
juncfiles_list = []

base_dir = os.path.dirname(workflow.snakefile)
command = " ".join(args)
# Replace snakefile path with the absolute path
command = command.replace(workflow.snakefile, os.path.join(base_dir, workflow.snakefile))
# Replace config yaml path with the absolute path
configfile_path= args[args.index("--configfile") + 1]
# Absolute path of the directory containing the config yaml file
configfile_dir_path = Path(configfile_path).resolve().parent
command = command.replace(configfile_path, os.path.join(str(configfile_dir_path), configfile_path))

workdir: config["workdir"]
container: config["container"]

rule all:
    input:
        event_all = expand("events/EVENT_{sample}.txt", sample = ["SE", "FIVE", "THREE", "MXE", "RI", "MSE", "AFE", "ALE"]),
        PSI = expand("results/splicing/PSI_{sample}.txt", sample = ["SE", "FIVE", "THREE", "MXE", "RI", "MSE", "AFE", "ALE"]),
        summary = "plots/summary.html",
        barplot_splicing_summary = "plots/png/barplot_splicing_summary.png",
        tpm = "results/expression/TPM.txt",
        cpm = "results/expression/CPM.txt",
        counts = "results/expression/counts.txt",
        deg = "results/expression/DEG.txt",
        tpm_pca = "results/pca/tpm_pca.tsv",
        tpm_contribution = "results/pca/tpm_contribution.tsv",
        psi_pca = "results/pca/psi_pca.tsv",
        psi_contribution = "results/pca/psi_contribution.tsv",
        report = "report.json"

rule generate_report:
    input:
        event_all = expand("events/EVENT_{sample}.txt", sample = ["SE", "FIVE", "THREE", "MXE", "RI", "MSE", "AFE", "ALE"]),
        PSI = expand("results/splicing/PSI_{sample}.txt", sample = ["SE", "FIVE", "THREE", "MXE", "RI", "MSE", "AFE", "ALE"]),
        summary = "plots/summary.html",
        barplot_splicing_summary = "plots/png/barplot_splicing_summary.png",
        tpm = "results/expression/TPM.txt",
        cpm = "results/expression/CPM.txt",
        counts = "results/expression/counts.txt",
        deg = "results/expression/DEG.txt",
        tpm_pca = "results/pca/tpm_pca.tsv",
        tpm_contribution = "results/pca/tpm_contribution.tsv",
        psi_pca = "results/pca/psi_pca.tsv",
        psi_contribution = "results/pca/psi_contribution.tsv"
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
        python -c 'from general import generate_report; generate_report("SnakeShiba", "{params.workdir}", "{params.version}", "{params.command}", "{params.experiment_table}")'
        """

rule bam2gtf:
    wildcard_constraints:
        sample = "|".join(experiment_dict)
    input:
        bam = lambda wildcards: experiment_dict[wildcards.sample]["bam"],
        gtf = config["gtf"]
    output:
        temp("annotation/{sample}.gtf")
    threads:
        8
    benchmark:
        "benchmark/bam2gtf/{sample}.txt"
    log:
        "log/bam2gtf/{sample}.log"
    params:
        longread_option = lambda wildcards: "-L" if experiment_dict[wildcards.sample]["technology"] == "long" else ""
    shell:
        """
        stringtie -p {threads} \
        -G {input.gtf} \
        -o {output} \
        {params.longread_option} \
        {input.bam} >& {log}
        """

rule merge_gtf:
    wildcard_constraints:
        sample = "|".join(experiment_dict)
    input:
        gtflist = expand("annotation/{sample}.gtf", sample = experiment_dict),
        gtf = config["gtf"]
    output:
        "annotation/assembled.gtf"
    threads:
        workflow.cores
    benchmark:
        "benchmark/merge_gtf.txt"
    log:
        "log/merge_gtf.log"
    shell:
        """
        stringtie --merge \
        -p {threads} \
        -G {input.gtf} \
        -o {output} \
        {input.gtflist} >& {log}
        """

rule gtf2event:
    input:
        assembled_gtf = "annotation/assembled.gtf",
        gtf = config["gtf"]
    output:
        events = directory("events"),
        events_all = expand("events/EVENT_{sample}.txt", sample = ["SE", "FIVE", "THREE", "MXE", "RI", "MSE", "AFE", "ALE"])
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
        -i {input.assembled_gtf} \
        -r {input.gtf} \
        -o {output.events} \
        -p {threads} \
        -v \
        >& {log}
        """

rule bam2junc:
    wildcard_constraints:
        sample = "|".join(experiment_dict)
    input:
        bam = lambda wildcards: experiment_dict[wildcards.sample]["bam"]
    output:
        junc = temp("junctions/{sample}.junc")
    threads:
        1
    benchmark:
        "benchmark/bam2junc/{sample}_regtools.txt"
    log:
        "log/bam2junc/{sample}_regtools.log"
    shell:
        """
        regtools junctions extract \
        -a {config[minimum_anchor_length]} \
        -m {config[minimum_intron_length]} \
        -M {config[maximum_intron_length]} \
        -s {config[strand]} \
        -o {output.junc} \
        {input.bam} >& {log}
        """

rule make_RI_saf:
    input:
        "events"
    output:
        temp("junctions/RI.saf")
    shell:
        """
        cat {input}/EVENT_RI.txt | \
        cut -f 6,7 | \
        sed -e 1d | \
        awk -F'\t' -v OFS='\t' '{{split($1,l,":"); split(l[2],m,"-"); print l[1]":"m[1]"-"m[1]+1,l[1],m[1],m[1]+1,$2; print l[1]":"m[2]-1"-"m[2],l[1],m[2]-1,m[2],$2}}' | \
        awk '!a[$0]++' > {output}
        """

rule bam2junc_RI:
    wildcard_constraints:
        sample = "|".join(experiment_dict)
    input:
        RI = "junctions/RI.saf",
        bam = lambda wildcards: experiment_dict[wildcards.sample]["bam"]
    output:
        junc = temp("junctions/{sample}_exon-intron.junc"),
        junc_summary = temp("junctions/{sample}_exon-intron.junc.summary")
    threads:
        8
    benchmark:
        "benchmark/bam2junc/{sample}_featureCounts_RI.txt"
    log:
        "log/bam2junc/{sample}_featureCounts_RI.log"
    params:
        base_dir = base_dir,
        longread_option = lambda wildcards: "-l" if experiment_dict[wildcards.sample]["technology"] == "long" else ""
    shell:
        """
        python {params.base_dir}/src/bam2junc_RI_snakemake.py \
        -b {input.bam} \
        -r {input.RI} \
        -o {output.junc} \
        -t {threads} \
        {params.longread_option} \
        -v \
        &> {log}
        """

rule merge_junc:
    input:
        exonexon = expand("junctions/{sample}.junc", sample = experiment_dict),
        exonintron = expand("junctions/{sample}_exon-intron.junc", sample = experiment_dict)
    output:
        "junctions/junctions.bed"
    benchmark:
        "benchmark/merge_junc.txt"
    log:
        "log/merge_junc.log"
    params:
        base_dir = base_dir
    shell:
        """
        python {params.base_dir}/src/merge_junc_snakemake.py \
        --exonexon {input.exonexon} \
        --exonintron {input.exonintron} \
        --output {output} \
        -v \
        >& {log}
        """

rule psi:
    input:
        junc = "junctions/junctions.bed",
        events_all = expand("events/EVENT_{sample}.txt", sample = ["SE", "FIVE", "THREE", "MXE", "RI", "MSE", "AFE", "ALE"])
    output:
        results = directory("results/splicing"),
        PSI = expand("results/splicing/PSI_{sample}.txt", sample = ["SE", "FIVE", "THREE", "MXE", "RI", "MSE", "AFE", "ALE"]),
        PSI_matrix_sample = "results/splicing/PSI_matrix_sample.txt"
    threads:
        1
    benchmark:
        "benchmark/psi.txt"
    log:
        "log/psi.log"
    params:
        base_dir = base_dir
    shell:
        """
        python {params.base_dir}/src/psi_snakemake.py \
        -p {threads} \
        -g {config[experiment_table]} \
        -f {config[fdr]} \
        -d {config[delta_psi]} \
        -m {config[minimum_reads]} \
        -r {config[reference_group]} \
        -a {config[alternative_group]} \
        -i {config[individual_psi]} \
        -t {config[ttest]} \
        --onlypsi False \
        --onlypsi-group False \
        --excel {config[excel]} \
        -v \
        {input.junc} \
        events \
        {output.results} >& {log}
        """

rule expression_featureCounts:
    wildcard_constraints:
        sample = "|".join(experiment_dict)
    input:
        bam = lambda wildcards: experiment_dict[wildcards.sample]["bam"],
        gtf = config["gtf"]
    output:
        counts = temp("results/expression/{sample}_counts.txt"),
        counts_summary = temp("results/expression/{sample}_counts.txt.summary")
    threads:
        8
    benchmark:
        "benchmark/expression/{sample}_featureCounts.txt"
    log:
        "log/expression/{sample}_featureCounts.log"
    params:
        base_dir = base_dir,
        longread_option = lambda wildcards: "-l" if experiment_dict[wildcards.sample]["technology"] == "long" else ""
    shell:
        """
        python {params.base_dir}/src/expression_featureCounts_snakemake.py \
        -b {input.bam} \
        -g {input.gtf} \
        -o {output.counts} \
        -t {threads} \
        {params.longread_option} \
        -v \
        &> {log}
        """

rule expression_tpm:
    input:
        counts = expand("results/expression/{sample}_counts.txt", sample = experiment_dict)
    output:
        tpm = "results/expression/TPM.txt",
        cpm = "results/expression/CPM.txt",
        counts = "results/expression/counts.txt"
    benchmark:
        "benchmark/expression/TPM_CPM.txt"
    log:
        "log/expression/TPM_CPM.log"
    params:
        base_dir = base_dir
    shell:
        """
        python {params.base_dir}/src/tpm_snakemake.py \
        --countfiles {input.counts} \
        --output results/expression/ \
        --excel {config[excel]} \
        -v \
        &> {log}
        """

rule deseq2:
    input:
        counts = "results/expression/counts.txt"
    output:
        deseq2 = "results/expression/DEG.txt"
    benchmark:
        "benchmark/expression/DESeq2.txt"
    log:
        "log/expression/DESeq2.log"
    params:
        base_dir = base_dir
    shell:
        """
        python {params.base_dir}/src/deseq2_snakemake.py \
        --count {input.counts} \
        --experiment-table {config[experiment_table]} \
        --reference {config[reference_group]} \
        --alternative {config[alternative_group]} \
        --output {output.deseq2} \
        -v \
        &> {log}
        """

rule pca:
    input:
        tpm = "results/expression/TPM.txt",
        PSI_matrix_sample = "results/splicing/PSI_matrix_sample.txt"
    output:
        tpm_pca = "results/pca/tpm_pca.tsv",
        tpm_contribution = "results/pca/tpm_contribution.tsv",
        psi_pca = "results/pca/psi_pca.tsv",
        psi_contribution = "results/pca/psi_contribution.tsv"
    benchmark:
        "benchmark/pca/pca.txt"
    log:
        "log/pca.log"
    params:
        base_dir = base_dir
    shell:
        """
        python {params.base_dir}/src/pca.py \
        --input-tpm {input.tpm} \
        --input-psi {input.PSI_matrix_sample} \
        -g 3000 \
        -o results/pca \
        -v \
        &> {log}
        """

rule plots:
    input:
        PSI = expand("results/splicing/PSI_{sample}.txt", sample = ["SE", "FIVE", "THREE", "MXE", "RI", "MSE", "AFE", "ALE"]),
        tpm_pca = "results/pca/tpm_pca.tsv",
        tpm_contribution = "results/pca/tpm_contribution.tsv",
        psi_pca = "results/pca/psi_pca.tsv",
        psi_contribution = "results/pca/psi_contribution.tsv"
    output:
        summary = "plots/summary.html",
        barplot_splicing_summary = "plots/png/barplot_splicing_summary.png"
    benchmark:
        "benchmark/plots.txt"
    log:
        "log/plots.log"
    params:
        base_dir = base_dir,
        command = command
    shell:
        """
        python {params.base_dir}/src/plots.py \
        -i results \
        -e {config[experiment_table]} \
        -s "{params.command}" \
        -o plots \
        -v \
        >& {log}
        """
