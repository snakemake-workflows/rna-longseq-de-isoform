import os


localrules:
    compress_nplot,
    compress_nplot_all,
    compress_map_qc,


configfile: "config/config.yml"


inputdir = config["inputdir"]  # "/lustre/project/m2_zdvhpc/transcriptome_data/"


# QC and metadata with NanoPlot

if config["summary"] == "None":
    sample_QC = (
        (expand("QC/NanoPlot/{sample}.tar.gz", sample=samples["sample"]),),
        "QC/NanoPlot/all_samples.tar.gz",
    )
else:
    sample_QC = "QC/NanoPlot/summary.tar.gz"


if config["summary"] == "None":

    rule plot_samples:
        input:
            fastq=lambda wildcards: get_mapped_reads_input(
                samples["sample"][wildcards.sample]
            ),
        output:
            directory("NanoPlot/{sample}"),
        log:
            "logs/NanoPlot/{sample}.log",
        resources:
            cpus_per_task=min(len({input}), config["max_cpus"]),  #problem with max(len(input.fastq),39)
        conda:
            "../envs/env.yml"
        shell:
            "mkdir {output}; "
            "NanoPlot -t {resources.cpus_per_task} --tsv_stats -f svg "
            "--fastq {input.fastq} -o {output} 2> {log}"

    rule plot_all_samples:
        input:
            aggregate_input(samples["sample"]),
        output:
            directory("NanoPlot/all_samples"),
        log:
            "logs/NanoPlot/all_samples.log",
        conda:
            "../envs/env.yml"
        shell:
            "mkdir {output}; "
            "NanoPlot -t {resources.cpus_per_task} --tsv_stats -f svg "
            "--fastq {input} -o {output} 2> {log}"

    rule compress_nplot:
        input:
            samples=rules.plot_samples.output,
        output:
            "QC/NanoPlot/{sample}.tar.gz",
        log:
            "logs/NanoPlot/compress_{sample}.log",
        conda:
            None
        shell:
            "tar zcvf {output} {input} &> {log}"

    rule compress_nplot_all:
        input:
            all_samples=rules.plot_all_samples.output,
        output:
            "QC/NanoPlot/all_samples.tar.gz",
        log:
            "logs/NanoPlot/compress_all_samples.log",
        conda:
            None
        shell:
            "tar zcvf {output} {input} &> {log}"

else:

    rule plot_samples:
        input:
            summary=config["summary"],
        output:
            directory("NanoPlot"),
        log:
            "logs/NanoPlot/NanoPlot.log",
        conda:
            "../envs/env.yml"
        shell:
            "mkdir {output}; "
            "NanoPlot -t {resources.cpus_per_task} --barcoded --tsv_stats "
            "--summary {input.summary} -o {output} 2> {log}"

    rule compress_nplot:
        input:
            samples=rules.plot_samples.output,
        output:
            "QC/NanoPlot/summary.tar.gz",
        log:
            "logs/NanoPlot/compress_summary.log",
        conda:
            None
        shell:
            "tar zcvf {output} {input} &> {log}"


rule sam_sort:
    input:
        sam="alignments/{sample}.sam",
    output:
        "sorted_alignments/{sample}_sorted.bam",
    log:
        "logs/samtools/samsort_{sample}.log",
    conda:
        "../envs/env.yml"
    shell:
        "samtools sort -@ {resources.cpus_per_task} {input.sam} -o {output} -O bam &> {log}"


rule map_qc:
    input:
        sorted_bam=rules.sam_sort.output,
    output:
        directory("QC/qualimap/{sample}"),
    log:
        "logs/qualimap/{sample}.log",
    conda:
        "../envs/env.yml"
    wrapper:
        "v3.12.1/bio/qualimap/bamqc"


rule compress_map_qc:
    input:
        map_qc=rules.map_qc.output,
    output:
        "QC/qualimap/{sample}.tar.gz",
    log:
        "logs/qualimap/compress_{sample}.log",
    conda:
        None
    shell:
        "tar zcvf {output} {input} &> {log}"


rule sam_stats:
    input:
        bam="sorted_alignments/{sample}.bam",
    output:
        "QC/samstats/{sample}.txt",
    log:
        "logs/samtools/samstats_{sample}.log",
    conda:
        "../envs/env.yml"
    shell:
        """
            samtools stats -@ {resources.cpus_per_task} {input.bam} > {output} 2> {log}
        """
