import os


localrules:
    alignment_qa_report,


# QC and metadata with NanoPlot
rule sample_qa_plot:
    input:
        fastq=lambda wildcards: get_mapped_reads_input(
            samples["sample"][wildcards.sample]
        ),
    output:
        scatter=report(
            "NanoPlot/{sample}/NanoPlot-report.html",
            category="Quality control",
            subcategory="NanoPlot",
            caption="../report/nanoplot_sample_report.rst",
            labels={
                "model": "NanoPlot",
                "figure": "{sample}",
            },
        ),
    params:
        outdir=lambda wildcards: f"NanoPlot/{wildcards.sample}",
    log:
        "logs/NanoPlot/{sample}.log",
    conda:
        "../envs/nanoplot.yml"
    threads: 4
    resources:
        mem_mb=lambda wildcards, input: max(
            1800, int(((os.path.getsize(input[0]) >> 20) * 0.2))
        ),
    shell:
        "NanoPlot --threads {threads} --tsv_stats --format svg "
        "--fastq {input.fastq} --outdir {params.outdir} 2> {log}"


rule total_sample_qa_plot:
    input:
        aggregate_input(samples["sample"]),
    output:
        scatter=report(
            "NanoPlot/NanoPlot-report.html",
            category="Quality control",
            subcategory="NanoPlot",
            caption="../report/nanoplot_all_samples_report.rst",
            labels={
                "model": "NanoPlot",
                "figure": "All samples",
            },
        ),
    # This parameter is in line with the Snakemake docs 8.20.3 guideline on how to avoid having parameters as output prefixes
    params:
        outdir=lambda wildcards, output: output[0][:-21],
    log:
        "logs/NanoPlot/all_samples.log",
    conda:
        "../envs/nanoplot.yml"
    shell:
        "NanoPlot --threads {threads} --tsv_stats --format svg "
        "--fastq {input} --outdir {params.outdir} 2> {log}"


rule alignment_qa:
    input:
        bam="sorted_alignments/{sample}_sorted.bam",
    output:
        directory("QC/qualimap/{sample}"),
    log:
        "logs/qualimap/{sample}.log",
    wrapper:
        "v4.4.0/bio/qualimap/bamqc"


# this is a dummy rule to create input for the report because the QualiMap wrapper only accepts directories as valid output
rule alignment_qa_report:
    input:
        map_qc=rules.alignment_qa.output,
    output:
        qm_report=report(
            "qualimap/{sample}/qualimapReport.html",
            category="Quality control",
            subcategory="QualiMap",
            caption="../report/qualimap.rst",
            labels={
                "model": "QualiMap",
                "figure": "{sample}",
            },
        ),
    log:
        "logs/qualimap/{sample}_report.log",
    conda:
        "../envs/base.yml"
    shell:
        "cp -a QC/qualimap/{wildcards.sample} qualimap/ 2> {log}"


rule bam_stats:
    input:
        bam="alignments/{sample}.bam",
    output:
        "QC/bamstats/{sample}.txt",
    log:
        "logs/samtools/bamstats_{sample}.log",
    params:
        extra=config["samtools"]["bamstats_opts"],
    resources:
        mem_mb_per_cpu=lambda wildcards, input, threads: max(
            1800, int(((os.path.getsize(input[0]) >> 20) * 0.2) / threads)
        ),
    wrapper:
        "v3.13.4/bio/samtools/stats"
