import os


localrules:
    merge_read_counts,
    transcriptid_to_gene,


rule count_reads:
    input:
        bam="alignments/{sample}.bam",
        trs="transcriptome/corrected_transcriptome.fa",
        annotation="references/standardized_genomic.gff",
    output:
        quant="counts/{sample}/{sample}.quant",
    log:
        "logs/count_reads/{sample}.log",
    conda:
        "../envs/oarfish.yml"
    threads: 8
    resources:
        mem_mb_per_cpu=lambda wildcards, input, threads: max(
            1800, int(((os.path.getsize(input[0]) >> 20) * 2) / threads)
        ),
    params:
        outdir=lambda wildcards, output: os.path.splitext(output.quant)[0],
        seqtech=config["quant"]["oarfish_seqtech"],
    shell:
        """
        oarfish -j {threads} -a {input.bam} -o {params.outdir} --seq-tech {params.seqtech} --quiet &>{log}
        """
        #"""
        #salmon --no-version-check quant -p {threads} {params.longreads} \
        #    -t {input.trs} -l {params.libtype} -a {input.bam} -o {params.outdir} 2>{log}
        #"""


rule merge_read_counts:
    input:
        count_tsvs=expand("counts/{sample}/{sample}.quant", sample=samples["sample"]),
    output:
        "merged/all_counts.tsv",
    log:
        "logs/merge_count.log",
    conda:
        "../envs/pandas.yml"
    script:
        "../scripts/merge_count_tsvs.py"


rule transcriptid_to_gene:
    input:
        all_counts="merged/all_counts.tsv",
        annotation="references/standardized_genomic.gff",
    output:
        all_counts="merged/all_counts_gene.tsv",
        plot=report(
            "merged/transcriptid_to_gene_plot.svg",
            category="Quality control",
            subcategory="Transcript Naming",
            caption="../report/name_stats.rst",
            labels={
                "model": "Matplotlib",
                "figure": "Naming Stats",
            },
        ),
    log:
        "logs/transcriptid_to_gene.log",
    conda:
        "../envs/pydeseq2.yml"
    script:
        "../scripts/transcriptid_to_gene.py"
