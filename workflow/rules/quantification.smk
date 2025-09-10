import os


localrules:
    merge_read_counts,


rule count_reads:
    input:
        bam="alignments/{sample}.bam",
        trs="transcriptome/corrected_transcriptome.fa",
    output:
        tsv="counts/{sample}_salmon/quant.sf",
    params:
        outdir=lambda wildcards: f"counts/{wildcards.sample}_salmon",
        libtype=config["quant"]["salmon_libtype"],
    log:
        "logs/salmon/{sample}.log",
    conda:
        "../envs/salmon.yml"
    threads: 8
    resources:
        mem_mb_per_cpu=lambda wildcards, input, threads: max(
            1800, int(((os.path.getsize(input[0]) >> 20) * 2) / threads)
        ),
    shell:
        """
        salmon --no-version-check quant --ont -p {threads} \
        -t {input.trs} -l {params.libtype} -a {input.bam} -o {params.outdir} 2> {log}
        """


rule merge_read_counts:
    input:
        count_tsvs=expand("counts/{sample}_salmon/quant.sf", sample=samples["sample"]),
    output:
        "merged/all_counts.tsv",
    log:
        "logs/merge_count.log",
    conda:
        "../envs/pandas.yml"
    script:
        "../scripts/merge_count_tsvs.py"
