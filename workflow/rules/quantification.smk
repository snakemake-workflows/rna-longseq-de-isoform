localrules:
    merge_counts,


rule count_reads:
    input:
        bam="sorted_alignments/{sample}.bam",
        trs="transcriptome/transcriptome.fa",
    output:
        tsv=directory("counts/{sample}_salmon"),
    params:
        libtype=config["salmon_libtype"],
    log:
        "logs/salmon/{sample}.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        salmon --no-version-check quant --ont -p {resources.cpus_per_task} \
        -t {input.trs} -l {params.libtype} -a {input.bam} -o {output.tsv} 2> {log}
    """


rule merge_counts:
    input:
        count_tsvs=expand("counts/{sample}_salmon/quant.sf", sample=samples["sample"]),
    output:
        "merged/all_counts.tsv",
    log:
        "logs/merge_count.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/merge_count_tsvs.py"
