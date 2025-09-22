localrules:
    download_genome,
    download_annotation,
    get_annotation,
    get_genome,


rule download_genome:
    output:
        temp("references/ncbi_dataset_genome.zip"),
    params:
        accession=config["ref"]["accession"],
    log:
        "logs/refs/download_genome.log",
    conda:
        "../envs/reference.yml"
    shell:
        """
        datasets download genome accession {params.accession} --include genome &> {log} && mv ncbi_dataset.zip {output}
        """


rule download_annotation:
    output:
        temp("references/ncbi_dataset_annotation.zip"),
    params:
        accession=config["ref"]["accession"],
    log:
        "logs/refs/download_annotation.log",
    conda:
        "../envs/reference.yml"
    shell:
        """
        datasets download genome accession {params.accession} --include gff3 &> {log} && mv ncbi_dataset.zip {output}
        """


rule get_genome:
    input:
        lambda wildcards: get_reference_files(config).get("genome"),
    output:
        temp("references/genomic.fa"),
    params:
        accession=config["ref"]["accession"],
    log:
        "logs/refs/get_genome.log",
    conda:
        "../envs/reference.yml"
    script:
        "../scripts/extract_refs.py"


rule get_annotation:
    input:
        lambda wildcards: get_reference_files(config).get("annotation"),
    output:
        temp("references/genomic.gff"),
    params:
        accession=config["ref"]["accession"],
    log:
        "logs/refs/get_annotation.log",
    conda:
        "../envs/reference.yml"
    script:
        "../scripts/extract_refs.py"