localrules:
    download_ncbi_genome,
    download_ncbi_annotation,
    download_ensembl_genome,
    download_ensembl_annotation,
    get_annotation,
    get_genome,


rule download_ncbi_genome:
    output:
        "references/ncbi_dataset_genome.zip",
    log:
        "logs/refs/download_ncbi_genome.log",
    # retries: 3
    conda:
        "../envs/reference.yml"
    params:
        accession=config["ref"]["accession"],
    shell:
        """
        datasets download genome accession {params.accession} --include genome --filename {output} --no-progressbar &>{log}
        """


rule download_ncbi_annotation:
    output:
        "references/ncbi_dataset_annotation.zip",
    log:
        "logs/refs/download_ncbi_annotation.log",
    retries: 3
    conda:
        "../envs/reference.yml"
    params:
        accession=config["ref"]["accession"],
    shell:
        """
        datasets download genome accession {params.accession} --include gff3 --filename {output} --no-progressbar &>{log}
        """


rule download_ensembl_genome:
    output:
        "references/ensembl_genome.fa",
    log:
        "logs/refs/download_ensembl_genome.log",
    retries: 3
    params:
        species=config["ref"]["ensembl_species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    wrapper:
        "v7.5.0/bio/reference/ensembl-sequence"


rule download_ensembl_annotation:
    output:
        "references/ensembl_annotation.gff3",
    log:
        "logs/refs/download_ensembl_annotation.log",
    retries: 3
    params:
        species=config["ref"]["ensembl_species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    wrapper:
        "v7.5.0/bio/reference/ensembl-annotation"


rule get_genome:
    input:
        lambda wildcards: get_reference_files(config).get("genome"),
    output:
        "references/genomic.fa",
    log:
        "logs/refs/get_genome.log",
    retries: 3
    conda:
        "../envs/reference.yml"
    params:
        accession=config["ref"]["accession"],
    script:
        "../scripts/extract_refs.py"


rule get_annotation:
    input:
        lambda wildcards: get_reference_files(config).get("annotation"),
    output:
        "references/genomic.gff",
    log:
        "logs/refs/get_annotation.log",
    retries: 3
    conda:
        "../envs/reference.yml"
    params:
        accession=config["ref"]["accession"],
    script:
        "../scripts/extract_refs.py"
