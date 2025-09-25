localrules:
    download_ncbi_genome,
    download_ncbi_annotation,
    download_ensembl_genome,
    download_ensembl_annotation,
    get_annotation,
    get_genome,


rule download_ncbi_genome:
    output:
        temp("references/ncbi_dataset_genome.zip"),
    params:
        accession=config["ref"]["accession"],
    log:
        "logs/refs/download_ncbi_genome.log",
    conda:
        "../envs/reference.yml"
    shell:
        """
        datasets download genome accession {params.accession} --include genome &> {log} && mv ncbi_dataset.zip {output}
        """


rule download_ncbi_annotation:
    output:
        temp("references/ncbi_dataset_annotation.zip"),
    params:
        accession=config["ref"]["accession"],
    log:
        "logs/refs/download_ncbi_annotation.log",
    conda:
        "../envs/reference.yml"
    shell:
        """
        datasets download genome accession {params.accession} --include gff3 &> {log} && mv ncbi_dataset.zip {output}
        """


rule download_ensembl_genome:
    output:
        temp("references/ensembl_genome.fa"),
    params:
        species=config["ref"]["ensembl_species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    log:
        "logs/references/download_ensembl_genome.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v7.5.0/bio/reference/ensembl-sequence"


rule download_ensembl_annotation:
    output:
        temp("references/ensembl_annotation.gff3"),
    params:  
        species=config["ref"]["ensembl_species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    log:
        "logs/references/download_ensembl_annotation.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v7.5.0/bio/reference/ensembl-annotation"


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