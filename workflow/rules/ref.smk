localrules:
    download_ncbi_genome,
    download_ncbi_annotation,
    download_ensembl_genome,
    download_ensembl_annotation,
    get_annotation,
    get_genome,

# the download rules may fail frequently due to network issues. The ruleorder
# directive ensures that Snakemake try those rules in sucession until one works.
ruleorder: download_ncbi_genome > download_ensembl_genome > get_genome
ruleorder: download_ncbi_annotation > download_ensembl_annotation > get_annotation

rule download_ncbi_genome:
    output:
        "references/ncbi_dataset_genome.zip"
    retries: 3
    cache: "omit-software"
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
        "references/ncbi_dataset_annotation.zip"
    retries: 3
    cache: "omit-software"
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
        "references/ensembl_genome.fa",
    retries: 3
    cache: "omit-software"
    params:
        species=config["ref"]["ensembl_species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    log:
        "logs/refs/download_ensembl_genome.log",
    wrapper:
        "v7.5.0/bio/reference/ensembl-sequence"


rule download_ensembl_annotation:
    output:
        "references/ensembl_annotation.gff3"
    retries: 3
    params:
        species=config["ref"]["ensembl_species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    log:
        "logs/refs/download_ensembl_annotation.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v7.5.0/bio/reference/ensembl-annotation"


rule get_genome:
    input:
        lambda wildcards: get_reference_files(config).get("genome"),
    retries: 3
    output:
        "references/genomic.fa"
    cache: "omit-software"
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
    retries: 3
    output:
        "references/genomic.gff"
    cache: "omit-software"
    params:
        accession=config["ref"]["accession"],
    log:
        "logs/refs/get_annotation.log",
    conda:
        "../envs/reference.yml"
    script:
        "../scripts/extract_refs.py"
