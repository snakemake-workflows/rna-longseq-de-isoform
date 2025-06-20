localrules:
    generate_gene_query,
    get_indexed_db,
    get_protein_names,


rule get_indexed_db:
    output:
        "protein_annotation/index/UniRef.lba.gz",
    params:
        ref=f'{config["protein_annotation"]["uniref"]}',
    log:
        "logs/lambda/get_indexed_db.log",
    conda:
        "../envs/wget.yml"
    shell:
        """
        wget -nv -q -O UniRef.lba.gz {params.ref} && \
        mv UniRef.lba.gz {output} 2> {log}
        """


rule generate_gene_query:
    input:
        sorted_lfc_counts="de_analysis/sorted_normalized_counts.csv",
        transcriptome="transcriptome/transcriptome.fa",
    output:
        "protein_annotation/de_genes.fa",
    log:
        "logs/lambda/generate_gene_query.log",
    conda:
        "../envs/biopython.yml"
    script:
        "../scripts/get_de_genes.py"


rule blast_genes:
    input:
        indexed_db="protein_annotation/index/UniRef.lba.gz",
        query="protein_annotation/de_genes.fa",
    output:
        lambda_results=report(
            "protein_annotation/blast_results.m8",
            category="Protein Annotation Results",
            subcategory="Lambda Results",
            caption="../report/lambda_results.rst",
            labels={
                "list": "lambda_results",
            },
        ),
    params:
        num_matches=f'{config["protein_annotation"]["num_matches"]}',
    log:
        "logs/lambda/blast_genes.log",
    conda:
        "../envs/lambda3.yml"
    shell:
        "lambda3 searchp -q {input.query} -i {input.indexed_db} -o {output} -n {params.num_matches} 2> {log}"


rule get_protein_names:
    input:
        "protein_annotation/blast_results.m8",
    output:
        protein_names=report(
            "protein_annotation/proteins.csv",
            category="Protein Annotation Results",
            subcategory="Identified Proteins",
            caption="../report/protein_annotation.rst",
            labels={
                "list": "protein_names",
            },
        ),
    log:
        "logs/lambda/get_protein_names.log",
    conda:
        "../envs/biopython.yml"
    script:
        "../scripts/query_uniref_ids.py"
