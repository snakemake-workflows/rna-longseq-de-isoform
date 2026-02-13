localrules:
    generate_gene_query,
    get_indexed_protein_db,
    get_protein_names,


rule get_indexed_protein_db:
    output:
        temp("protein_annotation/index/UniRef.lba.gz"),
    params:
        ref=f'{config["protein_annotation"]["uniref"]}',
    log:
        "logs/lambda/get_indexed_db.log",
    conda:
        "../envs/wget.yml"
    shell:
        """
        mkdir -p $(dirname {output}) && \
        wget -nv -O {output} {params.ref} 2> {log}
        """


rule generate_gene_query:
    input:
        sorted_lfc_counts="de_analysis/{factor}_{prop_a}_vs_{prop_b}_l2fc.tsv",
        transcriptome="transcriptome/corrected_transcriptome.fa",
    output:
        temp("protein_annotation/{factor}_{prop_a}_vs_{prop_b}_de_genes.fa"),
    log:
        "logs/lambda/generate_gene_query_{factor}_{prop_a}_vs_{prop_b}.log",
    conda:
        "../envs/biopython.yml"
    script:
        "../scripts/get_de_genes.py"


rule lambda_gene_annotation:
    input:
        indexed_db="protein_annotation/index/UniRef.lba.gz",
        query="protein_annotation/{factor}_{prop_a}_vs_{prop_b}_de_genes.fa",
    output:
        lambda_results=report(
            "protein_annotation/blast_results_{factor}_{prop_a}_vs_{prop_b}.m8",
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
        "logs/lambda/blast_genes_protein_annotation/{factor}_{prop_a}_vs_{prop_b}_de_genes.fa.log",
    conda:
        "../envs/lambda3.yml"
    shell:
        "lambda3 searchp -q {input.query} -i {input.indexed_db} -o {output} -n {params.num_matches} 2> {log}"


rule get_protein_names:
    input:
        "protein_annotation/blast_results_{factor}_{prop_a}_vs_{prop_b}.m8",
    output:
        protein_names=report(
            "protein_annotation/proteins_{factor}_{prop_a}_vs_{prop_b}.csv",
            category="Protein Annotation Results",
            subcategory="Identified Proteins",
            caption="../report/protein_annotation.rst",
            labels={
                "list": "protein_names",
            },
        ),
    log:
        "logs/lambda/get_protein_names_{factor}_{prop_a}_vs_{prop_b}.log",
    conda:
        "../envs/biopython.yml"
    script:
        "../scripts/query_uniref_ids.py"
