localrules:
    genome_to_transcriptome,
    standardize_gff,


rule standardize_gff:
    input:
        "references/genomic.gff",
    output:
        "references/standardized_genomic.gff",
    log:
        "logs/agat.log",
    conda:
        "../envs/agat.yml"
    message:
        "Standardizing GFF format for isoform analysis compatibility"
    shell:
        """
        agat_convert_sp_gxf2gxf.pl --gff {input} -o {output} &> {log}
        """


rule genome_to_transcriptome:
    input:
        genome="references/genomic.fa",
        annotation="references/standardized_genomic.gff",
    output:
        transcriptome="transcriptome/transcriptome.fa",
    log:
        "logs/gffread/genome_to_transcriptome.log",
    conda:
        "../envs/gffread.yml"
    threads: 1
    shell:
        """
        gffread -t {threads} -w {output} -g {input.genome} {input.annotation} &> {log}    
        """


rule filter_reads:
    input:
        fastq=lambda wildcards: get_mapped_reads_input(
            samples["sample"][wildcards.sample]
        ),
    output:
        "filter/{sample}_filtered.fq",
    message:
        f"Filtering with read length >= {config['min_length']}."
    log:
        "logs/filter_reads/{sample}.log",
    conda:
        "../envs/biopython.yml"
    script:
        "../scripts/read_filter.py"
