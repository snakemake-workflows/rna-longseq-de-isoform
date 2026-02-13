localrules:
    genome_to_transcriptome,
    standardize_gff,
    correct_transcriptome,


rule standardize_gff:
    input:
        "references/genomic.gff",
    output:
        temp("references/standardized_genomic.gff"),
    log:
        "logs/agat.log",
    conda:
        "../envs/agat.yml"
    message:
        "Standardizing GFF format for isoform analysis compatibility"
    shell:
        """
        agat_convert_sp_gxf2gxf.pl --gff {input} -o {output} &> {log};
        if [ -f genomic.agat.log ]; then
           cat genomic.agat.log >> {log} && rm genomic.agat.log
        fi
        """


rule genome_to_transcriptome:
    input:
        genome="references/genomic.fa",
        annotation="references/standardized_genomic.gff",
    output:
        transcriptome=temp("transcriptome/transcriptome.fa"),
        index=temp("references/genomic.fa.fai"),
    log:
        "logs/gffread/genome_to_transcriptome.log",
    conda:
        "../envs/gffread.yml"
    threads: 1
    shell:
        """
        gffread -w {output.transcriptome} -g {input.genome} {input.annotation} &> {log}
        """


rule correct_transcriptome:
    input:
        "transcriptome/transcriptome.fa",
    output:
        temp("transcriptome/corrected_transcriptome.fa"),
    log:
        "logs/gffreadcorrect_transcriptome.log",
    threads: 1
    conda:
        "../envs/gffread.yml"
    shell:
        """
        sed 's/ /_/g' {input} > {output} 2> {log}
        """


rule filter_reads:
    input:
        fastq=lambda wildcards: get_mapped_reads_input(
            samples["sample"][wildcards.sample]
        ),
    output:
        temp("filter/{sample}_filtered.fq"),
    message:
        f"Filtering with read length >= {config['read_filter']['min_length']}."
    log:
        "logs/filter_reads/{sample}.log",
    conda:
        "../envs/biopython.yml"
    script:
        "../scripts/read_filter.py"
