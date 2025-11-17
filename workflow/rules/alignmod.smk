rule sam_to_bam:
    input:
        sam="alignments/{sample}.sam",
    output:
        temp("alignments/{sample}.bam"),
    log:
        "logs/samtools/samtobam_{sample}.log",
    params:
        extra=f'{config["samtools"]["samtobam_opts"]}',
    wrapper:
        "v7.6.0/bio/samtools/view"


rule bam_sort:
    input:
        bam="alignments/{sample}.bam",
    output:
        temp("sorted_alignments/{sample}_sorted.bam"),
    log:
        "logs/samtools/bamsort_{sample}.log",
    params:
        extra=f'{config["samtools"]["bamsort_opts"]}',
    wrapper:
        "v7.6.0/bio/samtools/sort"


# flair bam2bed needs index files
rule bam_index:
    input:
        sbam="sorted_alignments/{sample}_sorted.bam",
    output:
        bai=temp("sorted_alignments/{sample}_sorted.bam.bai"),
    log:
        "logs/samtools/samindex_{sample}.log",
    params:
        extra=f'{config["samtools"]["bamindex_opts"]}',
    wrapper:
        "v7.6.0/bio/samtools/index"
