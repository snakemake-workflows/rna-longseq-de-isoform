rule sam_to_bam:
    input:
        sam="alignments/{sample}.sam",
    output:
        "alignments/{sample}.bam",
    log:
        "logs/samtools/samtobam_{sample}.log",
    params:
        extra=f'{config["samtobam_opts"]}',
    wrapper:
        "v3.13.4/bio/samtools/view"


rule bam_sort:
    input:
        bam="alignments/{sample}.bam",
    output:
        "sorted_alignments/{sample}_sorted.bam",
    log:
        "logs/samtools/bamsort_{sample}.log",
    params:
        extra=f'{config["bamsort_opts"]}',
    wrapper:
        "v3.13.4/bio/samtools/sort"


# flair bam2bed needs index files
rule bam_index:
    input:
        sbam="sorted_alignments/{sample}_sorted.bam",
    output:
        bai="sorted_alignments/{sample}_sorted.bam.bai",
    log:
        "logs/samtools/samindex_{sample}.log",
    params:
        extra=f'{config["bamindex_opts"]}',
    wrapper:
        "v4.5.0/bio/samtools/index"
