rule build_minimap_index:  ## build minimap2 index
    input:
        target="transcriptome/corrected_transcriptome.fa",
    output:
        index=temp("index/transcriptome_index.mmi"),
    log:
        "logs/minimap2/index.log",
    threads: 4
    params:
        extra=config["minimap2"]["index_opts"],
    wrapper:
        "v7.6.0/bio/minimap2/index"


# mapping reads with minimap2
rule map_reads:
    input:
        target="index/transcriptome_index.mmi",
        query="filter/{sample}_filtered.fq",
    output:
        temp("alignments/{sample}.sam"),
    log:
        "logs/minimap2/mapping_{sample}.log",
    threads: 32
    params:
        extra=f"-p {config['minimap2']['secondary_score_ratio']} -N {config['minimap2']['maximum_secondary']} {config['minimap2']['opts']}",
    wrapper:
        "v7.6.0/bio/minimap2/aligner"
