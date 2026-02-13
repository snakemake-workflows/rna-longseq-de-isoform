localrules:
    pathway_annotation,


if config["enrichment"]["perform_enrichment"] == True:

    rule pathway_annotation:
        input:
            # Collect all DE result tables (no wildcards) by expanding the
            # precomputed `contrasts` list from `commons.smk`.
            tsvfile=[
                expand("de_analysis/{factor}_{prop_a}_vs_{prop_b}_l2fc.tsv", **c)[0]
                for c in contrasts
            ],
        output:
            legend="pathway_images/pathway_legend.svg",
            plots=directory("pathway_images"),
            list_of_files="pathway_images/list_of_files.txt",
        conda:
            "../envs/colorize_pathways.yml"
        params:
            species=config["ref"]["species"],
            min_genes=config["enrichment"]["min_genes"],
        log:
            "logs/enrichment/pathway_annotation.log",
        threads: 1  # this script will remain single-threaded
        script:
            "../scripts/colorize_pathways.R"
