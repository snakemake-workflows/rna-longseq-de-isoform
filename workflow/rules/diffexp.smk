localrules:
    deseq2_init,
    deseq2,


# TODO: add mincount from config to discard loci with fewer counts.
rule deseq2_init:
    input:
        all_counts="merged/all_counts.tsv",
        samples=samples_file,
    output:
        "de_analysis/all.rds",
        "de_analysis/normcounts.tsv",
    log:
        "logs/deseq2-init.log",
    conda:
        "../envs/deseq2.yml"
    script:
        "../scripts/deseq2-init.R"


rule deseq2:
    input:
        "de_analysis/all.rds",
    output:
        table="de_analysis/{factor}_{num}_vs_{den}_l2fc.tsv",
        ma_plot=report(
            f"de_analysis/{{factor}}_{{num}}_vs_{{den}}_MA_plot.svg",
            category="DGE Results",
            caption="../report/ma_graph.rst",
            labels={
                "figure": "MA plot",
            },
        ),
        sample_heatmap=report(
            f"de_analysis/{{factor}}_{{num}}_vs_{{den}}_sample_heatmap.svg",
            category="DGE Results",
            caption="../report/correlation_matrix.rst",
            labels={
                "figure": "Correlation matrix",
            },
        ),
        count_heatmap=report(
            f"de_analysis/{{factor}}_{{num}}_vs_{{den}}_count_heatmap.svg",
            category="DGE Results",
            caption="../report/heatmap.rst",
            labels={
                "figure": "Gene heatmap",
            },
        ),
        top_count_heatmap=report(
            f"de_analysis/{{factor}}_{{num}}_vs_{{den}}_top_count_heatmap.svg",
            category="DGE Results",
            caption="../report/heatmap_top.rst",
            labels={
                "figure": "Top gene heatmap",
            },
        ),
        dispersion_plot=report(
            f"de_analysis/{{factor}}_{{num}}_vs_{{den}}_dispersion_plot.svg",
            category="DGE Results",
            caption="../report/dispersion_graph.rst",
            labels={
                "figure": "Dispersion graph",
            },
        ),
    params:
        factor="{factor}",
        numerator="{num}",
        denominator="{den}",
        colormap=config["deseq2"]["colormap"],
        alpha=config["deseq2"]["alpha"],
        lfc_null=config["deseq2"]["lfc_null"],
        alt_hypothesis=config["deseq2"]["alt_hypothesis"],
        threshold_plot=config["deseq2"]["threshold_plot"],
    log:
        "logs/deseq2_{factor}_{num}_vs_{den}.log",
    conda:
        "../envs/deseq2.yml"
    script:
        "../scripts/deseq2.R"
