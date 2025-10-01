localrules:
    deseq2_init,
    deseq2,


# rule diffexp_analysis:
#     input:
#         all_counts="merged/all_counts.tsv",
#     output:
#         dispersion_graph=report(
#             f"de_analysis/dispersion_graph.{config['deseq2']['figtype']}",
#             category="DGE Results",
#             caption="../report/dispersion_graph.rst",
#             labels={
#                 "figure": "Dispersion graph",
#             },
#         ),
#         ma_graph=report(
#             f"de_analysis/ma_graph.{config['deseq2']['figtype']}",
#             category="DGE Results",
#             caption="../report/ma_graph.rst",
#             labels={
#                 "figure": "MA plot",
#             },
#         ),
#         de_heatmap=report(
#             f"de_analysis/heatmap.{config['deseq2']['figtype']}",
#             category="DGE Results",
#             caption="../report/heatmap.rst",
#             labels={
#                 "figure": "Gene heatmap",
#             },
#         ),
#         correlation_matrix=report(
#             f"de_analysis/correlation_matrix.{config['deseq2']['figtype']}",
#             category="DGE Results",
#             caption="../report/correlation_matrix.rst",
#             labels={
#                 "figure": "Correlation matrix",
#             },
#         ),
#         normalized_counts=report("de_analysis/normalized_counts.csv"),
#         de_top_heatmap=report(
#             f"de_analysis/heatmap_top.{config['deseq2']['figtype']}",
#             category="DGE Results",
#             caption="../report/heatmap_top.rst",
#             labels={
#                 "figure": "Top gene heatmap",
#             },
#         ),
#         sorted_normalized_counts=report("de_analysis/sorted_normalized_counts.csv"),
#         lfc_analysis="de_analysis/lfc_analysis.csv",
#         volcano_plot=report(
#             f"de_analysis/volcano_plot.{config['deseq2']['figtype']}",
#             category="DGE Results",
#             caption="../report/volcano_plot.rst",
#             labels={
#                 "figure": "Volcano plot",
#             },
#         ),
#     params:
#         samples=samples,
#     log:
#         "logs/de_analysis.log",
#     threads: 8
#     conda:
#         "../envs/pydeseq2.yml"
#     script:
#         "../scripts/de_analysis.py"


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
        ma_plot="de_analysis/{factor}_{num}_vs_{den}_MA_plot.svg",
        sample_heatmap="de_analysis/{factor}_{num}_vs_{den}_sample_heatmap.svg",
        count_heatmap="de_analysis/{factor}_{num}_vs_{den}_count_heatmap.svg",
    params:
        factor="{factor}",
        numerator="{num}",
        denominator="{den}",
        colormap=config["deseq2"]["colormap"],
    log:
        "logs/deseq2_{factor}_{num}_vs_{den}.log",
    conda:
        "../envs/deseq2.yml"
    script:
        "../scripts/deseq2.R"
