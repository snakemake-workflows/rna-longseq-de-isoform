import os
import sys
from pathlib import Path
from itertools import combinations

import matplotlib

matplotlib.use("Agg")  # suppress creating of interactive plots
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "sans-serif"
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
from bioinfokit import visuz


from snakemake.exceptions import WorkflowError

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

def perform_pairwise_de_analysis(dds, metadata, factor_name, output_prefix, config):
    """
    Perform pairwise differential expression analysis for all combinations 
    of levels in a categorical factor.
    
    Parameters:
    -----------
    dds : DeseqDataSet
        The fitted DESeq2 dataset
    metadata : pd.DataFrame
        Sample metadata
    factor_name : str
        Name of the categorical factor to analyze
    output_prefix : str
        Prefix for output files
    config : dict
        Configuration dictionary
    
    Returns:
    --------
    dict: Dictionary containing results for each comparison
    """
    results = {}
    
    # Get unique levels of the factor
    levels = metadata[factor_name].unique()
    
    if len(levels) < 2:
        raise WorkflowError(f"Factor '{factor_name}' must have at least 2 levels for comparison.")
    
    # Generate all pairwise combinations
    for level_a, level_b in combinations(levels, 2):
        comparison_name = f"{factor_name}_{level_a}_vs_{level_b}"
        print(f"Performing DE analysis for: {comparison_name}")
        
        # Create DeseqStats object for this comparison
        stat_res = DeseqStats(dds, contrast=[factor_name, level_a, level_b])
        
        # Run statistical tests
        if ("lfc_null" in config["deseq2"] or "alt_hypothesis" in config["deseq2"]):
            summary = stat_res.summary(
                lfc_null=config["deseq2"].get("lfc_null", 0),
                alt_hypothesis=config["deseq2"].get("alt_hypothesis", None),
            )
        else:
            summary = stat_res.summary()
        
        # Perform LFC shrinkage
        try:
            stat_res.lfc_shrink(coeff=f"{factor_name}_{level_a}_vs_{level_b}")
        except KeyError:
            try:
                stat_res.lfc_shrink(coeff=f"{factor_name}_{level_b}_vs_{level_a}")
            except KeyError:
                print(f"Warning: Could not perform LFC shrinkage for {comparison_name}")
        
        # Store results
        results[comparison_name] = {
            'stat_res': stat_res,
            'summary': summary,
            'levels': (level_a, level_b)
        }
        
        # Save results to CSV
        output_file = f"{output_prefix}_{comparison_name}_lfc_analysis.csv"
        stat_res.results_df.to_csv(output_file)
        
        # Generate MA plot
        ma_plot_file = f"{output_prefix}_{comparison_name}_ma_plot.{config['deseq2']['figtype']}"
        stat_res.plot_MA(
            s=config["deseq2"]["point_width"],
            save_path=ma_plot_file,
        )
        
        # Generate volcano plot if alpha is stringent enough
        if config["deseq2"]["alpha"] < 0.9:
            volcano_file = f"{output_prefix}_{comparison_name}_volcano.{config['deseq2']['figtype']}"
            visuz.GeneExpression.volcano(
                df=stat_res.results_df.fillna(1),
                lfc="log2FoldChange",
                pv="padj",
                lfc_thr=(
                    config["deseq2"]["lfc_null"],
                    config["deseq2"]["lfc_null"],
                ),
                pv_thr=(
                    config["deseq2"]["alpha"],
                    config["deseq2"]["alpha"],
                ),
                sign_line=True,
                gstyle=2,
                show=False,
                plotlegend=True,
                legendpos="upper right",
                legendanchor=(1.46, 1),
                figtype=config["deseq2"]["figtype"],
            )
            os.rename(
                f"volcano.{config['deseq2']['figtype']}",
                volcano_file,
            )
    
    return results

def generate_factor_combinations(metadata, factors):
    """
    Generate all unique combinations of factor levels across multiple factors.
    
    Parameters:
    -----------
    metadata : pd.DataFrame
        Sample metadata
    factors : list
        List of factor names to combine
        
    Returns:
    --------
    list: List of unique factor combinations found in the data
    """
    # Get all unique combinations that actually exist in the data
    existing_combinations = metadata[factors].drop_duplicates()
    
    # Convert to list of tuples for easier handling
    combinations_list = []
    for _, row in existing_combinations.iterrows():
        combination = tuple(row[factor] for factor in factors)
        combinations_list.append(combination)
    
    return combinations_list, factors

def perform_cross_factor_analysis(dds, metadata, factors, output_prefix, config):
    """
    Perform differential expression analysis for all pairwise combinations
    of factor level combinations (e.g., condition_A+treatment_X vs condition_B+treatment_Y).
    
    Parameters:
    -----------
    dds : DeseqDataSet
        The fitted DESeq2 dataset
    metadata : pd.DataFrame
        Sample metadata
    factors : list
        List of factor names to combine
    output_prefix : str
        Prefix for output files
    config : dict
        Configuration dictionary
    
    Returns:
    --------
    dict: Dictionary containing results for each combination comparison
    """
    results = {}
    
    # Generate all existing combinations
    combinations_list, factor_names = generate_factor_combinations(metadata, factors)
    
    if len(combinations_list) < 2:
        print(f"Warning: Only {len(combinations_list)} unique combinations found for factors {factors}")
        return results
    
    print(f"Found {len(combinations_list)} unique factor combinations: {combinations_list}")
    
    # Create a combined factor column for DESeq2
    combined_factor_name = "_".join(factors)
    combined_metadata = metadata.copy()
    
    # Create combined factor labels
    combined_labels = []
    for _, row in metadata.iterrows():
        label_parts = [str(row[factor]) for factor in factors]
        combined_label = "_".join(label_parts)
        combined_labels.append(combined_label)
    
    combined_metadata[combined_factor_name] = combined_labels
    
    # Get unique combined labels that actually exist
    unique_labels = combined_metadata[combined_factor_name].unique()
    
    # Generate all pairwise comparisons
    for i, label_a in enumerate(unique_labels):
        for label_b in unique_labels[i+1:]:
            comparison_name = f"{combined_factor_name}_{label_a}_vs_{label_b}"
            print(f"Performing cross-factor DE analysis for: {comparison_name}")
            
            try:
                # Create DeseqStats object for this comparison
                # We need to create a new DDS with the combined factor
                temp_dds = DeseqDataSet(
                    counts=dds.counts,
                    metadata=combined_metadata,
                    design_factors=config["deseq2"]["design_factors"] + [combined_factor_name],
                    continuous_factors=config["deseq2"].get("continuous_factors", None),
                    refit_cooks=True,
                    inference=dds.inference,
                    fit_type=config["deseq2"]["fit_type"],
                )
                
                # Fit the model
                temp_dds.deseq2(fit_type=config["deseq2"]["fit_type"])
                
                # Create stats object
                stat_res = DeseqStats(temp_dds, contrast=[combined_factor_name, label_a, label_b])
                
                # Run statistical tests
                if ("lfc_null" in config["deseq2"] or "alt_hypothesis" in config["deseq2"]):
                    summary = stat_res.summary(
                        lfc_null=config["deseq2"].get("lfc_null", 0),
                        alt_hypothesis=config["deseq2"].get("alt_hypothesis", None),
                    )
                else:
                    summary = stat_res.summary()
                
                # Perform LFC shrinkage
                try:
                    stat_res.lfc_shrink(coeff=f"{combined_factor_name}_{label_a}_vs_{label_b}")
                except KeyError:
                    try:
                        stat_res.lfc_shrink(coeff=f"{combined_factor_name}_{label_b}_vs_{label_a}")
                    except KeyError:
                        print(f"Warning: Could not perform LFC shrinkage for {comparison_name}")
                
                # Store results
                results[comparison_name] = {
                    'stat_res': stat_res,
                    'summary': summary,
                    'labels': (label_a, label_b),
                    'factor_combination': factors
                }
                
                # Save results to CSV
                output_file = f"{output_prefix}_{comparison_name}_lfc_analysis.csv"
                stat_res.results_df.to_csv(output_file)
                
                # Generate MA plot
                ma_plot_file = f"{output_prefix}_{comparison_name}_ma_plot.{config['deseq2']['figtype']}"
                stat_res.plot_MA(
                    s=config["deseq2"]["point_width"],
                    save_path=ma_plot_file,
                )
                
                # Generate volcano plot if alpha is stringent enough
                if config["deseq2"]["alpha"] < 0.9:
                    volcano_file = f"{output_prefix}_{comparison_name}_volcano.{config['deseq2']['figtype']}"
                    visuz.GeneExpression.volcano(
                        df=stat_res.results_df.fillna(1),
                        lfc="log2FoldChange",
                        pv="padj",
                        lfc_thr=(
                            config["deseq2"]["lfc_null"],
                            config["deseq2"]["lfc_null"],
                        ),
                        pv_thr=(
                            config["deseq2"]["alpha"],
                            config["deseq2"]["alpha"],
                        ),
                        sign_line=True,
                        gstyle=2,
                        show=False,
                        plotlegend=True,
                        legendpos="upper right",
                        legendanchor=(1.46, 1),
                        figtype=config["deseq2"]["figtype"],
                    )
                    os.rename(
                        f"volcano.{config['deseq2']['figtype']}",
                        volcano_file,
                    )
                    
            except Exception as e:
                print(f"Error in cross-factor analysis for {comparison_name}: {str(e)}")
                continue
    
    return results

def handle_unknown_values(metadata, factor_name):
    """
    Handle unknown/missing values in a factor by creating three analysis scenarios.
    
    Parameters:
    -----------
    metadata : pd.DataFrame
        Sample metadata
    factor_name : str
        Name of the factor to process
        
    Returns:
    --------
    dict: Dictionary with three scenarios for handling unknowns
    """
    unknown_indicators = ["unknown", "not available", "na", "n/a", ""]
    
    # Get the factor values
    factor_values = metadata[factor_name].astype(str).str.lower()
    
    # Identify unknown samples
    unknown_mask = factor_values.isin(unknown_indicators)
    unknown_samples = metadata[unknown_mask]
    known_samples = metadata[~unknown_mask]
    
    if len(unknown_samples) == 0:
        # No unknown values, return original data
        return {f"{factor_name}_original": metadata}
    
    # Get known factor levels
    known_levels = known_samples[factor_name].unique()
    
    if len(known_levels) < 2:
        print(f"Warning: Factor '{factor_name}' has fewer than 2 known levels after removing unknowns")
        return {f"{factor_name}_original": metadata}
    
    scenarios = {}
    
    # Scenario 1: Group unknowns with first factor level
    scenario1 = metadata.copy()
    scenario1.loc[unknown_mask, factor_name] = known_levels[0]
    scenarios[f"{factor_name}_unknown_as_{known_levels[0]}"] = scenario1
    
    # Scenario 2: Group unknowns with second factor level
    if len(known_levels) >= 2:
        scenario2 = metadata.copy()
        scenario2.loc[unknown_mask, factor_name] = known_levels[1]
        scenarios[f"{factor_name}_unknown_as_{known_levels[1]}"] = scenario2
    
    # Scenario 3: Exclude unknowns entirely
    scenario3 = known_samples.copy()
    scenarios[f"{factor_name}_exclude_unknown"] = scenario3
    
    return scenarios

def analyze_factor_with_unknowns(dds, metadata, factor_name, output_prefix, config):
    """
    Analyze a factor handling unknown values with multiple scenarios.
    """
    results = {}
    
    # Get scenarios for handling unknowns
    scenarios = handle_unknown_values(metadata, factor_name)
    
    for scenario_name, scenario_metadata in scenarios.items():
        print(f"\n=== Analyzing {scenario_name} ===")
        
        # Check if we need to create a new DDS for this scenario (if samples were excluded)
        if len(scenario_metadata) != len(metadata):
            # Need to subset the counts and create new DDS
            scenario_samples = scenario_metadata.index.tolist()
            scenario_counts = dds.counts.loc[scenario_samples]
            
            # Create new DDS for this scenario
            scenario_dds = DeseqDataSet(
                counts=scenario_counts,
                metadata=scenario_metadata,
                design_factors=config["deseq2"]["design_factors"],
                continuous_factors=config["deseq2"].get("continuous_factors", None),
                refit_cooks=True,
                inference=dds.inference,
                fit_type=config["deseq2"]["fit_type"],
            )
            scenario_dds.deseq2(fit_type=config["deseq2"]["fit_type"])
        else:
            # Same samples, different factor assignments - can reuse DDS structure
            scenario_dds = DeseqDataSet(
                counts=dds.counts,
                metadata=scenario_metadata,
                design_factors=config["deseq2"]["design_factors"],
                continuous_factors=config["deseq2"].get("continuous_factors", None),
                refit_cooks=True,
                inference=dds.inference,
                fit_type=config["deseq2"]["fit_type"],
            )
            scenario_dds.deseq2(fit_type=config["deseq2"]["fit_type"])
        
        # Perform analysis for this scenario
        scenario_results = perform_pairwise_de_analysis(
            scenario_dds, scenario_metadata, factor_name, 
            f"{output_prefix}_{scenario_name}", config
        )
        results[scenario_name] = scenario_results
        
    return results

def analyze_multiple_factors(dds, metadata, factors, output_prefix, config):
    """
    Analyze multiple categorical factors independently and in combination,
    with special handling for unknown values.
    
    Parameters:
    -----------
    dds : DeseqDataSet
        The fitted DESeq2 dataset
    metadata : pd.DataFrame
        Sample metadata
    factors : list
        List of factor names to analyze
    output_prefix : str
        Prefix for output files
    config : dict
        Configuration dictionary
    
    Returns:
    --------
    dict: Dictionary containing results for each factor and factor combinations
    """
    all_results = {}
    
    # Check if unknown handling is enabled
    handle_unknowns = config["deseq2"].get("handle_unknown_values", False)
    
    # 1. Individual factor analysis with unknown handling
    for factor in factors:
        if factor in metadata.columns:
            if handle_unknowns:
                print(f"\n=== Analyzing factor with unknown handling: {factor} ===")
                factor_results = analyze_factor_with_unknowns(
                    dds, metadata, factor, f"{output_prefix}_{factor}", config
                )
                all_results[factor] = factor_results
            else:
                print(f"\n=== Analyzing individual factor: {factor} ===")
                factor_results = perform_pairwise_de_analysis(
                    dds, metadata, factor, f"{output_prefix}_{factor}", config
                )
                all_results[factor] = {f"{factor}_original": factor_results}
        else:
            print(f"Warning: Factor '{factor}' not found in metadata columns: {metadata.columns.tolist()}")
    
    # 2. Cross-factor combinations (new behavior)
    cross_factor_mode = config["deseq2"].get("cross_factor_analysis", False)
    if cross_factor_mode and len(factors) >= 2:
        print(f"\n=== Performing cross-factor analysis ===")
        
        # Generate all combinations of 2 or more factors
        from itertools import combinations as iter_combinations
        
        # For cross-factor analysis with unknowns, we'll use the "exclude_unknown" scenarios
        if handle_unknowns:
            print("Note: Cross-factor analysis will use 'exclude_unknown' scenarios for factors with unknowns")
        
        # Analyze 2-factor combinations
        for factor_pair in iter_combinations(factors, 2):
            if all(f in metadata.columns for f in factor_pair):
                print(f"\n--- Analyzing factor combination: {factor_pair} ---")
                
                # If handling unknowns, use cleaned metadata for cross-factor analysis
                if handle_unknowns:
                    # Create metadata with unknowns excluded for both factors
                    clean_metadata = metadata.copy()
                    unknown_indicators = ["unknown", "not available", "na", "n/a", ""]
                    
                    for factor in factor_pair:
                        factor_values = clean_metadata[factor].astype(str).str.lower()
                        unknown_mask = factor_values.isin(unknown_indicators)
                        clean_metadata = clean_metadata[~unknown_mask]
                    
                    if len(clean_metadata) < 4:  # Need minimum samples for analysis
                        print(f"Warning: Too few samples after removing unknowns for {factor_pair}")
                        continue
                else:
                    clean_metadata = metadata
                
                cross_results = perform_cross_factor_analysis(
                    dds, clean_metadata, list(factor_pair), 
                    f"{output_prefix}_cross_{'_'.join(factor_pair)}", config
                )
                all_results[f"cross_{'_'.join(factor_pair)}"] = cross_results
        
        # Optionally analyze 3+ factor combinations if requested
        max_factors = config["deseq2"].get("max_cross_factors", 2)
        if max_factors > 2 and len(factors) >= 3:
            for n_factors in range(3, min(max_factors + 1, len(factors) + 1)):
                for factor_combo in iter_combinations(factors, n_factors):
                    if all(f in metadata.columns for f in factor_combo):
                        print(f"\n--- Analyzing {n_factors}-factor combination: {factor_combo} ---")
                        
                        # Clean metadata for multi-factor analysis
                        if handle_unknowns:
                            clean_metadata = metadata.copy()
                            unknown_indicators = ["unknown", "not available", "na", "n/a", ""]
                            
                            for factor in factor_combo:
                                factor_values = clean_metadata[factor].astype(str).str.lower()
                                unknown_mask = factor_values.isin(unknown_indicators)
                                clean_metadata = clean_metadata[~unknown_mask]
                            
                            if len(clean_metadata) < 4:
                                print(f"Warning: Too few samples after removing unknowns for {factor_combo}")
                                continue
                        else:
                            clean_metadata = metadata
                        
                        cross_results = perform_cross_factor_analysis(
                            dds, clean_metadata, list(factor_combo), 
                            f"{output_prefix}_cross_{'_'.join(factor_combo)}", config
                        )
                        all_results[f"cross_{'_'.join(factor_combo)}"] = cross_results

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

ncpus = snakemake.threads
samples = snakemake.params.samples

inference = DefaultInference(n_cpus=ncpus)

metadata = samples.loc[:, samples.columns != "samples"]

counts_df = pd.read_csv(f"{snakemake.input.all_counts}", sep="\t", header=0)
# we have a header line containing "Reference" as attribute, hence the following line
# otherwise, we would add an index row, with which we cannot work
counts_df.set_index("Reference", inplace=True)
counts_df = counts_df.T

# next we filter out counts, with counts lower than 10
genes_to_keep = counts_df.columns[
    counts_df.sum(axis=0) >= snakemake.config["deseq2"]["mincount"]
]
counts_df = counts_df[genes_to_keep]

dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design_factors=snakemake.config["deseq2"]["design_factors"],
    continuous_factors=snakemake.config["deseq2"].get("continuous_factors", None),
    refit_cooks=True,
    inference=inference,
    fit_type=snakemake.config["deseq2"]["fit_type"],
)
# compute normalization factors
# dds.fit_size_factors()
# fit dispersion(s) and LFCs
#  this - fits the size factors
#       - the genewise dispersion
#       - the prior dispersion variance
#       - MAP dispersion
#       - log2fold change
#       - cooks distances
dds.deseq2(fit_type=snakemake.config["deseq2"]["fit_type"])

# fitting genewise dispersions
dds.fit_genewise_dispersions()

dds.plot_dispersions(save_path=f"{snakemake.output.dispersion_graph}")

# compute p-values for our dispersion
# Check if we should perform multi-factor analysis
if snakemake.config["deseq2"].get("multi_factor_analysis", False):
    # Get factors to analyze from config or use all categorical columns
    factors_to_analyze = snakemake.config["deseq2"].get("factors_to_analyze", None)
    if factors_to_analyze is None:
        # Auto-detect categorical factors (excluding continuous factors)
        continuous_factors = snakemake.config["deseq2"].get("continuous_factors", [])
        factors_to_analyze = [col for col in metadata.columns 
                            if col not in continuous_factors and col != "samples"]
    
    print(f"Performing multi-factor analysis for factors: {factors_to_analyze}")
    
    # Perform analysis for multiple factors
    all_results = analyze_multiple_factors(
        dds, metadata, factors_to_analyze, 
        snakemake.output.lfc_analysis.replace('.csv', ''), 
        snakemake.config
    )
    
    # Create a summary report
    summary_results = []
    for factor, factor_results in all_results.items():
        for comparison, result in factor_results.items():
            stat_res = result['stat_res']
            significant_genes = stat_res.results_df[
                stat_res.results_df['padj'] < snakemake.config["deseq2"]["alpha"]
            ]
            summary_results.append({
                'factor': factor,
                'comparison': comparison,
                'total_genes': len(stat_res.results_df),
                'significant_genes': len(significant_genes),
                'upregulated': len(significant_genes[significant_genes['log2FoldChange'] > 0]),
                'downregulated': len(significant_genes[significant_genes['log2FoldChange'] < 0])
            })
    
    # Save summary
    summary_df = pd.DataFrame(summary_results)
    summary_file = snakemake.output.lfc_analysis.replace('.csv', '_summary.csv')
    summary_df.to_csv(summary_file, index=False)
    
    # Use the first comparison for downstream analysis (heatmaps, etc.)
    # Or use the comparison with the most significant genes
    if summary_results:
        best_comparison = max(summary_results, key=lambda x: x['significant_genes'])
        factor_name = best_comparison['factor']
        comparison_name = best_comparison['comparison']
        stat_res = all_results[factor_name][comparison_name]['stat_res']
        
        # Save the best comparison as the main output for compatibility
        stat_res.results_df.to_csv(snakemake.output.lfc_analysis)
        stat_res.plot_MA(
            s=snakemake.config["deseq2"]["point_width"],
            save_path=f"{snakemake.output.ma_graph}",
        )
    else:
        # Fallback: create empty files
        pd.DataFrame().to_csv(snakemake.output.lfc_analysis)
        Path(snakemake.output.ma_graph).touch()

else:
    # Original binary comparison logic
    stat_res = DeseqStats(dds)
    # conditions:
    conditions = list()
    for condition in samples["condition"]:
        if condition not in conditions:
            conditions.append(condition)

    if len(conditions) != 2:
        raise WorkflowError(
            "Only binary conditions are allowed. Make sure your samples.csv only has 2 conditions."
        )
    a_condition = conditions[1]  # this order ensures the the lfc shrink condition is met
    b_condition = conditions[0]

    # run Wald test and plot, perform optional threshold tests, if wanted
    if (
        "lfc_null" in snakemake.config["deseq2"]
        or "alt_hypothesis" in snakemake.config["deseq2"]
    ):
        summary = stat_res.summary(
            lfc_null=snakemake.config["deseq2"].get("lfc_null", 0),
            alt_hypothesis=snakemake.config["deseq2"].get("alt_hypothesis", None),
        )
    else:
        summary = stat_res.summary()
    # performing LFC shrinkage - we try both combination, because, we
    # have no foreknowledge of which conditions comes first
    try:
        stat_res.lfc_shrink(coeff=f"condition_{a_condition}_vs_{b_condition}")
    except KeyError:
        stat_res.lfc_shrink(coeff=f"condition_{b_condition}_vs_{a_condition}")

    stat_res.results_df.to_csv(snakemake.output.lfc_analysis)

    stat_res.plot_MA(
        s=snakemake.config["deseq2"]["point_width"],
        save_path=f"{snakemake.output.ma_graph}",
    )

# create a clustermap, based on normalized counts
# dds_df = dds.to_df()
# ds_df.to_csv('dds_df.csv'
# getting and applying the scaling factors
sf = dds.obsm["size_factors"]

normalized = counts_df.T * sf

# shorthand for log2fold and pvalue columns
log2foldchange = np.abs(stat_res.results_df["log2FoldChange"])
# 'pvalue' is a pandas series, linear, of length(number of aligned loci)
pvalue = stat_res.results_df["pvalue"]
padj = stat_res.results_df["padj"]

normalized = normalized.join(log2foldchange)
# normalized = normalized.join(np.array(pvalue[1]))
normalized = normalized.join(padj)

# delete rows, which do not meet our p-value criterion
# the comparison operator is >= because we drop all values >= our desired alpha
# normalized.drop(normalized[padj >= snakemake.config["deseq2"]["alpha"]].index, inplace=True)
normalized.to_csv(snakemake.output.normalized_counts)

# warning: dropping rows before writing might remove the numerical of the selection column
#          the reason is unclear.
sorted = normalized.copy()

sorted.sort_values(by=["log2FoldChange", "padj"], ascending=False, inplace=True)
sorted.drop(sorted[padj >= snakemake.config["deseq2"]["alpha"]].index, inplace=True)
# throw away these columns
# sorted.drop("log2FoldChange", axis=1, inplace=True)
# sorted.drop("padj", axis=1, inplace=True)
sorted.to_csv(snakemake.output.sorted_normalized_counts)

# throw away these columns
normalized.drop("log2FoldChange", axis=1, inplace=True)
normalized.drop("padj", axis=1, inplace=True)
sorted.drop("log2FoldChange", axis=1, inplace=True)
sorted.drop("padj", axis=1, inplace=True)


normalized.dropna(inplace=True)

# precompute linkages, to prevent missing values crashing the script
row_dism = 1 - sorted.T.corr()
row_linkage = hc.linkage(sp.distance.squareform(row_dism), method="complete")
col_dism = 1 - sorted.corr()
col_linkage = hc.linkage(sp.distance.squareform(col_dism), method="complete")

# we calculate the calculation matrix once, with filling NaNs as 0
correlation_matrix = sorted.corr().fillna(0)
# next, a triangle mask is needed, to avoid redundant plotting of the matrix
# in a square
mask = np.triu(np.ones_like(correlation_matrix))

# TODO: add condition labels (e.g. male/female to the map)
cluster = sns.clustermap(
    correlation_matrix,
    cmap=snakemake.config["deseq2"]["colormap"],
    linewidths=0,
)  # , xticklables = metadata.index.to_list())#, yticklabels = sta)
values = cluster.ax_heatmap.collections[0].get_array().reshape(correlation_matrix.shape)
new_values = np.ma.array(values, mask=mask)
cluster.ax_heatmap.collections[0].set_array(new_values)
cluster.ax_col_dendrogram.set_visible(False)
plt.savefig(snakemake.output.correlation_matrix)

# TODO: add condition labels (e.g. male/female to the map)
sns.clustermap(
    normalized.fillna(0),
    cmap=snakemake.config["deseq2"]["colormap"],
    linewidths=0,
    norm=LogNorm(),
)
plt.savefig(snakemake.output.de_heatmap)

n = snakemake.config["deseq2"]["threshold_plot"]
sns.clustermap(
    sorted.fillna(0).iloc[:n],
    cmap=snakemake.config["deseq2"]["colormap"],
    linewidths=0,
    norm=LogNorm(),
)
plt.savefig(snakemake.output.de_top_heatmap)

# our test case has no significant values
# in our CI test, we have no significant data, hence:
if snakemake.config["deseq2"]["alpha"] < 0.9:
    # Note: stat_res should be defined from either multi-factor or binary analysis above
    if 'stat_res' in locals():
        visuz.GeneExpression.volcano(
            df=stat_res.results_df.fillna(1),
            lfc="log2FoldChange",
            pv="padj",
            lfc_thr=(
                snakemake.config["deseq2"]["lfc_null"],
                snakemake.config["deseq2"]["lfc_null"],
            ),
            pv_thr=(
                snakemake.config["deseq2"]["alpha"],
                snakemake.config["deseq2"]["alpha"],
            ),
            sign_line=True,
            gstyle=2,
            show=False,
            plotlegend=True,
            legendpos="upper right",
            legendanchor=(1.46, 1),
            figtype=snakemake.config["deseq2"]["figtype"],
        )
        os.rename(
            f"volcano.{snakemake.config['deseq2']['figtype']}",
            snakemake.output.volcano_plot,
        )
    else:
        Path(snakemake.output.volcano_plot).touch()
else:
    Path(snakemake.output.volcano_plot).touch()
