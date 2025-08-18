import os
from glob import glob
from pathlib import Path
import sys
from itertools import product

import pandas as pd
from snakemake.remote import FTP
from snakemake.utils import validate
from snakemake.exceptions import WorkflowError

# global list of valid suffixes
exts = (".fastq", ".fq", ".fastq.gz", ".fq.gz")

validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(
        os.path.join(
            os.path.realpath(os.path.dirname(workflow.configfiles[0])),
            config["samples"],
        ),
        sep=r"\s+",
        dtype={"sample": str, "condition": str, "condition2": str, "batch": str},
        header=0,
        comment="#",
    )
    .set_index("sample", drop=False)
    .sort_index()
)
# validation of config files
validate(samples, schema="../schemas/samples.schema.yaml")
validate(config, schema="../schemas/config.schema.yaml")

# flair uses the values of condition1 in its file naming scheme, therefore we extract them as wildcards from samples
condition_val = samples["condition"].unique().tolist()
condition_value1, condition_value2 = condition_val[0], condition_val[1]
condition_samples = {
    cond: samples[samples["condition"] == cond]["sample"].tolist()
    for cond in condition_val
}
if config["isoform_analysis"]["FLAIR"]:
    if len(condition_val) != 2:
        raise ValueError(
            "If you want to perform differential isoform analysis, 'condition' in samples.csv must have exactly two distinct values."
        )


def get_reference_files(config):
    """
    Get reference files from config file and validate them.
    """
    ref = config.get("ref", {})
    genome_exts = (".fa", ".fna", ".fasta")
    annotation_exts = (".gtf", ".gff")
    # Validate genome and annotation files
    genome = (
        ref.get("genome")
        if Path(ref["genome"]).exists()
        and Path(ref["genome"]).suffix.lower() in genome_exts
        else None
    )
    annotation = (
        ref.get("annotation")
        if Path(ref["annotation"]).exists()
        and Path(ref["annotation"]).suffix.lower() in annotation_exts
        else None
    )
    if genome and annotation:
        return {"genome": genome, "annotation": annotation}

    accession = ref.get("accession")
    files = {}
    if genome:
        files["genome"] = genome
    else:
        if accession:
            files["genome"] = "references/ncbi_dataset_a.zip"

    if annotation:
        files["annotation"] = annotation
    else:
        if accession:
            files["annotation"] = "references/ncbi_dataset_b.zip"

    # ValueError: If reference configuration is invalid or missing
    if not files:
        raise ValueError("No valid reference files or accession number provided.")
    return files


def get_mapped_reads_input(sample):
    path = Path(os.path.join(config["inputdir"], sample))
    for extension in exts:
        # we need to append the extension with +, because
        # path.with_suffix might consider everything after a . in
        # the file name a suffix!
        if os.path.exists(str(path) + extension):
            return str(path) + extension

    raise WorkflowError(
        f"No valid sample found for sample: '{sample}' with possible extension '{exts}'"
    )


def aggregate_input(samples):
    # possible extensions:
    valids = list()
    for sample, ext in product(samples, exts):
        path = Path(os.path.join(config["inputdir"], sample))
        # we need to append the extension with +, because
        # path.with_suffix might consider everything after a . in
        # the file name a suffix!
        if os.path.exists(str(path) + ext):
            valids.append(str(path) + ext)

    if not len(valids):
        raise WorkflowError(f"no valid samples found, allowed extensions are: '{exts}'")
    return valids


def get_multifactor_outputs():
    """
    Generate expected output files for multi-factor analysis with unknown value handling.
    """
    outputs = []
    
    # Check if multi-factor analysis is enabled
    if not config["deseq2"].get("multi_factor_analysis", False):
        return outputs
    
    # Get factors to analyze
    factors_to_analyze = config["deseq2"].get("factors_to_analyze", None)
    if factors_to_analyze is None:
        # Auto-detect categorical factors
        continuous_factors = config["deseq2"].get("continuous_factors", [])
        metadata_columns = samples.columns.tolist()
        factors_to_analyze = [col for col in metadata_columns 
                            if col not in continuous_factors and col not in ["sample", "samples"]]
    
    handle_unknowns = config["deseq2"].get("handle_unknown_values", False)
    cross_factor_analysis = config["deseq2"].get("cross_factor_analysis", False)
    figtype = config["deseq2"]["figtype"]
    
    # Generate outputs for individual factors
    for factor in factors_to_analyze:
        if factor in samples.columns:
            if handle_unknowns:
                # Check if factor has unknown values
                unknown_indicators = ["unknown", "not available", "na", "n/a", ""]
                factor_values = samples[factor].astype(str).str.lower()
                has_unknowns = factor_values.isin(unknown_indicators).any()
                
                if has_unknowns:
                    # Get known levels for this factor
                    known_samples = samples[~factor_values.isin(unknown_indicators)]
                    known_levels = known_samples[factor].unique()
                    
                    if len(known_levels) >= 2:
                        # Three scenarios for handling unknowns
                        scenarios = [
                            f"{factor}_unknown_as_{known_levels[0]}",
                            f"{factor}_unknown_as_{known_levels[1]}" if len(known_levels) >= 2 else None,
                            f"{factor}_exclude_unknown"
                        ]
                        scenarios = [s for s in scenarios if s is not None]
                        
                        for scenario in scenarios:
                            # For each scenario, generate outputs for all pairwise comparisons
                            if scenario.endswith("_exclude_unknown"):
                                scenario_levels = known_levels
                            else:
                                scenario_levels = known_levels  # unknowns are reassigned to existing levels
                            
                            for i, level_a in enumerate(scenario_levels):
                                for level_b in scenario_levels[i+1:]:
                                    comparison = f"{factor}_{level_a}_vs_{level_b}"
                                    prefix = f"de_analysis_{scenario}_{comparison}"
                                    outputs.extend([
                                        f"{prefix}_lfc_analysis.csv",
                                        f"{prefix}_ma_plot.{figtype}",
                                        f"{prefix}_volcano.{figtype}"
                                    ])
                else:
                    # No unknowns, standard analysis
                    levels = samples[factor].unique()
                    for i, level_a in enumerate(levels):
                        for level_b in levels[i+1:]:
                            comparison = f"{factor}_{level_a}_vs_{level_b}"
                            prefix = f"de_analysis_{factor}_{comparison}"
                            outputs.extend([
                                f"{prefix}_lfc_analysis.csv",
                                f"{prefix}_ma_plot.{figtype}",
                                f"{prefix}_volcano.{figtype}"
                            ])
            else:
                # Standard multi-factor analysis without unknown handling
                levels = samples[factor].unique()
                for i, level_a in enumerate(levels):
                    for level_b in levels[i+1:]:
                        comparison = f"{factor}_{level_a}_vs_{level_b}"
                        prefix = f"de_analysis_{factor}_{comparison}"
                        outputs.extend([
                            f"{prefix}_lfc_analysis.csv",
                            f"{prefix}_ma_plot.{figtype}",
                            f"{prefix}_volcano.{figtype}"
                        ])
    
    # Generate outputs for cross-factor analysis
    if cross_factor_analysis and len(factors_to_analyze) >= 2:
        from itertools import combinations
        max_factors = config["deseq2"].get("max_cross_factors", 2)
        
        for n_factors in range(2, min(max_factors + 1, len(factors_to_analyze) + 1)):
            for factor_combo in combinations(factors_to_analyze, n_factors):
                if all(f in samples.columns for f in factor_combo):
                    # For cross-factor analysis, we use cleaned data (excluding unknowns)
                    clean_samples = samples.copy()
                    if handle_unknowns:
                        unknown_indicators = ["unknown", "not available", "na", "n/a", ""]
                        for factor in factor_combo:
                            factor_values = clean_samples[factor].astype(str).str.lower()
                            unknown_mask = factor_values.isin(unknown_indicators)
                            clean_samples = clean_samples[~unknown_mask]
                    
                    if len(clean_samples) >= 4:  # Minimum samples needed
                        # Generate combined factor labels
                        combined_labels = []
                        for _, row in clean_samples.iterrows():
                            label_parts = [str(row[factor]) for factor in factor_combo]
                            combined_label = "_".join(label_parts)
                            combined_labels.append(combined_label)
                        
                        unique_labels = list(set(combined_labels))
                        
                        # Generate pairwise comparisons
                        for i, label_a in enumerate(unique_labels):
                            for label_b in unique_labels[i+1:]:
                                factor_name = "_".join(factor_combo)
                                comparison = f"{factor_name}_{label_a}_vs_{label_b}"
                                prefix = f"de_analysis_cross_{factor_name}_{comparison}"
                                outputs.extend([
                                    f"{prefix}_lfc_analysis.csv",
                                    f"{prefix}_ma_plot.{figtype}",
                                    f"{prefix}_volcano.{figtype}"
                                ])
    
    # Add summary file
    if outputs:
        outputs.append("de_analysis_summary.csv")
    
    return outputs


def rule_all_input():
    all_input = list()
    all_input.extend(
        expand("NanoPlot/{sample}/NanoPlot-report.html", sample=samples["sample"])
    )
    all_input.append("NanoPlot/NanoPlot-report.html")
    all_input.extend(expand("QC/bamstats/{sample}.txt", sample=samples["sample"]))
    all_input.extend(
        expand("qualimap/{sample}/qualimapReport.html", sample=samples["sample"])
    )
    all_input.extend(
        expand("counts/{sample}_salmon/quant.sf", sample=samples["sample"])
    )
    all_input.append("merged/all_counts.tsv")
    
    # Standard DE analysis outputs
    all_input.append(f"de_analysis/dispersion_graph.{config['deseq2']['figtype']}")
    all_input.append(f"de_analysis/ma_graph.{config['deseq2']['figtype']}")
    all_input.append(f"de_analysis/heatmap.{config['deseq2']['figtype']}")
    all_input.append("de_analysis/lfc_analysis.csv")
    
    # Multi-factor analysis outputs (including unknown value handling)
    multifactor_outputs = get_multifactor_outputs()
    all_input.extend(multifactor_outputs)
    
    if config["isoform_analysis"]["FLAIR"] == True:
        all_input.extend(
            expand(
                "iso_analysis/diffexp/genes_deseq2_{condition_value1}_v_{condition_value2}.tsv",
                condition_value1=[condition_value1],
                condition_value2=[condition_value2],
            )
        )
        all_input.append("iso_analysis/report/isoforms")
        all_input.append("iso_analysis/report/usage")
    if config["protein_annotation"]["lambda"] == True:
        all_input.append("protein_annotation/proteins.csv")
    return all_input
