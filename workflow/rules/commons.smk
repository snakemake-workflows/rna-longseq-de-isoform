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

# full path to samples.csv relative to config.yml
samples_file = os.path.join(
    os.path.dirname(os.path.realpath(workflow.configfiles[0])),  # config.yaml dir
    config["samples"],  # 'samples.csv'
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
    annotation_exts = (".gtf", ".gff", ".gff3")
    # Validate genome and annotation files
    genome_path = ref.get("genome")
    genome = (
        genome_path
        if genome_path
        and Path(genome_path).exists()
        and Path(genome_path).suffix.lower() in genome_exts
        else None
    )
    annotation_path = ref.get("annotation")
    annotation = (
        annotation_path
        if annotation_path
        and Path(annotation_path).exists()
        and Path(annotation_path).suffix.lower() in annotation_exts
        else None
    )
    if genome and annotation:
        return {"genome": genome, "annotation": annotation}

    accession = ref.get("accession")
    ensembl_species = ref.get("ensembl_species")
    files = {}
    if genome:
        files["genome"] = genome
    elif ensembl_species:
        files["genome"] = "references/ensembl_genome.fa"
    elif accession:
        files["genome"] = "references/ncbi_dataset_genome.zip"
    else:
        raise ValueError(
            "No valid genome source: provide local genome path, Ensembl parameters, or NCBI accession."
        )

    if annotation:
        files["annotation"] = annotation
    elif ensembl_species:
        files["annotation"] = "references/ensembl_annotation.gff3"
    elif accession:
        files["annotation"] = "references/ncbi_dataset_annotation.zip"
    else:
        raise ValueError(
            "No valid annotation source: provide local annotation path, Ensembl parameters, or NCBI accession."
        )

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


# Obtain all pairwise contrasts from samples.csv for deseq2
# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#contrasts
contrasts = []
for factor in config["deseq2"]["design_factors"]:
    # Extract all unique conditions for a factor
    levels = samples[factor].unique().tolist()
    # For each combination of two conditions create a contrast dict
    # e.g.: [{"factor": "condition", "prop_a": "male", "prop_b": "female"}]
    for i in range(len(levels)):
        for j in range(i + 1, len(levels)):
            contrasts.append(
                {"factor": factor, "prop_a": levels[i], "prop_b": levels[j]}
            )


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
    all_input.extend(
        [
            expand("de_analysis/{factor}_{prop_a}_vs_{prop_b}_l2fc.tsv", **c)[0]
            for c in contrasts
        ]
    )
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
        all_input.extend(
            [
                expand(
                    "protein_annotation/proteins_{factor}_{prop_a}_vs_{prop_b}.csv", **c
                )[0]
                for c in contrasts
            ]
        )
    return all_input
