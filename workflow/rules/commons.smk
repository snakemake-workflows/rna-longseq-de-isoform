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


def rule_all_input():
    all_input = list()
    all_input.append("versions.txt")
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
    all_input.append(f"de_analysis/dispersion_graph.{config['deseq2']['figtype']}")
    all_input.append(f"de_analysis/ma_graph.{config['deseq2']['figtype']}")
    all_input.append(f"de_analysis/heatmap.{config['deseq2']['figtype']}")
    all_input.append("de_analysis/lfc_analysis.csv")
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
