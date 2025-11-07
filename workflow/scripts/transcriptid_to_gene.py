import re
import sys
import pandas as pd

# Start logging
log_file = open(snakemake.log[0], "w")
sys.stderr = sys.stdout = log_file


counts = snakemake.input.all_counts
annotation = snakemake.input.annotation
output = snakemake.output[0]


# Load count data
counts = pd.read_csv(counts, sep="\t", dtype=str)


# Parse gff attributes
def parse_attributes(attr_str):
    attrs = {}
    for part in attr_str.strip().split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            attrs[key] = value
    return attrs


# Build mapping dictionairy {ID: gene}
id_to_gene = {}
with open(annotation) as gff:
    for line in gff:
        if line.startswith("#"):
            continue
        # The 9th position in gff3 contains attributes
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        attrs = parse_attributes(parts[8])
        if "ID" in attrs and "gene" in attrs:
            id_to_gene[attrs["ID"]] = attrs["gene"]


# Clean transcript ID for matching
def clean_id(raw_id):
    return re.sub(r"_CDS=.*", "", raw_id)


# Map transcript IDs to gene names
def map_id(original_ref):
    clean_ref = clean_id(original_ref)
    # Look for matches
    if original_ref in id_to_gene:
        return f"{id_to_gene[original_ref]}::{original_ref}"
    if clean_ref in id_to_gene:
        return f"{id_to_gene[clean_ref]}::{original_ref}"
    # Fallback: return original
    return original_ref


# Replace transcript IDs with ID mapping results
counts["Reference"] = counts["Reference"].apply(map_id)

# Write output
counts.to_csv(output, sep="\t", index=False)

log_file.close()
