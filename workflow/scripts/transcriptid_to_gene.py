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


# Parse gff attributes into a {key: value} dictionairy
def parse_attributes(attr_str):
    attrs = {}
    for part in attr_str.strip().split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            attrs[key] = value
    return attrs


# Build mapping dictionairy {ID: gene}
# Parses gff line by line, extracts attributes field (9th),
# Includes entries that contain both "ID" and "gene" keys
id_to_gene = {}
with open(annotation) as gff:
    for line in gff:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        attrs = parse_attributes(parts[8])
        if "ID" in attrs and "gene" in attrs:
            id_to_gene[attrs["ID"]] = attrs["gene"]


# Clean transcript ID from region-specific suffixes
#   - "gene-Dmel_CR6900" matches directly within annotation
#   - "rna-Dmel_CG34063_CDS=1-1024" must be cleaned first
# The pattern "_CDS=1-1024" is stripped before matching
def clean_id(raw_id):
    return re.sub(r"_CDS=.*", "", raw_id)


# Map a transcript ID to its gene name
# Try exact match using transcript ID
# If not found, try using cleaned ID
# If still not found, return transcript ID as fallback
def map_id(original_ref):
    clean_ref = clean_id(original_ref)
    if original_ref in id_to_gene:
        return f"{id_to_gene[original_ref]}::{original_ref}"
    if clean_ref in id_to_gene:
        return f"{id_to_gene[clean_ref]}::{original_ref}"
    return original_ref


# Replace transcript IDs in "Reference" column with mapping results
counts["Reference"] = counts["Reference"].apply(map_id)

# Write output
counts.to_csv(output, sep="\t", index=False)

log_file.close()
