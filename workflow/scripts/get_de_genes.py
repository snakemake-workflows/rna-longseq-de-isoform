#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from Bio import SeqIO
import pandas as pd

# Start logging
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

sorted_lfc_counts = snakemake.input.sorted_lfc_counts
transcriptome = snakemake.input.transcriptome
gene_list = snakemake.output[0]

# Remove genes with l2fc below lfc threshold from diffexp results
df = pd.read_csv(sorted_lfc_counts)
df.drop(
    df[df["log2FoldChange"].abs() <= snakemake.config["deseq2"]["lfc_null"]].index,
    inplace=True,
)

# Create diffexp gene IDs
gene_names = set(df["Reference"].str.strip())

# Obtain gene records matching diffexp genes
filtered_records = []
for record in SeqIO.parse(transcriptome, "fasta"):
    record_id = record.id.split()[0]
    if record_id in gene_names:
        filtered_records.append(record)

# Error If no matching genes were found
if not filtered_records:
    raise ValueError(
        "No matching gene records found in the transcriptome for differentially expressed genes. "
        "The output FASTA would be empty, which causes lambda to fail."
    )

# Write diffexp genes with sequence to new fasta file
with open(gene_list, "w") as out:
    SeqIO.write(filtered_records, out, "fasta")
