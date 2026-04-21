#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO
import gzip
import shutil

# start logging
sys.stderr = sys.stdout = open(snakemake.log[0], "wt")


def is_zipped(fname):
    return fname.endswith(".gz")


if not isinstance(snakemake.config["read_filter"]["min_length"], int):
    raise TypeError(
        "Your configuration has no valid minimum read length, allowed are: ('int'), check the 'min_length' parameter, please."
    )

if snakemake.config["read_filter"]["min_length"] <= 0:
    if is_zipped(snakemake.input[0]):
        with gzip.open(snakemake.input[0], "rt") as f_in, open(
            snakemake.output[0], "wt"
        ) as f_out:
            shutil.copyfileobj(f_in, f_out)
    else:
        shutil.copy2(snakemake.input[0], snakemake.output[0])
else:
    open_fh = gzip.open if is_zipped(snakemake.input[0]) else open
    with open_fh(snakemake.input[0], "rt") as sample:
        input_iterator = SeqIO.parse(sample, "fastq")

        filter_iterator = (
            read
            for read in input_iterator
            if len(read.seq) >= snakemake.config["read_filter"]["min_length"]
        )

        SeqIO.write(filter_iterator, snakemake.output[0], "fastq")
