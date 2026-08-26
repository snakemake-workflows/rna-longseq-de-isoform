#!/usr/bin/env python
# -*- coding: utf-8 -*-

from os import path
import sys
import pandas as pd
import numpy as np
from functools import reduce

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

dfs = {x: pd.read_csv(x, sep="\t") for x in snakemake.input}

ndfs = []
for x, df in dfs.items():
    # Transform counts to integers:
    df = df.rename(columns={"tname": "Reference","len": "length", "num_reads": "count"})
    df["count"] = np.array(df["count"], dtype=int)
    # Take only non-zero counts:
    df = df[df["count"] > 0]
    df = df[["Reference", "count"]]
    df = df.sort_values(by=["count"], ascending=False)
    name = path.dirname(x).rsplit("/", 1)[1].split("_salmon")[0]
    df = df.rename(columns={"count": name})
    ndfs.append(df)
dfs = ndfs

df_merged = reduce(
    lambda left, right: pd.merge(left, right, on="Reference", how="outer"), dfs
)
if True:
    df_merged = df_merged.fillna(0)

df_merged.to_csv(snakemake.output[0], sep="\t", index=False)
