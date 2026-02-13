import os
from pathlib import Path
import re
import sys


log_file = open(snakemake.log[0], "w")

sys.stderr = sys.stdout = log_file

samples_df = snakemake.params.samples[["sample", "condition", "batch"]]

exts = (".fastq", ".fq", ".fastq.gz", ".fq.gz")


def get_sample_path(sample_name, exts):
    # Check for 'raw' directory first, otherwise traverse all directories
    base_path = Path.cwd()
    raw_dir = base_path / "raw"
    sample_regex = re.compile(
        rf"^{re.escape(sample_name)}(?![a-zA-Z0-9]).*({'|'.join([re.escape(ext) for ext in exts])})$"
    )

    search_path = raw_dir if raw_dir.exists() else base_path

    for root, dirs, files in os.walk(search_path):
        for file in files:
            if sample_regex.match(file):
                return os.path.join(root, file)
    raise FileNotFoundError(
        f"No file found for sample '{sample_name}' "
        f"under the search path {search_path} "
        f"with extensions {exts}"
    )


# remove underscores from the name field because flair does not accept them
samples_df["sample_clean"] = samples_df["sample"].str.replace("_", "", regex=False)
# Verify no duplicate sample names were created by the cleaning
if samples_df["sample_clean"].duplicated().any():
    raise ValueError(
        """Exchanging '_' to '' in sample names created duplicates.
           Please ensure original sample names 
           will remain unique after removing underscores.
        """
    )

# get the absolute  filepath for each sample
samples_df["sample_path"] = samples_df["sample"].apply(
    lambda x: get_sample_path(x, exts)
)

samples_df[["sample_clean", "condition", "batch", "sample_path"]].to_csv(
    snakemake.output[0], sep="\t", index=False, header=False
)

log_file.close()
