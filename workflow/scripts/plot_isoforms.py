import sys
import os
import pandas as pd
import subprocess
from multiprocessing import Pool, cpu_count

log_file = open(snakemake.log[0], "w")
sys.stderr = sys.stdout = log_file

de_gene_list = snakemake.input.genes[0]
isoforms_bed = snakemake.input.isob
counts_matrix = snakemake.input.counts_matrix
out_dir = snakemake.output[0]

os.makedirs(out_dir, exist_ok=True)


def get_gene_names(de_gene_list):
    try:
        df = pd.read_csv(de_gene_list, sep="\t")
        if df.empty:
            raise ValueError("Empty gene list file")
        return (gene for gene in df.iloc[:,0])
    except (pd.errors.EmptyDataError, pd.errors.ParserError) as e:
        raise ValueError(f"Failed to parse gene list file: {e}")
    except FileNotFoundError:
        raise FileNotFoundError(f"Gene list file not found: {de_gene_list}")


def run_plot_script(isoforms_bed, counts_matrix, gene_name, out_dir):
    try:
        result = subprocess.run(
            [
                "plot_isoform_usage",
                isoforms_bed,
                counts_matrix,
                gene_name,
                "-o",
                f"{out_dir}/{gene_name}",
            ],
            check=True,
            text=True,
            capture_output=True,
        )
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        return (gene_name, "Success", None)
    except subprocess.CalledProcessError as e:
        error_msg = f"Failed to run plot_isoform_usage for gene {gene_name}: {e}"
        return (gene_name, "Failed", error_msg)


def _run_plot_script_worker(args):
    """Wrapper function for multiprocessing to handle a single gene."""
    isoforms_bed, counts_matrix, gene_name, out_dir = args
    try:
        result = subprocess.run(
            [
                "plot_isoform_usage",
                isoforms_bed,
                counts_matrix,
                gene_name,
                "-o",
                f"{out_dir}/{gene_name}",
            ],
            check=True,
            text=True,
            capture_output=True,
        )
        if result.stderr:
            print(result.stderr)
        return (gene_name, "Success", None)
    except subprocess.CalledProcessError as e:
        error_msg = f"Failed to run plot_isoform_usage for gene {gene_name}: {e}"
        return (gene_name, "Failed", error_msg)


# Get the number of worker from Snakemake's threads setting.
num_workers = snakemake.threads if snakemake.threads > 0 else cpu_count()

# Prepare arguments for each gene
genes = list(get_gene_names(de_gene_list))
args_list = [
    (isoforms_bed, counts_matrix, gene, out_dir) for gene in genes
]

# Run in parallel
with Pool(num_workers) as pool:
    results = pool.map(_run_plot_script_worker, args_list)

# Log results
failed = list()
for gene_name, status, error_msg in results:
    if status == "Success":
        print(f"{gene_name}: Successfully processed")
    else:
        print(f"{gene_name}: {error_msg}", file=sys.stderr)
        failed.append(gene_name)

if failed:
    failed_genes = ", ".join(failed)
    print(f"Failed to process {len(failed)} gene(s): {failed_genes}", file=sys.stderr)

log_file.close()

if failed:
    sys.exit(1)

