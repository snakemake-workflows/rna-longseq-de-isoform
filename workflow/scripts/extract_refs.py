import sys
import os
import zipfile
import glob
import shutil

log_file = open(snakemake.log[0], "w")
sys.stderr = sys.stdout = log_file

ref = snakemake.input[0]
output = snakemake.output[0]
accession = snakemake.params.accession

genome_exts = {".fa", ".fna", ".fasta"}
annotation_exts = {".gff", ".gtf", "gff3"}

is_genome = output.endswith(".fa")
is_annotation = output.endswith(".gff")


def extract_ref(ref_path, out_path, ext, temp_suffix):
    temp_dir = f"temp_{accession}_{temp_suffix}"
    with zipfile.ZipFile(ref_path, "r") as zf:
        zf.extractall(temp_dir)

    for fname in glob.glob(f"{temp_dir}/ncbi_dataset/data/{accession}/*{ext}"):
        shutil.copyfile(fname, out_path)
        shutil.rmtree(temp_dir)
        return
    shutil.rmtree(temp_dir)
    raise FileNotFoundError(
        f"No file with {ext} in {ref_path} for accession number: {accession}"
    )


if is_genome and any(ref.endswith(ext) for ext in genome_exts):
    shutil.copyfile(ref, output)
elif is_genome:
    extract_ref(ref, output, ".fna", "genome")

if is_annotation and any(ref.endswith(ext) for ext in annotation_exts):
    shutil.copyfile(ref, output)
elif is_annotation:
    extract_ref(ref, output, ".gff", "annotation")
