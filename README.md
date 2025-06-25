[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://img.shields.io/github/actions/workflow/status/snakemake-workflows/transcriptome-differential-expression/.github%2Fworkflows%2Fmain.yml?branch=main
)](https://github.com/snakemake-workflows/transcriptome-differential-expression/actions?query=branch%3Amain+workflow%3ATests)
[![Conventional Commits](https://img.shields.io/badge/Conventional%20Commits-1.0.0-%23FE5196?logo=conventionalcommits&logoColor=white)](https://conventionalcommits.org)

This workflows performs differenential gene expression analysis on Nanopore long reads using `minimap`, `salmon` and `pyDeseq2`. Optionally, splicing detection and isoform quantification is performed using `flair`. If an un-annotated reference transcriptome has to be used, the workflow is capable of comparing transcripts using `lambda` (a `BLAST` alternative) for annotation (using the UniRef database).







