# Snakemake workflow: rna-longseq-de-isoform

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://img.shields.io/github/actions/workflow/status/snakemake-workflows/transcriptome-differential-expression/.github%2Fworkflows%2Fmain.yml?branch=main
)](https://github.com/snakemake-workflows/transcriptome-differential-expression/actions?query=branch%3Amain+workflow%3ATests)
[![Conventional Commits](https://img.shields.io/badge/Conventional%20Commits-1.0.0-%23FE5196?logo=conventionalcommits&logoColor=white)](https://conventionalcommits.org)

This workflow performs differential gene expression and isoform splicing analysis on Nanopore long read RNA-Seq data using minimap2, salmon and pyDeseq2. Isoform detection, quantification and differential expression analysis is integrated via FLAIR. We plan to incorporate auto-annotation when ontological data is missing.








