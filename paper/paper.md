---
title: 'A Snakemake workflow for differential expression analyis with alternative splicing detection using long read data'
tags:
  - Snakemake
  - Nanopore
  - HPC
  - differential gene expression
  - alterternative splicing detection
authors:
  - name: Yannic Eising
    orcid: 0009-0003-9103-5689
    affiliation: "1, 2"
  - name: Sören Lukas Hellmann
    orcid: 0000-0003-4958-1419
    affiliation: 2
  - name: Christiane Krämer
    orcid: XXX ?
    affiliation: 2
  - name: Christian Meesters
    corresponding: true
    orcid: 0000-0003-2408-7588
    affiliation: 1
affiliations:
 - name: Nucleic Acids Core Facility, Johannes Gutenberg-University Mainz, Germany
   index: 1
 - name: NHR-SouthWest / High Performance Computing Group, Johannes Gutenberg-University Mainz, Germany
   index: 2
   
date: 04 April 2025 <- update
bibliography: paper.bib

---

# Summary

Long-read RNA sequencing technologies enable the characterization of full-length transcripts and complex splicing patterns. While offering new opportunities for transcriptomic analysis, these data come with substantial computational demands, especially when scaling to multiple samples, replicates, and experimental conditions.

We present a modular, reproducible workflow tailored for differential expression and alternative splicing analysis from long-read RNA sequencing data.
The workflow is designed for use on high-performance computing (HPC) or cloud systems, enabling efficient parallel execution of computationally intensive steps such as read alignment, quantification, and isoform detection.

It supports quality filtering, statistical analysis of gene expression across conditions, and isoform-level splicing analysis. For under-annotated or novel genomes, it includes an optional annotation step based on local similarity searches to assign putative gene functions.
Reference data can be supplied via local files or retrieved automatically using accession numbers.

It is well-suited for researchers working with large datasets and complex experimental designs who require transparent, reproducible, and HPC-compatible analysis pipelines.

# Statement of Need

Long-read sequencing technologies, such as Oxford Nanopore Technologies (ONT), have revolutionized transcriptomic studies by enabling direct detection of full-length RNA molecules [@delahaye_sequencing_2021].
This advancement facilitates more accurate analyses of differential gene expression [@dong_long_2021] and alternative splicing events, both of which are essential for understanding transcriptomic complexity and functional genomics.
However, analyzing long-read transcriptomic data remains technically challenging due to the intricacies of read preprocessing, isoform-level quantification, and the need for reproducible and scalable computational workflows.

Several existing tools—such as FLAIR [@tang_full-length_2020-1] and TALON [@Nicolai2020]—provide frameworks for analyzing long-read transcriptomic data.
While these tools offer powerful features, they often rely on manual configuration, may not fully support reproducible execution across computing environments, and frequently lack integration with high-performance computing (HPC) infrastructure.
Furthermore, in our evaluation of available tools, we identified limitations in the alignment and quantification steps, prompting us to re-implement and optimize the alignment procedures originally provided by FLAIR to improve performance and maintainability within our workflow.

To address these gaps, we present a Snakemake-based workflow that automates the analysis of Nanopore long-read sequencing data with a focus on differential gene expression and alternative splicing detection.
While other workflows exist that support either differential expression analysis or isoform-level analysis, our workflow integrates both in a modular and reproducible pipeline designed for scalability across local machines, HPC clusters, and cloud environments.

A distinctive feature of our workflow is its capability to operate on ill-annotated or completely unannotated genomes.
To support these cases, the workflow includes optional local alignments using tools such as BLAST [@altschul_basic_1990;@camacho_blast_2009] or lambda [@hauswedell_lambda_2014], enabling the functional annotation of transcripts by identifying putative gene functions.
This enhances interpretability in non-model organisms and supports exploratory analyses in less-characterized transcriptomes.

By leveraging Snakemake’s robust workflow management capabilities [@molder_sustainable_2021,
], our pipeline offers transparent provenance tracking, efficient resource handling, and reproducible execution.
It provides a flexible foundation for advanced long-read transcriptomic analyses and fills a critical gap in the ecosystem of accessible, reproducible, and extensible workflows for Nanopore RNA sequencing data.

## Implementation

## Input Data and Reference Handling

The workflow accepts raw ONT reads in FASTQ format, along with either user-specified or automatically downloaded reference data. Reference transcriptomes and genome assemblies can be provided as file paths, or alternatively, specified using NCBI accession numbers, in which case the required data are retrieved using `ncbi-datasets` [@oleary_exploring_2024].
This allows users to flexibly apply the workflow to well-characterized model organisms or newly sequenced, poorly annotated species.

## Quality Filtering and Assessment

Prioar downstream analysis, reads undergo a configurable quality control step. Users can specify minimum average read quality and read length thresholds. For this we make use of the BioPython library [@cock_biopython_2009] Quality statistics and read length distributions are assessed using NanoPlot [@de_coster_nanopack_2018], which generates interactive and publication-ready QC plots. Those are included in the workflow report and ensures high-confidence input for downstream expression and splicing analysis.

## Transcriptome Alignment and Differential Expression Analysis

Reads passing quality filters are aligned to the reference transcriptome minimap2 [@li_minimap2_2018]. Following alignment, read counts per transcript are computed and used for differential expression analysis using pyDESeq2 [@zhu_heavy-tailed_2019;@love_moderated_2014], a Python-native implementation of the DESeq2 method.

This enables statistical analysis of gene expression changes across experimental conditions while staying within a Python-based workflow ecosystem.

## Alternative Splicing Analysis

For isoform-level analysis, the pipeline integrates a modified version of the FLAIR toolkit [@tang_full-length_2020-1]. We adapted key alignment and postprocessing steps to improve compatibility with Snakemake and enhance robustness for HPC use. Isoforms are collapsed, quantified, and categorized to identify splicing patterns and events across conditions.
Optional Functional Annotation via Local Alignment

When reference data are incomplete, unannotated, or of uncertain quality, the workflow offers optional functional annotation. Transcripts or isoforms can be locally aligned against curated protein or nucleotide databases using BLAST or lambda. This provides putative gene product functions that support biological interpretation in non-model organisms or exploratory studies.

ADD rulegraph and caption, here.

# Usage


# Acknowledgements

Any?

# References
