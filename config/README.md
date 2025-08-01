# Transcriptome Differential Expression Workflow

This workflow facilitates the analysis of transcriptomic data obtained from Oxford Nanopore long-read sequencing, including alignment, quantification, differential expression analysis, alternative splicing analysis and quality control. Users should be aware that alignment parameters may need to be adjusted accordingly for optimal performance with alternative sequencing techniques.

## Configuration Files

To set up the workflow, modify the following files to reflect your dataset and analysis parameters:

- `config/samples.csv`: Contains sample information and experimental design.
- `config/config.yaml`: General workflow configuration and parameter settings.

## samples.csv

Each line in `samples.csv` represents a biological sample with associated metadata. The required columns are:
- `sample`: Unique identifies that matches the sample file in the input directory
- `condition`: Experimental condition or treatment group
- `batch`: Batch of samples
and the following columns forward optional metadata:
- `platform`: Sequencing platform
- `purity`: Purity value

## config.yml

The `config.yml` file contains the main configuration parameters for the workflow.

### General Workflow Parameters

- `inputdir`: Directory containing input samples.

### Reference Genome Parameters

Since Salmon requires transcriptomic alignments for quantification, a transcriptome is constructed using genomic and annotation reference data, these reference files can be provided locally or automatically retrieved from NCBI using an accession number.

- **ref**:
  - `species`: Name of the species.
  - `genome`: Path to the genome file (FASTA format, can be omitted if remote retrieval through accession number is prefered).
  - `annotation`: Path to the annotation file (GFF or GTF format, can be omitted if remote retrieval through accession number is prefered).
  - `accession`: NCBI accession number (used if local files are not provided).

### Read Filtering

As this workflow is designed for long-read sequencing, a custom Python script is used to filter out short reads that may represent contamination or sequencing artifacts, ensuring that only reads meeting the specified length threshold are used for analysis.

- **read_filter**:
  - `min_length`: Minimum read length to retain. Can be left at 0 to consider all reads.

### Alignment (minimap2)

Alignment performed using Minimap2. A comprehensive explanation of its parameters can be found in the [Minimap2 documentation](https://lh3.github.io/minimap2/minimap2.html#10).

- **minimap2**:
  - `index_opts`: Used to define additional options for indexing.
  - `opts`: Used to define additional mapping options.
  - `maximum_secondary`: Maximum number of secondary alignments, `-N` within the minimap2 documentation.
  - `secondary_score_ratio`: Score ratio for secondary alignments, `-p` within the minimap2 documentation.

### Alignment Processing (Samtools)

Since multiple downstream analysis tools require different Alignment formats, Samtools is used to convert SAM files into BAM format for Salmon quantification or sorted SAM files used for FLAIR splice-isoform analysis. Additionally, Samtools generates alignment statistics that serve as quality control metrics. More details can be found in the [Samtools documentation](http://www.htslib.org/doc/samtools.html).

- **samtools**:
  - `samtobam_opts`: Additional options for SAM to BAM conversion.
  - `bamsort_opts`: Additional options for sorting BAM files.
  - `bamindex_opts`: Additional options for indexing sorted BAM files.
  - `bamstats_opts`: Additional options for generating Alignment statistics.

### Quantification (Salmon)

Transcripts are quantified using Salmon in alignment-based mode. TO ensure accurate quantification, Salmon requires information about the strandedness if the sequencing reads. More details can be found in the [Salmon documentation](https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype).

- **quant**:
  - `salmon_libtype`: Library type for Salmon quantification.

### Differential Expression Analysis (DESeq2)

Differential expression analysis is performed using PyDESeq2 to model raw read counts wtih a negative binomial distribution, estimating dispersion parameters to identify differentially expressed genes. See the [PyDESeq2 documentation](https://pydeseq2.readthedocs.io/en/stable/index.html) for more details.

- **deseq2**:
  - `fit_type`: Type of fitting of dispersions to the mean intensity. `parametric`: fit a dispersion-mean relation via a robust gamma-family GLM. `mean`: use the mean of gene-wise dispersion estimates. Will set the fit type for the DEA and the vst transformation. If needed, it can be set separately for each method.
  - `design_factors`: List of design factors for the analysis.
  - `lfc_null`: The (log2) log fold change under the null hypothesis for Wald test.
  - `alt_hypothesis`: The alternative hypothesis for computing wald p-values. By default, the normal Wald test assesses deviation of the estimated log fold change from the null hypothesis, as given by `lfc_null`. The alternative hypothesis corresponds to what the user wants to find rather than the null hypothesis.
  - `point_width`: Marker size for MA-plot
  - `mincount`: Minimum count threshold, genes below the threshold will be removed from analysis.
  - `alpha`: Type I error cutoff value.
  - `threshold_plot`: Number of top differentially expressed genes to plot in additional heatmap.
  - `colormap`: Colormap for heatmaps.
  - `figtype`: Figure output format (e.g., `png`).

### Isoform Analysis (FLAIR)

FLAIR is used to identify alternative splice isoforms in full-length transcripts obtained from long-read sequencing. It then quantifies these transcripts and performs differential expression analysis on the corresponding genes with splice-isoforms. More details can be found in the [FLAIR documentation](https://flair.readthedocs.io/en/latest/index.html).

- **isoform_analysis**:
  - `FLAIR`: Enable FLAIR alternative isoform analysis (`true` or `false`).
  - `qscore`: Minimum MAPQ for read alignment. `--quality` for FLAIR modules.
  - `exp_thresh`: Minimum read count expression threshold for differential expression analysis. Genes with less counts will be removed form analysis.
  - `col_opts`: Additional options for flair collapse module.

### Protein Annotation (lambda)

Lambda aligns sequences of differentially expressed genes or transcripts against indexed protein databases (e.g. UniProt). This process is similar to BLAST, enabling identification of similar proteins and functional annotation of transcripts.

- **Protein Annotation**:
  - `lambda`: Enable lambda Sequence alignment (`true` or `false`).
  - `uniref`: URL of the indexed UniRef database from the lambda wiki.
  - `num_matches`: Maximum number of proteins that have been identififed per sequence.

