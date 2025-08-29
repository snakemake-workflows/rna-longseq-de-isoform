 This workflow performs differential expression analysis of RNA-seq data obtained from Oxford Nanopore long-read sequencing technology for the species {{snakemake.config["ref"]["species"]}}. 
 First a transcriptome FASTA is constructed using `gffread <https://github.com/gpertea/gffread>`_. Reads are then mapped to the transcriptome with the long-read optimized alignment tool `minimap2 <https://github.com/lh3/minimap2>`_. 
 Next quantification is performed using `salmon <https://github.com/COMBINE-lab/salmon>`_ before normalization and differential expression analysis are conducted by `PyDESeq2 <https://github.com/owkin/PyDESeq2>`_.
 The workflow can optionally analyze splice-isoforms through integrating the `FLAIR <https://github.com/BrooksLabUCSC/flair>`_ workflow.
 For protein searches based on the identified upregulated transcripts, `lambda <https://github.com/seqan/lambda>`_ is used to align sequences against the UniProt reference.
 Additionaly, `NanoPlot <https://github.com/wdecoster/NanoPlot>`_ is employed to analyze initial sequencing data and `QualiMap <https://github.com/EagleGenomics-cookbooks/QualiMap>`_ is used to evaluate mapping results.
