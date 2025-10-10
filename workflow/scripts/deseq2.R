# Start logging
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(cli)
library("DESeq2")
library("pheatmap")
library("RColorBrewer")

# Enable parallelization if multiple threads are available
parallel <- FALSE
if (snakemake@threads > 1){
    library("BiocParallel")
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# Load deseq2-init results
dds <- readRDS(snakemake@input[[1]])

# Get contrasts from Snakemake wildcards
factor <- snakemake@params[["factor"]]
numerator <- snakemake@params[["numerator"]]
denominator <- snakemake@params[["denominator"]]
contrast <- c(factor, numerator, denominator)

# Compute diffexp results
res <- results(
    dds,
    contrast = contrast,
    parallel = parallel,
    alpha = snakemake@params[["alpha"]],
    lfcThreshold = as.numeric (snakemake@params[["lfc_null"]]),
    altHypothesis = snakemake@params[["alt_hypothesis"]]
)

# Apply lfc shrinkage to results
res <- lfcShrink(
    dds,
    contrast = contrast,
    res = res,
    type = "ashr",
)

# Sort results by significance
res <- res[order(res$padj), ]

# Generate MA plot
svg(snakemake@output[["ma_plot"]])
plotMA(res, ylim = c(-2, 2))
dev.off()

# Export diffexp results
write.table(
    data.frame(
        gene = rownames(res),
        res
    ),
    file = snakemake@output[["table"]],
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
)

# Variance Stabilizing Transformation
# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#extracting-transformed-values
vsd <- vst(dds, blind = FALSE)
# Obtain sample to smaple dists
sampleDists <- dist(t(assay(vsd)))
# Generate sample dist matrix
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)

# Create sample to sample dist heatmap
svg(snakemake@output[["sample_heatmap"]])
pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  col = colorRampPalette(rev(brewer.pal(9, snakemake@params[["colormap"]])))(255)
)
dev.off()

# Select top 500 genes by condition variance
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 500)

# Create gene heatmap
svg(snakemake@output[["count_heatmap"]])
pheatmap(
    assay(vsd)[topVarGenes, ],
    scale ="row",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    show_rownames = FALSE,
    col = colorRampPalette(rev(brewer.pal(9, snakemake@params[["colormap"]])))(255)
)
dev.off()

# Select top config["deseq2"]["threshold_plot"] by condition variance
toptopVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), as.numeric(snakemake@params[["threshold_plot"]]))

# Create top gene heatmap including gene legend
svg(snakemake@output[["top_count_heatmap"]])
pheatmap(
    assay(vsd)[toptopVarGenes, ],
    scale ="row",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    show_rownames = TRUE,
    col = colorRampPalette(rev(brewer.pal(9, snakemake@params[["colormap"]])))(255)
)
dev.off()

# Create dispersion plot
svg(snakemake@output[["dispersion_plot"]])
plotDispEsts(dds)
dev.off()