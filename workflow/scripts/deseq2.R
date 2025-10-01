log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(cli)
library("DESeq2")
library("pheatmap")
library("RColorBrewer")

parallel <- FALSE
if (snakemake@threads > 1){
    library("BiocParallel")
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

dds <- readRDS(snakemake@input[[1]])

factor <- snakemake@params[["factor"]]
numerator <- snakemake@params[["numerator"]]
denominator <- snakemake@params[["denominator"]]
contrast <- c(factor, numerator, denominator)

res <- results(
    dds,
    contrast = contrast,
    parallel = parallel
)

res <- lfcShrink(
    dds,
    contrast = contrast,
    res = res,
    type = "ashr"
)

res <- res[order(res$padj), ]

svg(snakemake@output[["ma_plot"]])
plotMA(res, ylim = c(-2, 2))
dev.off()

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

vsd <- vst(dds, blind = FALSE)

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)

svg(snakemake@output[["sample_heatmap"]])
pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
)
dev.off()