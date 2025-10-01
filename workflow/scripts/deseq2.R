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