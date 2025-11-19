log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library("DESeq2")

# load deseq2 data
dds <- readRDS(snakemake@input[[1]])

# compute regularized log-transformed counts for PCA
counts <- rlog(dds, blind=FALSE)

# Validate that the variable exists in colData
variable <- snakemake@wildcards[["variable"]]
if (!variable %in% colnames(colData(dds))) {
  stop(sprintf("Variable '%s' not found in DESeqDataSet colData. Available variables: %s",
               variable, paste(colnames(colData(dds)), collapse=", ")))
}

svg(snakemake@output[[1]])
plotPCA(counts, intgroup = variable)
dev.off()
