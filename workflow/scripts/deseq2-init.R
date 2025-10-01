log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

library(stringr)
library("DESeq2")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

counts_data <- read.table(
  snakemake@input[["all_counts"]],
  header = TRUE,
  row.names = "Reference",
  check.names = FALSE
)
counts_data <- counts_data[, order(names(counts_data))]

col_data <- read.table(
  snakemake@input[["samples"]],
  header = TRUE,
  row.names = "sample",
  check.names = FALSE
)
col_data <- col_data[order(row.names(col_data)), , drop = FALSE]

for (defa in snakemake@config[["deseq2"]][["design_factors"]]) {
  col_data[[defa]] <- factor(col_data[[defa]])
}

batch_effect <- snakemake@config[["deseq2"]][["batch_effect"]]
for (effect in batch_effect) {
    if (str_length(effect) > 0) {
        col_data[[effect]] <- factor(col_data[[effect]])
    }
}

design_formula <- snakemake@config[["deseq2"]][["fit_type"]]

if (str_length(design_formula) == 0) {
  batch_effect <- str_flatten(batch_effect, " + ")
  if (str_length(batch_effect) > 0) {
    batch_effect <- str_c(batch_effect, " + ")
  }
  defa_interactions <- str_flatten(
    snakemake@config[["deseq2"]][["design_factors"]],
    " * "
  )
  design_formula <- str_c("~", batch_effect, defa_interactions)
}

dds <- DESeqDataSetFromMatrix(
    countData = counts_data,
    colData = col_data,
    design = as.formula(design_formula)
)

dds <- dds[rowSums(counts(dds)) > 1, ]

dds <- DESeq(dds, parallel = parallel)

saveRDS(dds, file = snakemake@output[[1]])

norm_counts <- counts(dds, normalized = TRUE)
write.table(
    data.frame(
        "Reference" = rownames(norm_counts),
        norm_counts
    ),
    file = snakemake@output[[2]],
    sep = "\t",
    row.names = FALSE
)