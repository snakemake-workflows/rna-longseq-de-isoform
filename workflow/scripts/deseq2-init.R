# Start logging
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

library(stringr)
library("DESeq2")

# Enable parallelization if multiple threads are available
parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# Load merged count input data
counts_data <- read.table(
  snakemake@input[["all_counts"]],
  header = TRUE,
  row.names = "Reference",
  check.names = FALSE
)
counts_data <- counts_data[, order(names(counts_data))]

# Load sample metadata
col_data <- read.table(
  snakemake@input[["samples"]],
  header = TRUE,
  row.names = "sample",
  check.names = FALSE
)
col_data <- col_data[order(row.names(col_data)), , drop = FALSE]

# Convert sample metadata columns to deseq2 factors
for (defa in snakemake@config[["deseq2"]][["design_factors"]]) {
  col_data[[defa]] <- factor(col_data[[defa]])
}

# Convert batch effect varaibles into factors
batch_effect <- snakemake@config[["deseq2"]][["batch_effect"]]
for (effect in batch_effect) {
    if (str_length(effect) > 0) {
        col_data[[effect]] <- factor(col_data[[effect]])
    }
}

# Building the model.
design_formula <- snakemake@config[["deseq2"]][["fit_type"]]

# If fit type is empty autp-generate a model
# for all design_factors and batch_effect
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

# Create deseq2 object
dds <- DESeqDataSetFromMatrix(
    countData = counts_data,
    colData = col_data,
    design = as.formula(design_formula)
)

# Remove low count genes
dds <- dds[rowSums(counts(dds)) > snakemake@config[["deseq2"]][["mincount"]], ]

# Run deseq2 normalization
dds <- DESeq(dds, parallel = parallel)

# Save results as .dds
saveRDS(dds, file = snakemake@output[[1]])

# Export normalized counts
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