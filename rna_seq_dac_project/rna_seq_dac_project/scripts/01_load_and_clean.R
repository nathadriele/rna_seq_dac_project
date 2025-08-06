library(DESeq2)
library(dplyr)
library(readr)

counts <- read_csv("data/raw/counts_matrix.csv")
metadata <- read_csv("data/raw/sample_metadata.csv")

count_mat <- counts %>%
  column_to_rownames(var = "gene_id") %>%
  as.matrix()

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = metadata,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) > 10, ]

saveRDS(dds, "data/processed/dds_clean.rds")