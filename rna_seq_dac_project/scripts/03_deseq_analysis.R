library(DESeq2)
library(dplyr)
library(readr)
library(ggplot2)

dds <- readRDS("data/processed/dds_clean.rds")

dds <- DESeq(dds)

res <- results(dds)
res <- lfcShrink(dds, coef = 2, type = "apeglm", res = res)

res_ordered <- res[order(res$padj), ]

write_csv(as.data.frame(res_ordered), "results/deseq2_results.csv")
saveRDS(res_ordered, "data/processed/deseq2_results.rds")

res_df <- as.data.frame(res_ordered)
res_df$gene <- rownames(res_ordered)
res_df$significant <- res_df$padj < 0.05

volcano <- ggplot(res_df, aes(log2FoldChange, -log10(padj), color = significant)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10(p-adj)")

ggsave("results/volcano_plot.png", volcano, width = 6, height = 5)