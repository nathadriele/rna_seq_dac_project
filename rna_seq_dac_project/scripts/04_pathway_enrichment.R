# 04_pathway_enrichment.R
# ------------------------
# Functional enrichment using clusterProfiler

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(readr)

res <- readRDS("data/processed/deseq2_results.rds")
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  pull(gene)

entrez_ids <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene = entrez_ids$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                readable = TRUE)

write_csv(as.data.frame(ego), "results/go_enrichment.csv")

barplot <- barplot(ego, showCategory = 10)
ggsave("results/go_barplot.png", barplot, width = 8, height = 5)
