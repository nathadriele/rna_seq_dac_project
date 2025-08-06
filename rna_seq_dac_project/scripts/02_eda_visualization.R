# 02_eda_visualization.R
# -----------------------
# Exploratory Data Analysis (EDA) for RNA-Seq data

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)

dds <- readRDS("data/processed/dds_clean.rds")

vsd <- vst(dds, blind = FALSE)

saveRDS(vsd, "data/processed/vsd_transformed.rds")

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of Samples") +
  theme_minimal()

ggsave("results/pca_plot.png", plot = p, width = 6, height = 5)

topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 30)
mat <- assay(vsd)[topVarGenes, ]
mat <- mat - rowMeans(mat)

annotation <- as.data.frame(colData(vsd)[, "condition", drop = FALSE])

png("results/heatmap_top30_genes.png", width = 1000, height = 1000)
pheatmap(mat,
         annotation_col = annotation,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         fontsize_row = 8,
         main = "Top 30 Variable Genes",
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))
dev.off()
