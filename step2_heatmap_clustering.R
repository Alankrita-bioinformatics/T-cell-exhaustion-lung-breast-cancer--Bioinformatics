library(pheatmap)
library(ggplot2)

# 1. Set your working directory (if needed)
# setwd("/your/path/here")
image_dir <- getwd() 
save_path <- "/your/path/here/heatmap_top_50_genes.png"
# 2. Load your save file
load("expression_data_step1.RData")

# 3. Check that everything is back
ls() 
# You will see "expression_data", "metadata", "subtype" listed.
# You can now continue immediately to the next step.

# Heatmap of top 100 most variable genes

top_var_genes <- order(apply(expression_data, 1, var), decreasing = TRUE)[1:50]
heatmap_data <- expression_data[top_var_genes, ]
annotation_col <- data.frame(Subtype = subtype)
rownames(annotation_col) <- colnames(expression_data)

pheatmap_plot <- pheatmap(heatmap_data,
                          annotation_col = annotation_col,
                          show_rownames = FALSE,
                          show_colnames = FALSE,
                          main = "Heatmap of Top 50 Most Variable Genes",
                          filename = save_path)
png("heatmap_top_50_genes.png", width = 10, height = 8, units = "in", res = 300)
print(pheatmap_plot)
dev.off()
# Hierarchical clustering
dist_matrix <- dist(t(expression_data))
hclust_result <- hclust(dist_matrix, method = "complete")
png("hierarchical_clustering_dendrogram.png", width = 12, height = 8, units = "in", res = 300)
plot(hclust_result, labels = colnames(expression_data), cex = 0.8, las = 2)
dev.off()

