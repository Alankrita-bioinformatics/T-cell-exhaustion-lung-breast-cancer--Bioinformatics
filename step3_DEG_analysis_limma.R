# ---------------------------------------------------------
# Step 3 & 4: Differential Expression + Annotation
# ---------------------------------------------------------
# Make sure to install the required libraries first

library(limma)
library(openxlsx)
library(hgu133plus2.db) 
library(ggplot2)
library(ggrepel)
library(dplyr)

# 1. Load Data
# setwd("/your/path/here")
load("expression_data_step1.RData")
load("expression_data_step2_geneLevel.RData")

# 2. Design & Contrast
design <- model.matrix(~0 + subtype)
colnames(design) <- levels(subtype)
colnames(design)
fit <- lmFit(expression_data, design)
cont.matrix <- makeContrasts(Tumor - Normal, levels=design)
#contrast_name <- colnames(cont.matrix)[1]
#head(contrast_name)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# ---------------------------------------------------------
# Step 3: SAFE Class Labeling (Works for any dataset)
# ---------------------------------------------------------

# 1. Run the Limits / Bayes Statistics
all_genes <- topTable(fit2, coef=1, adjust="fdr", number=Inf) 

# 2. Extract the comparison names from YOUR contrast matrix
# Example output: "Tumor - Normal"
contrast_name <- colnames(cont.matrix)[1] 
print(paste("Your Contrast is:", contrast_name))

# Split the string to find out which is (+) and which is (-)
# "Tumor - Normal" -> Group1="Tumor" (Positive), Group2="Normal" (Negative)
groups <- unlist(strsplit(contrast_name, " - "))
positive_group <- groups[1]  # The one on the Left
negative_group <- groups[2]  # The one on the Right

# 3. Assign Class labels dynamically
all_genes$Class <- ifelse(all_genes$logFC > 0, positive_group, negative_group)

# 4. Save
write.csv(all_genes, "DEGs.csv")

cat("Success! Labeled genes based on contrast:", contrast_name, "\n")
cat("Positive LogFC (>0) = ", positive_group, "\n")
cat("Negative LogFC (<0) = ", negative_group, "\n")

# 5. Prepare Plot Data (Labeling)
#    Define Thresholds
PVAL_CUTOFF <- 0.05
LOGFC_CUTOFF <- 0.5  # Changed to 1.0 to reduce list size (0.5 is often too loose)

plot_data <- all_genes
plot_data$Significant <- "Not Significant"
plot_data$Significant[plot_data$adj.P.Val < PVAL_CUTOFF & plot_data$logFC > LOGFC_CUTOFF] <- "Upregulated"
plot_data$Significant[plot_data$adj.P.Val < PVAL_CUTOFF & plot_data$logFC < -LOGFC_CUTOFF] <- "Downregulated"

# 6. VOLCANO PLOT
top_genes <- plot_data %>%
  arrange(adj.P.Val) %>%
  head(10)
top_genes$ID <- rownames(top_genes)

volcano_plot <- ggplot(plot_data, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Downregulated" = "blue", "Not Significant" = "grey", "Upregulated" = "red")) +
  geom_vline(xintercept = c(-LOGFC_CUTOFF, LOGFC_CUTOFF), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(PVAL_CUTOFF), linetype = "dashed", color = "black") +
  geom_label_repel(data = top_genes, aes(label = ID), box.padding = 0.5, max.overlaps = Inf) +
  labs(title = "Volcano Plot: Tumor vs Control", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_minimal()

ggsave("Volcano_Plot_Colorful.png", plot = volcano_plot, width = 8, height = 6)

# 7. MA PLOT
ma_plot <- ggplot(plot_data, aes(x = AveExpr, y = logFC, color = Significant)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c("Downregulated" = "blue", "Not Significant" = "grey", "Upregulated" = "red")) +
  geom_hline(yintercept = 0, color = "black", linetype="dashed") +
  labs(title = "MA Plot: Tumor vs Control", x = "Average Expression", y = "Log2 Fold Change") +
  theme_minimal()

ggsave("MA_Plot_Colorful.png", plot = ma_plot, width = 8, height = 6)
print("Created plots.")

# 8. SAVE STRICT SIGNIFICANT LIST
#    Filter specifically for P < 0.05 AND |LogFC| > 1.0
significant_strict <- plot_data[plot_data$Significant != "Not Significant", ]
write.csv(significant_strict, "Significant.csv")

cat("Saved Significant.csv with", nrow(significant_strict), "genes.\n")