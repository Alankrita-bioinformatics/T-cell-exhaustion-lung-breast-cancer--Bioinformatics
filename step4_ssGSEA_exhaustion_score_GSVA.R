# Load and install all required libraries
library(GSVA)
library(qusage)
library(org.Hs.eg.db)
library(limma)
library(dplyr)
library(ggplot2)
library(ggforce)
library(tibble)
library(magrittr)
library(illuminaHumanv4.db)
library(GEOquery)

# ===============================================
# STEP 1: Load gene-level expression data
# ===============================================
setwd("/your/path/here/GSE70947")
load("expression_data_step2_geneLevel.RData")  # loads 'expr_matrix'

cat("Expression matrix loaded. Dimensions:", dim(expr_matrix), "\n")

# ===============================================
# STEP 2: Define TEX gene set
# ===============================================
tex_genes <- c(
  "PDCD1","LAG3","HAVCR2","TIGIT","CTLA4","CD39","CD244",
  "TOX","TOX2","EOMES","BATF","IRF4","NFAT","NR4A1","IKZF2",
  "PRDM1","EGR2","TCF7","RORC","CXCL13","LAYN","CXCR5","SLAMF6",
  "TNFRSF14","GZMB","GZMA","PRF1","IFNG","TNF","IL2","LEF1","BACH2",
  "IL7R","CCR7","BCL6","SELL","TNFRSF9","ICOS","CD28","CD27","DNMT3A",
  "DNMT3B","UHRF1","KMT2D","HDAC9","SLC7A5","SLC7A1","LDHA","PFKFB3",
  "ITGA4","ITGB1","CDH1","CD3E","CD8A","CD8B","ZAP70","BCL2","MCL1","FAS",
  "CRTAM","CLEC2D","KLRB1"
)

tex_list <- list(T_Cell_Exhaustion = intersect(tex_genes, rownames(expr_matrix)))
cat("TEX genes available in dataset:", length(tex_list[[1]]), "\n")

# ===============================================
# STEP 5: Run ssGSEA (GSVA)
# ===============================================
param <- ssgseaParam(exprData = expr_matrix, geneSets = tex_list)
gsva_results <- gsva(param, verbose = TRUE)

tex_df <- data.frame(
  SampleID = colnames(gsva_results),
  TEX_Score = as.numeric(gsva_results[1, ])
)

# ===============================================
# STEP 6: MANUAL SUBTYPE SELECTION 
# ===============================================
cat("Fetching phenotype data from GEO...\n")
gset <- getGEO("GSE70947", GSEMatrix = TRUE, getGPL = FALSE)
pheno_data <- pData(gset[[1]])

# --- CHECK LABELS HERE ---
# This prints the unique titles so you can see exactly what to type below
print(unique(pheno_data$title))

# --- CHOOSE YOUR LABELS ---
# Replace the text inside the quotes with the exact labels you see in the printout above
tumor_label  <- "tumor"   
normal_label <- "normal"  

subtype_mapping <- data.frame(
  SampleID = rownames(pheno_data),
  # Using an exact match based on your choices above
  Subtype = ifelse(grepl(normal_label, pheno_data$title, ignore.case = TRUE), "Normal", 
                   ifelse(grepl(tumor_label, pheno_data$title, ignore.case = TRUE), "Tumor", "Other"))
)

# Join and filter out 'Other' if necessary
tex_df2 <- tex_df %>%
  left_join(subtype_mapping, by = "SampleID") %>%
  filter(Subtype != "Other")

cat("\nFinal Sample distribution:\n")
print(table(tex_df2$Subtype))

# ===============================================
# STEP 7: CREATE THE SINA PLOT
# ===============================================
tex_df2$Subtype <- factor(tex_df2$Subtype, levels = c("Normal", "Tumor"))

sina_plot <- ggplot(tex_df2, aes(x = Subtype, y = TEX_Score, color = Subtype)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, fatten = 2) +
  geom_sina(alpha = 0.6, size = 2.5) +
  scale_color_manual(values = c("Normal" = "#3182bd", "Tumor" = "#de2d26")) +
  theme_bw() +
  labs(title = "T Cell Exhaustion Scores: Tumor vs Normal (GSE70947)",
       x = "Sample Group", y = "TEX Score", color = "Group") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

ggsave("TEX_Sina_Plot_GSE70947.png", plot = sina_plot, width = 10, height = 8, dpi = 600)
cat("Sina Plot saved successfully!\n")
