# ==============================================================================
# STEP 0: INITIALIZATION & LIBRARIES
# ==============================================================================
library(WGCNA)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(pheatmap)
library(igraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(STRINGdb)

# Enable multi-threading
enableWGCNAThreads(nThreads = 12)
options(stringsAsFactors = FALSE)
allowWGCNAThreads() 

# CRITICAL FIX: Ensure WGCNA uses its own correlation function globally
cor <- WGCNA::cor 

# ==============================================================================
# STEP 1: LOAD DATA
# ==============================================================================
# setwd("/your/path/here/GSE10072")
load("expression_data_step2_geneLevel.RData")
load("expression_data_step1.RData")

input_mat <- t(expr_matrix)
nSamples <- nrow(input_mat)
rm(expr_matrix) 
gc() 

# ==============================================================================
# STEP 2: SOFT THRESHOLD SELECTION
# ==============================================================================
powers <- c(1:10, seq(12, 20, by=2))
sft <- pickSoftThreshold(input_mat, powerVector=powers, verbose=5, networkType="signed")

if(!dir.exists("plots")) dir.create("plots")
png("plots/WGCNA_SoftThreshold.png", width=1000, height=500)
par(mfrow=c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     type="n", xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, signed R^2", main="Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=0.9, col="red")
abline(h=0.85, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], type="n", xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", main="Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9, col="red")
dev.off()

# ==============================================================================
# STEP 3: NETWORK CONSTRUCTION
# ==============================================================================
picked_power <- 8 # Using your specified power

netwk <- blockwiseModules(input_mat,
                          power=picked_power,
                          networkType="signed",
                          deepSplit=2,
                          minModuleSize=30,
                          maxBlockSize=5000,    
                          mergeCutHeight=0.25,
                          numericLabels=TRUE,
                          pamRespectsDendro=FALSE,
                          saveTOMs=TRUE,        
                          saveTOMFileBase="TOM",
                          verbose=3)

mergedColors <- labels2colors(netwk$colors)

# ==============================================================================
# STEP 4: PLOT DENDROGRAM
# ==============================================================================
pdf("plots/WGCNA_GeneDendrogram.pdf", width=12, height=8)
plotDendroAndColors(netwk$dendrograms[[1]], 
                    mergedColors[netwk$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels=FALSE, hang=0.03,
                    addGuide=TRUE, guideHang=0.05,
                    main="Gene Dendrogram and Module Colors")
dev.off()

# ==============================================================================
# STEP 5: TRAIT CORRELATION
# ==============================================================================
metadata_matched <- metadata[match(rownames(input_mat), rownames(metadata)), ]
head(metadata_matched$title,150)
tumor_trait  <- ifelse(grepl("tumor", metadata_matched$title, ignore.case=TRUE), 1, 0)
normal_trait <- ifelse(grepl("normal", metadata_matched$title, ignore.case=TRUE), 1, 0)

tex_df <- read.csv("TEX_Scores.csv") 
tex_trait <- tex_df$TEX_Score[match(rownames(input_mat), tex_df$SampleID)]

traits <- data.frame(Tumor=tumor_trait, Normal=normal_trait, TEX_Score=tex_trait)
rownames(traits) <- rownames(input_mat)

MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes
module_eigengenes <- orderMEs(MEs0)

moduleTraitCor <- WGCNA::cor(module_eigengenes, traits, use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# ==============================================================================
# STEP 6: HEATMAP (Module-Trait Relationships)
# ==============================================================================
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep="")
dim(textMatrix) <- dim(moduleTraitCor)

png("plots/WGCNA_ModuleTrait_Heatmap.png", width=800, height=1000, res=100)
par(mar=c(6, 8.5, 3, 3))
labeledHeatmap(Matrix=moduleTraitCor, 
               xLabels=colnames(traits),
               yLabels=names(module_eigengenes), 
               ySymbols=names(module_eigengenes), 
               colorLabels=FALSE,
               colors=blueWhiteRed(50), 
               textMatrix=textMatrix,
               setStdMargins=FALSE, 
               cex.text=0.6, 
               zlim=c(-1,1),
               main="Module-Trait Relationships")
dev.off()

# ==============================================================================
# STEP 7: SAVE WORKSPACE
# ==============================================================================
save(input_mat, netwk, mergedColors, traits, module_eigengenes, 
     moduleTraitCor, moduleTraitPvalue, metadata_matched, 
     file = "WGCNA_output.RData")

# ==============================================================================
# STEP 8: AUTOMATIC TOP MODULE IDENTIFICATION & MM vs GS PLOT
# ==============================================================================

# 1. FIND THE TOP MODULE (Dynamic for any dataset)
# Calculate absolute correlation with TEX_Score column in the moduleTraitCor matrix
tex_column <- "TEX_Score" # Ensure this matches your column name exactly
tex_correlations <- abs(moduleTraitCor[, tex_column])

# Automatically find the row name (e.g., "MEblue") of the highest correlation
target_ME <- rownames(moduleTraitCor)[which.max(tex_correlations)]

# Remove the "ME" prefix to get just the color name (e.g., "blue")
module_name <- gsub("ME", "", target_ME)

print(paste("The TOP module for this dataset is:", module_name))

# 2. Get ALL Genes in this Top Module
wgcna_module_genes <- colnames(input_mat)[tolower(mergedColors) == tolower(module_name)]

# 3. Calculate Stats for the Plot
gene_trait_cor <- WGCNA::cor(input_mat, traits$TEX_Score, use = "p")
names(gene_trait_cor) <- colnames(input_mat)

gene_module_membership <- WGCNA::cor(input_mat[, wgcna_module_genes],
                                     module_eigengenes[, target_ME], use = "p")

module_gene_info <- data.frame(
  Gene = wgcna_module_genes,
  ModuleMembership = abs(as.numeric(gene_module_membership)),
  GeneSignificance = abs(as.numeric(gene_trait_cor[wgcna_module_genes])),
  stringsAsFactors = FALSE
)

# 4. Calculate Correlation and P-value for the Plot Legend
fit <- cor.test(module_gene_info$ModuleMembership, module_gene_info$GeneSignificance)
r_value <- round(fit$estimate, 2)
p_value <- signif(fit$p.value, 2)

# 5. Generate the Plot (Replicating Paper Figure 1F)
png(paste0("plots/MM_vs_GS_", module_name, "_PaperStyle.png"), width = 800, height = 600, res = 100)
par(mar = c(5, 5, 4, 2))
plot(module_gene_info$ModuleMembership, 
     module_gene_info$GeneSignificance,
     xlab = paste("Module Membership (MM) in", module_name), 
     ylab = "Gene Significance (GS) for TEX",
     main = paste("MM vs GS in", module_name, "module"), 
     pch = 20, 
     col = module_name, 
     cex = 1.2, cex.lab = 1.2)

# Add Regression Line
abline(lm(GeneSignificance ~ ModuleMembership, data = module_gene_info), col = "red", lwd = 2)

# Add Legend with Stats
legend("bottomright", legend = paste0("cor = ", r_value, "\np = ", p_value),
       bty = "n", cex = 1.2, text.col = "black")
dev.off()

# ==============================================================================
# STEP 9: EXPORT LISTS FOR VENNY & MODULE ASSIGNMENTS
# ==============================================================================

# 1. Get Significant DEGs from your file
deg_data_full <- read.csv("DEGs.csv")
significant_deg_data <- deg_data_full[deg_data_full$P.Value < 0.05 & abs(deg_data_full$logFC) > 0.5, ]
deg_list <- unique(na.omit(as.character(significant_deg_data$Symbol)))

# 2. Create the 2-column list for Venny
max_len <- max(length(wgcna_module_genes), length(deg_list))
output_df <- data.frame(
  WGCNA_Total_Module_Genes = c(wgcna_module_genes, rep(NA, max_len - length(wgcna_module_genes))),
  Significant_DEGs         = c(deg_list,           rep(NA, max_len - length(deg_list)))
)
write.csv(output_df, "Lists_For_Venny_TotalWGCNA_vs_DEGs.csv", row.names = FALSE, na = "")

# 3. Save the full module color assignments file
module_df <- data.frame(Gene = colnames(input_mat), ModuleColor = mergedColors)
write.csv(module_df, "WGCNA_Modules.csv", row.names = FALSE)

cat("\n--- DONE ---\n")
cat("Identified Module:", module_name, "\n")
cat("Module Genes:", length(wgcna_module_genes), "\n")
cat("Significant DEGs:", length(deg_list), "\n")
cat("\n======================================================\n")
cat("                FINAL ANALYSIS SUMMARY                 \n")
cat("======================================================\n")
cat("1. Top Identified Module:          ", module_name, "\n")
cat("2. Module-Trait Correlation:       ", round(max(tex_correlations), 2), "\n")
cat("3. Validation Plot Stats:          cor =", r_value, ", p =", p_value, "\n")
cat("------------------------------------------------------\n")
cat("4. Total Genes in Module:          ", length(wgcna_module_genes), "(Use for Venny Col A)\n")
cat("5. Significant DEG Count:          ", length(deg_list), "(Use for Venny Col B)\n")
cat("------------------------------------------------------\n")
cat("Generated File for Intersection:    Lists_For_Venny_TotalWGCNA_vs_DEGs.csv\n")
