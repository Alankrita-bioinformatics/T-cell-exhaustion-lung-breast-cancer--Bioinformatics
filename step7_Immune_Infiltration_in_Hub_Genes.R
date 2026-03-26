# ==============================================================================
# PER-DATASET TISSUE ANALYSIS (INDIVIDUAL GEO DATASET IMMUNE INFILTRATION)
# ==============================================================================
# 0. Install following libraries and run only once.
# if (!require("GSVA")) BiocManager::install("GSVA")
# if (!require("pheatmap")) install.packages("pheatmap")
# if (!require("corrplot")) install.packages("corrplot")
# if (!require("ggplot2")) install.packages("ggplot2")
# if (!require("reshape2")) install.packages("reshape2")
# if (!require("Hmisc")) install.packages("Hmisc") 

# 1. SETUP LIBRARIES
library(GSVA)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(corrplot)
library(Hmisc)

# ==============================================================================
# 2. USER CONFIGURATION
# ==============================================================================
lung_path   <- "/your/path/here/lung_data/"
breast_path <- "/your/path/here/breast_data/"

hub_genes <- c("CD8A", "CD4", "CCL5", "FCGR2B", "FCGR2A", 
               "CD38", "CCR7", "PDCD1", "CTLA4", "HAVCR2")

cell_markers <- list(
  Activated_B_cell = c("CD79A", "CD79B", "CD19"),
  Activated_CD4_T_cell = c("CD4", "IL7R", "CD40LG"),
  Activated_CD8_T_cell = c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1"),
  Activated_dendritic_cell = c("CD80", "CD83", "CD86", "LAMP3"),
  CD56dim_natural_killer_cell = c("KIR2DL3", "KIR3DL1", "KIR3DL2"),
  CD56bright_natural_killer_cell = c("NCAM1", "XCL1", "XCL2"),
  Central_memory_CD4_T_cell = c("CCR7", "SELL", "CD27"),
  Central_memory_CD8_T_cell = c("CCR7", "SELL", "IL7R"),
  Effector_memeory_CD8_T_cell = c("CX3CR1", "FGFBP2", "FCGR3A"),
  Effector_memeory_CD4_T_cell = c("CCR4", "CCR6", "CXCR3"),
  Eosinophil = c("CCR3", "IL5RA", "SIGLEC8"),
  Gamma_delta_T_cell = c("TRGC1", "TRGC2", "TRDC"),
  Immature_B_cell = c("CD19", "MS4A1", "CD22"),
  Immature_dendritic_cell = c("CD1A", "CD1B", "CD1C"),
  Macrophage = c("CD68", "CD163", "MRC1"),
  Mast_cell = c("TPSAB1", "TPSB2", "CPA3"),
  MDSC = c("CD33", "ITGAM", "HLA-DRA"),
  Memory_B_cell = c("CD27", "CD19", "MS4A1"),
  Natural_killer_cell = c("NCAM1", "NCR1", "NCR3"),
  Natural_killer_T_cell = c("NCAM1", "CD3D", "CD3E"),
  Neutrophil = c("CEACAM8", "FCGR3B", "FPR1", "SIGLEC5"),
  Plasmacytoid_dendritic_cell = c("CLEC4C", "LILRA4", "PACSIN1"),
  Regulatory_T_cell = c("FOXP3", "IL2RA", "IKZF2"),
  T_follicular_helper_cell = c("CXCL13", "BCL6", "ICOS"),
  Type_1_T_helper_cell = c("TBX21", "IFNG", "CXCR3"),
  Type_17_T_helper_cell = c("RORC", "IL17A", "CCR6"),
  Type_2_T_helper_cell = c("GATA3", "IL4", "CCR4"),
  Monocyte = c("CD14", "CD163", "CD68")
)

# ==============================================================================
# 3. HELPER FUNCTIONS 
# ==============================================================================
# Standardize messy GEO labels into strictly "Tumor" or "Normal"
standardize_labels <- function(labels) {
  new_labels <- rep(NA, length(labels))
  labels_lower <- tolower(as.character(labels))
  
  normal_patterns <- c("normal", "control", "adjacent", "healthy", "unremarkable")
  tumor_patterns  <- c("tumor", "tumour", "cancer", "carcinoma", "idc", "tum", "adenocarcinoma")
  
  for(i in seq_along(labels_lower)) {
    lbl <- labels_lower[i]
    if (any(sapply(normal_patterns, function(p) grepl(p, lbl)))) {
      new_labels[i] <- "Normal"
    } else if (any(sapply(tumor_patterns, function(p) grepl(p, lbl)))) {
      new_labels[i] <- "Tumor"
    }
  }
  return(new_labels)
}

# The main analysis function (Now runs on a single dataset)
run_analysis <- function(expr, meta, tissue_name, ds_name, output_dir) {
  print(paste("   Running GSVA Analysis for:", ds_name))
  
  # Prefix for all output files (e.g., LUNG_GSE10072_)
  prefix <- paste(tissue_name, ds_name, sep="_")
  
  # 1. Run GSVA
  params <- ssgseaParam(expr, cell_markers, normalize = TRUE)
  ssgsea_res <- gsva(params, verbose=FALSE)
  
  # 2. Normalize scores (0-1) safely
  ssgsea_norm <- apply(ssgsea_res, 1, function(x) { 
    if(max(x) == min(x)) return(rep(0, length(x)))
    (x - min(x)) / (max(x) - min(x)) 
  })
  
  # Clean Data for Boxplot 
  plot_df <- as.data.frame(ssgsea_norm)
  plot_df$Group <- meta[rownames(plot_df), "Group"]
  
  # Remove missing groups
  plot_df <- plot_df[!is.na(plot_df$Group) & plot_df$Group != "", ]
  
  print(paste("   -> Total valid samples plotted:", nrow(plot_df)))
  
  # Format for ggplot
  plot_melt <- melt(plot_df, id.vars="Group")
  
  # --- PLOTTING ---
  
  # Figure A: Boxplot
  pdf(file.path(output_dir, paste0(prefix, "_Fig4A_Infiltration.pdf")), width=12, height=6)
  p1 <- ggplot(plot_melt, aes(x=variable, y=value, fill=Group)) +
    geom_boxplot(outlier.size=0.1) + theme_bw() +
    labs(title=paste(ds_name, "Immune Infiltration"), x="", y="ssGSEA Score") +
    theme(axis.text.x = element_text(angle=90, hjust=1)) +
    scale_fill_manual(values=c("Normal"="blue", "Tumor"="red")) 
  print(p1)
  dev.off()
  
  # Figure B: Interactions
  cor_cells <- cor(ssgsea_norm, method="spearman")
  pdf(file.path(output_dir, paste0(prefix, "_Fig4B_Network.pdf")), width=10, height=10)
  corrplot(cor_cells, method="circle", type="upper", tl.col="black", tl.cex=0.6, 
           title=paste(ds_name, "Immune Network"), mar=c(0,0,2,0))
  dev.off()
  
  # Figure C: Hub Genes Heatmap & P-Value Extraction
  valid_hubs <- intersect(rownames(expr), hub_genes)
  if(length(valid_hubs) > 0) {
    hub_expr <- expr[valid_hubs, ]
    
    cor_gene_cell <- cor(t(hub_expr), ssgsea_norm, method="spearman")
    
    pdf(file.path(output_dir, paste0(prefix, "_Fig4C_HubGenes.pdf")), width=8, height=5)
    pheatmap(cor_gene_cell, display_numbers=TRUE, fontsize_number=6,
             color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
             main=paste(ds_name, "Hub Genes vs Immune Cells"))
    dev.off()
    
    # --- EXTRACT EXACT P-VALUES ---
    combined_data <- cbind(t(hub_expr), ssgsea_norm)
    cor_results <- rcorr(as.matrix(combined_data), type="spearman")
    
    gene_names <- valid_hubs
    cell_names <- colnames(ssgsea_norm)
    
    final_r <- cor_results$r[gene_names, cell_names]
    final_p <- cor_results$P[gene_names, cell_names]
    
    write.csv(final_r, file.path(output_dir, paste0(prefix, "_Hub_vs_Immune_Corr_r.csv")))
    write.csv(final_p, file.path(output_dir, paste0(prefix, "_Hub_vs_Immune_Corr_Pvals.csv")))
  } else {
    print(paste("   -> WARNING: None of the hub genes were found in", ds_name))
  }
}

# Function to loop through tissue directories and process datasets ONE BY ONE
process_tissue_per_dataset <- function(tissue_path, tissue_name) {
  subfolders <- list.dirs(tissue_path, recursive = FALSE)
  subfolders <- subfolders[grep("GSE", subfolders)]
  
  for (folder in subfolders) {
    ds_name <- basename(folder)
    print(paste("--------------------------------------------------"))
    print(paste("Processing Dataset:", ds_name))
    
    rdata_path <- file.path(folder, "expression_data_step2_geneLevel.RData")
    scores_path <- file.path(folder, "TEX_Scores_With_Labels.csv")
    
    if (file.exists(rdata_path) && file.exists(scores_path)) {
      load(rdata_path) 
      if(exists("expr_matrix")) {
        current_expr <- as.matrix(expr_matrix)
        rm(expr_matrix) 
        
        current_scores <- read.csv(scores_path)
        current_meta <- unique(current_scores[, c("SampleID", "Group")])
        
        # APPLY THE LABEL STANDARDIZER
        current_meta$Group <- standardize_labels(current_meta$Group)
        
        # Remove any rows where group couldn't be determined
        current_meta <- current_meta[!is.na(current_meta$Group), ]
        rownames(current_meta) <- current_meta$SampleID
        
        common <- intersect(colnames(current_expr), rownames(current_meta))
        if (length(common) > 0) {
          expr <- current_expr[, common]
          meta <- current_meta[common, , drop=FALSE]
          
          # Run analysis and save outputs IN THIS SPECIFIC FOLDER
          run_analysis(expr, meta, tissue_name, ds_name, folder)
          
        } else {
          print(paste("   -> ERROR: No common samples found for", ds_name))
        }
      }
    } else {
      print(paste("   -> WARNING: Missing RData or CSV file in", ds_name))
    }
  }
}

# ==============================================================================
# 4. EXECUTION
# ==============================================================================
# --- PROCESS LUNG ---
print(">>> STARTING PER-DATASET LUNG ANALYSIS <<<")
process_tissue_per_dataset(lung_path, "LUNG")

# --- PROCESS BREAST ---
print(">>> STARTING PER-DATASET BREAST ANALYSIS <<<")
process_tissue_per_dataset(breast_path, "BREAST")
# Check individual GSE folders for the PDFs and CSV files."
