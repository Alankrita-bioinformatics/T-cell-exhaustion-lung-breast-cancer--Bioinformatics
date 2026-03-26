# Step 1: Data Inspection and Preprocessing
setwd("your/path/here")

# Load libraries
library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)

options(timeout = 10000)

# Load dataset, change given geoid to your geoid.
gse <- getGEO("GSE19804", 
              GSEMatrix = TRUE,
              AnnotGPL = FALSE)

# Extract the ExpressionSet
gse_data <- gse[[1]]

# Sample metadata
metadata <- pData(gse_data)

# Expression matrix
expression_data <- exprs(gse_data)
summary(expression_data)
dim(expression_data)
# Check for missing values
cat("Missing values in expression data:", sum(is.na(expression_data)), "\n")
# Remove rows with ANY NA values
#expression_data <- na.omit(expression_data)

# Boxplot before and after normalization
png(file = file.path("boxplot_before_normalization.png"))
boxplot(expression_data, main = "Before Normalization", las = 2, outline = FALSE)
dev.off()

expression_data <- normalizeBetweenArrays(expression_data, method = "quantile")
expression_data <- log2(expression_data + 1)

png(file = file.path("boxplot_after_normalization.png"))
boxplot(expression_data, main = "After Normalization", las = 2, outline = FALSE)
dev.off()

# Filter low-expressed genes
# More lenient filter (standard practice)
keep <- rowSums(expression_data > 5) >= 10  # Top 10% samples
expression_data <- expression_data[keep, ]

cat("Dimensions after filtering:", dim(expression_data), "\n")

write.csv(expression_data, "expressionData.csv")
colnames(metadata)
table(metadata$characteristics_ch1)


# Extract 
unique(metadata$characteristics_ch1)

# 1. Check the actual text first (always do this!)
print(table(metadata$characteristics_ch1))

# 2. Use specific keywords found in GSE19804 or your geoid 
# e.g., "invasive" matches "invasive ductal carcinoma" and "normal" matches "normal lung tissue"
subtype_clean <- ifelse(grepl("invasive", metadata$characteristics_ch1, ignore.case=TRUE), "Tumor",
                        ifelse(grepl("normal", metadata$characteristics_ch1, ignore.case=TRUE), "Normal", "Unknown"))

# 3. Convert to factor and remove any "Unknown" if they exist
subtype <- factor(subtype_clean, levels = c("Normal", "Tumor"))

# 4. Verify the counts
print(table(subtype))

# PCA plot
pca_result <- prcomp(t(expression_data), scale. = TRUE)

pca_df <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2], Subtype = subtype)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Subtype)) +
  geom_point() +
  ggtitle("PCA Plot: Subtype Clustering") +
  theme_minimal()

ggsave(file.path("pca_plot.png"), plot = pca_plot)

# Save processed data
save(expression_data, metadata, subtype, file = file.path("expression_data_step1.RData"))

