# Load required libraries
library(dplyr)
library(tidyr)
library(tibble)

# Read the expression matrix
expr_data <- read.delim("expression_matrix_clean_fixed.txt", 
                       header=TRUE,
                       sep="\t",
                       stringsAsFactors=FALSE,
                       check.names=FALSE)

# Set Gene_ID as row names
rownames(expr_data) <- expr_data$Gene_ID
expr_data$Gene_ID <- NULL

# Calculate mean expression for each condition
high_shear_cols <- grep("high_shear", colnames(expr_data))
low_shear_cols <- grep("low_shear", colnames(expr_data))

high_shear_mean <- rowMeans(as.matrix(expr_data[, high_shear_cols]))
low_shear_mean <- rowMeans(as.matrix(expr_data[, low_shear_cols]))

# Calculate log2 fold change and mean expression
epsilon <- 1  # Add small value to avoid log(0)
log2_fc <- log2((high_shear_mean + epsilon) / (low_shear_mean + epsilon))
mean_expr <- (high_shear_mean + low_shear_mean) / 2

# Create results data frame
results <- data.frame(
  Gene_ID = rownames(expr_data),
  High_Shear_Mean = high_shear_mean,
  Low_Shear_Mean = low_shear_mean,
  Log2_FC = log2_fc,
  Mean_Expression = mean_expr,
  stringsAsFactors = FALSE
)

# Filter for expressed genes (mean expression > 10)
results_filtered <- results[results$Mean_Expression > 10, ]

# Sort by absolute fold change
results_filtered$Abs_Log2_FC <- abs(results_filtered$Log2_FC)
results_sorted <- results_filtered[order(-results_filtered$Abs_Log2_FC), ]

# Get top 10 genes for each condition
top_high_shear <- results_sorted[results_sorted$Log2_FC > 0, ][1:10, ]
top_low_shear <- results_sorted[results_sorted$Log2_FC < 0, ][1:10, ]

# Print results
cat("\nTop 10 genes with higher expression in high shear force condition:\n")
print(top_high_shear[, c("Gene_ID", "High_Shear_Mean", "Low_Shear_Mean", "Log2_FC")])

cat("\nTop 10 genes with higher expression in low shear force condition:\n")
print(top_low_shear[, c("Gene_ID", "High_Shear_Mean", "Low_Shear_Mean", "Log2_FC")])

# Save results to files
write.table(top_high_shear, 
            file = "results/expression_atlas/top_high_shear_genes.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

write.table(top_low_shear, 
            file = "results/expression_atlas/top_low_shear_genes.txt", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# Print summary statistics
cat("\nSummary Statistics:\n")
cat("Total genes analyzed:", nrow(results), "\n")
cat("Number of expressed genes (mean expression > 10):", nrow(results_filtered), "\n")
cat("Number of genes with higher expression in high shear force:", sum(log2_fc > 1), "\n")
cat("Number of genes with higher expression in low shear force:", sum(log2_fc < -1), "\n") 