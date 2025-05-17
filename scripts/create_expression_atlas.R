# Script to create chicken RNA-seq expression atlas
# Based on Bush et al. 2018 (BMC Genomics)

# Load required libraries
suppressPackageStartupMessages({
    library(DESeq2)
    library(edgeR)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
    library(gridExtra)
})

# Function to read count data
read_counts <- function(file_path) {
    counts <- read.table(file_path, header=TRUE, row.names=1)
    # Remove length and other metadata columns
    counts <- counts[,6:ncol(counts)]
    return(counts)
}

# Function to normalize counts to TPM
counts_to_tpm <- function(counts, gene_lengths) {
    # Rate = reads / length
    rate <- counts / gene_lengths
    
    # TPM = rate / sum(rate) * 1e6
    tpm <- t(t(rate) / colSums(rate) * 1e6)
    return(tpm)
}

# Read all count files
count_files <- list.files(pattern="*_counts.txt")
sample_names <- gsub("_counts.txt", "", count_files)

# Create count matrix
counts_list <- lapply(count_files, read_counts)
counts_matrix <- do.call(cbind, counts_list)
colnames(counts_matrix) <- sample_names

# Get gene lengths from first count file
gene_info <- read.table(count_files[1], header=TRUE)
gene_lengths <- gene_info$Length

# Calculate TPM
tpm_matrix <- counts_to_tpm(counts_matrix, gene_lengths)

# Filter low expression genes (TPM > 1 in at least one sample)
expressed_genes <- rowSums(tpm_matrix > 1) > 0
filtered_tpm <- tpm_matrix[expressed_genes,]

# Calculate correlation matrix
cor_matrix <- cor(filtered_tpm)

# Generate heatmap
pdf("../expression_atlas/sample_correlation_heatmap.pdf")
pheatmap(cor_matrix,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Sample Correlation Heatmap")
dev.off()

# Save normalized expression data
write.table(filtered_tpm,
            file = "../expression_atlas/chicken_expression_atlas_TPM.txt",
            sep = "\t",
            quote = FALSE)

# Calculate PCA
pca <- prcomp(t(log2(filtered_tpm + 1)))

# Plot PCA
pdf("../expression_atlas/pca_plot.pdf")
plot(pca$x[,1:2],
     main = "PCA of Chicken RNA-seq Atlas",
     xlab = paste0("PC1 (", round(summary(pca)$importance[2,1] * 100, 1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca)$importance[2,2] * 100, 1), "%)"))
text(pca$x[,1:2], labels = rownames(pca$x), pos = 3, cex = 0.8)
dev.off()

# Generate summary statistics
stats <- data.frame(
    Sample = sample_names,
    Total_Genes = colSums(counts_matrix > 0),
    Expressed_Genes = colSums(filtered_tpm > 1)
)

write.table(stats,
            file = "../expression_atlas/sample_statistics.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# Print summary
cat("Expression atlas generation complete!\n")
cat("Number of samples processed:", length(sample_names), "\n")
cat("Number of genes after filtering:", nrow(filtered_tpm), "\n")
cat("Results saved in ../expression_atlas/\n") 