# Script to create chicken muscle shear force RNA-seq analysis
# Analysis of gene expression patterns in chicken breast muscles under different shear force conditions

# Load required libraries
suppressPackageStartupMessages({
    library(DESeq2)
    library(edgeR)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
    library(gridExtra)
    library(reshape2)
    library(ggrepel)
})

# Read expression matrix
expr_matrix <- read.table("../expression_matrix_clean.txt", header=TRUE, row.names=1, check.names=FALSE)

# Create condition metadata
sample_info <- data.frame(
    sample = colnames(expr_matrix),
    condition = ifelse(grepl("high_shear", colnames(expr_matrix)), "High Shear Force", "Low Shear Force"),
    replicate = gsub(".*_rep([0-9]+)$", "\\1", colnames(expr_matrix))
)

# Print sample information for verification
print("Sample Information:")
print(sample_info)

# Filter low expression genes (>1 count in at least 3 samples)
keep <- rowSums(expr_matrix > 1) >= 3
filtered_expr <- expr_matrix[keep,]

# Create DGEList object
dge <- DGEList(counts = filtered_expr)
dge$samples$group <- factor(sample_info$condition)

# Normalize
dge <- calcNormFactors(dge)
norm_expr <- cpm(dge, log=TRUE)

# Generate correlation heatmap
pdf("../plots/sample_correlation_heatmap.pdf", width=10, height=8)
pheatmap(cor(norm_expr),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Sample Correlation Heatmap - Shear Force Conditions",
         annotation_col = data.frame(Condition = sample_info$condition,
                                   row.names = sample_info$sample),
         annotation_colors = list(Condition = c("High Shear Force" = "#E41A1C",
                                              "Low Shear Force" = "#377EB8")))
dev.off()

# PCA plot with ggplot2
pca_res <- prcomp(t(norm_expr))
pca_df <- data.frame(PC1 = pca_res$x[,1],
                     PC2 = pca_res$x[,2],
                     sample = rownames(pca_res$x),
                     condition = sample_info$condition)

var_explained <- summary(pca_res)$importance[2,1:2] * 100

pdf("../plots/pca_plot.pdf", width=10, height=8)
ggplot(pca_df, aes(x=PC1, y=PC2, color=condition, label=sample)) +
    geom_point(size=3) +
    geom_text_repel(size=3, show.legend=FALSE) +
    xlab(sprintf("PC1 (%.1f%%)", var_explained[1])) +
    ylab(sprintf("PC2 (%.1f%%)", var_explained[2])) +
    ggtitle("PCA of Chicken Muscle Shear Force RNA-seq") +
    theme_bw() +
    scale_color_manual(values=c("High Shear Force"="#E41A1C", "Low Shear Force"="#377EB8")) +
    theme(legend.title=element_blank())
dev.off()

# MA plots
design <- model.matrix(~sample_info$condition)
fit <- lmFit(norm_expr, design)
fit <- eBayes(fit)

pdf("../plots/MA_plot.pdf", width=10, height=8)
limma::plotMA(fit, main="MA Plot: High vs Low Shear Force")
abline(h=0, col="red", lty=2)
dev.off()

# Generate replicate correlation plots
plot_replicate_correlation <- function(condition) {
    condition_samples <- sample_info$sample[sample_info$condition == condition]
    pairs <- combn(condition_samples, 2, simplify=FALSE)
    
    plots <- lapply(pairs, function(pair) {
        df <- data.frame(
            rep1 = norm_expr[,pair[1]],
            rep2 = norm_expr[,pair[2]]
        )
        correlation <- cor(df$rep1, df$rep2)
        
        ggplot(df, aes(x=rep1, y=rep2)) +
            geom_point(alpha=0.3, size=1) +
            geom_abline(intercept=0, slope=1, color="red", linetype="dashed") +
            ggtitle(sprintf("%s vs %s\nR = %.3f", pair[1], pair[2], correlation)) +
            theme_bw() +
            theme(plot.title = element_text(size=10))
    })
    
    pdf(sprintf("../plots/%s_replicate_correlations.pdf", tolower(gsub(" ", "_", condition))),
        width=12, height=8)
    do.call(grid.arrange, c(plots, ncol=2))
    dev.off()
}

plot_replicate_correlation("High Shear Force")
plot_replicate_correlation("Low Shear Force")

# Save normalized expression data
write.table(norm_expr,
            file = "../expression_matrix_normalized.txt",
            sep = "\t",
            quote = FALSE)

# Generate summary statistics
stats <- data.frame(
    Sample = colnames(norm_expr),
    Condition = sample_info$condition,
    Total_Genes = colSums(expr_matrix > 0),
    Expressed_Genes = colSums(expr_matrix > 1),
    Filtered_Genes = colSums(filtered_expr > 0)
)

write.table(stats,
            file = "../sample_statistics.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# Print summary
cat("\nExpression analysis complete!\n")
cat("Number of samples processed:", ncol(norm_expr), "\n")
cat("Number of genes after filtering:", nrow(filtered_expr), "\n")
cat("Results saved in plots/ directory\n") 