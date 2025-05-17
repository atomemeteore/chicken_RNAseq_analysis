#!/usr/bin/env python3

# Set matplotlib backend before importing
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

def load_and_clean_data(matrix_file):
    """Load and clean the expression matrix data."""
    print(f"Loading expression matrix from {matrix_file}")
    
    # Read the file with explicit delimiter and handle whitespace
    df = pd.read_csv(matrix_file, sep='\s+')
    print(f"Initial data shape: {df.shape}")
    
    # Set the first column as index
    df = df.set_index('Gene_ID')
    print(f"Data shape after setting index: {df.shape}")
    
    # Convert all values to numeric
    df = df.apply(pd.to_numeric)
    print(f"Data shape after conversion: {df.shape}")
    
    # Basic statistics
    print("\nBasic statistics:")
    print(df.describe())
    
    return df

def plot_sample_correlations(df, output_dir):
    """Generate correlation scatter plots between replicates."""
    print("Generating sample correlation plots...")
    
    # Calculate correlation matrix
    corr = df.corr()
    
    # Create correlation heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(corr, annot=True, cmap='coolwarm', center=0, fmt='.3f')
    plt.title('Sample Correlation Heatmap')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/sample_correlation_heatmap.png", dpi=300)
    plt.close()
    
    # Create pairwise scatter plots for replicates
    tissues = ['liver', 'muscle', 'brain', 'kidney']
    for tissue in tissues:
        rep1 = f"{tissue}_rep1"
        rep2 = f"{tissue}_rep2"
        
        plt.figure(figsize=(8, 8))
        plt.scatter(np.log2(df[rep1] + 1), np.log2(df[rep2] + 1), alpha=0.5, s=1)
        plt.xlabel(f"{rep1} (log2(counts + 1))")
        plt.ylabel(f"{rep2} (log2(counts + 1))")
        plt.title(f"{tissue.capitalize()} Replicate Comparison")
        
        # Add correlation coefficient
        corr_val = df[rep1].corr(df[rep2])
        plt.text(0.05, 0.95, f"r = {corr_val:.3f}", 
                transform=plt.gca().transAxes)
        
        # Add density to show where most points are
        plt.hist2d(np.log2(df[rep1] + 1), np.log2(df[rep2] + 1), 
                  bins=50, cmap='viridis', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{tissue}_replicate_scatter.png", dpi=300)
        plt.close()

def plot_ma(df, output_dir):
    """Generate MA plots for each tissue comparison."""
    print("Generating MA plots...")
    tissues = ['liver', 'muscle', 'brain', 'kidney']
    
    for tissue in tissues:
        rep1 = f"{tissue}_rep1"
        rep2 = f"{tissue}_rep2"
        
        # Calculate M (log ratio) and A (mean average)
        M = np.log2((df[rep1] + 1) / (df[rep2] + 1))
        A = 0.5 * np.log2((df[rep1] + 1) * (df[rep2] + 1))
        
        plt.figure(figsize=(8, 8))
        
        # Add density to show where most points are
        plt.hist2d(A, M, bins=50, cmap='viridis', alpha=0.3)
        
        # Add scatter plot with smaller points
        plt.scatter(A, M, alpha=0.5, s=1, c='black')
        plt.axhline(y=0, color='r', linestyle='--')
        plt.xlabel('A (average log expression)')
        plt.ylabel('M (log ratio)')
        plt.title(f"{tissue.capitalize()} MA Plot")
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{tissue}_MA_plot.png", dpi=300)
        plt.close()

def plot_pca(df, output_dir):
    """Generate PCA plot of all samples."""
    print("Generating PCA plot...")
    
    # Prepare data for PCA
    X = np.log2(df + 1)  # Log transform
    X = StandardScaler().fit_transform(X.T)  # Scale the data
    
    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(X)
    
    # Calculate variance explained
    var_explained = pca.explained_variance_ratio_ * 100
    
    # Create DataFrame for plotting
    pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
    pca_df['sample'] = df.columns
    pca_df['tissue'] = [x.split('_')[0] for x in df.columns]
    
    # Plot
    plt.figure(figsize=(10, 8))
    sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='tissue', style='tissue', s=100)
    plt.title('PCA Plot of Samples')
    plt.xlabel(f'PC1 ({var_explained[0]:.1f}% variance explained)')
    plt.ylabel(f'PC2 ({var_explained[1]:.1f}% variance explained)')
    
    # Add sample labels
    for idx, row in pca_df.iterrows():
        plt.annotate(row['sample'], (row['PC1'], row['PC2']))
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/pca_plot.png", dpi=300)
    plt.close()

def main():
    # Define input/output paths
    matrix_file = "chicken_atlas/expression_matrix_clean.txt"
    output_dir = "chicken_atlas/plots"
    
    # Create output directory if it doesn't exist
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    # Load and clean data
    df = load_and_clean_data(matrix_file)
    
    # Generate plots
    plot_sample_correlations(df, output_dir)
    plot_ma(df, output_dir)
    plot_pca(df, output_dir)
    
    print(f"Plots have been saved to {output_dir}")

if __name__ == "__main__":
    main() 