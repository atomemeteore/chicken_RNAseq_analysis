#!/bin/bash

# Exit on error
set -e

echo "Installing RNA-seq analysis pipeline dependencies..."

# Create conda environment
conda create -n chicken_rnaseq_env python=3.9 -y

# Activate environment
source activate chicken_rnaseq_env

# Install bioinformatics tools
conda install -c bioconda -y \
    hisat2 \
    samtools \
    fastqc \
    trimmomatic \
    multiqc \
    subread

# Install Python dependencies
pip install -r ../requirements.txt

# Install R and required packages
conda install -c conda-forge r-base r-essentials -y

# Install Bioconductor packages
R -e 'if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install(c("DESeq2", "edgeR"))'

# Install additional R packages
R -e 'install.packages(c("ggplot2", "pheatmap", "RColorBrewer", "gridExtra"), repos="http://cran.us.r-project.org")'

echo "Installation complete! You can now run the analysis pipeline." 