#!/bin/bash

# Exit on error
set -e

# Set up base directories
BASE_DIR="$(pwd)"
DATA_DIR="${BASE_DIR}/chicken_atlas/data/raw"
QC_DIR="${BASE_DIR}/chicken_atlas/results/qc/fastqc"
SCRIPT_DIR="${BASE_DIR}/chicken_atlas/scripts"

# Create output directory
mkdir -p "${QC_DIR}"

# Ensure conda and FastQC are available
if ! command -v conda &> /dev/null; then
    echo "conda not found. Please install conda first."
    exit 1
fi

# Create and activate conda environment with FastQC
echo "Setting up conda environment with FastQC..."
conda create -n qc_env -c bioconda fastqc multiqc -y || true
eval "$(conda shell.bash hook)"
conda activate qc_env

# Run FastQC on all samples
echo "Running FastQC on all samples..."
cd "${DATA_DIR}"

# Process all fastq.gz files
for file in *_1_*.fastq.gz; do
    echo "Processing $file..."
    fastqc -o "${QC_DIR}" -t 4 "$file"
done

# Run MultiQC to aggregate results
echo "Running MultiQC to create summary report..."
cd "${QC_DIR}"
multiqc .

echo "Quality control analysis complete!"
echo "FastQC reports are in: ${QC_DIR}"
echo "MultiQC report is in: ${QC_DIR}/multiqc_report.html" 