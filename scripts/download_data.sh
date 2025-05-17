#!/bin/bash

# Exit on error
set -e

# Data source information
# These RNA-seq data are from the study by Piórkowska et al. (2016) in Animal Genetics
# The SRA accessions were obtained from NCBI's Sequence Read Archive (SRA)
# BioProject: PRJNA297364
# Publication: https://doi.org/10.1111/age.12409
# Each accession corresponds to RNA-seq data from different chicken tissues:
# - SRR2554344, SRR2554345: Liver samples (2 replicates)
# - SRR2554362, SRR2554363: Muscle samples (2 replicates)
# - SRR2554364, SRR2554365: Brain samples (2 replicates)
# - SRR2554366, SRR2554367: Kidney samples (2 replicates)

# Set up base directories
BASE_DIR="$(pwd)"
DATA_DIR="${BASE_DIR}/chicken_atlas/data/raw"
SCRIPT_DIR="${BASE_DIR}/chicken_atlas/scripts"

# Create necessary directories
mkdir -p "${DATA_DIR}"

# Ensure conda and SRA toolkit are available
if ! command -v conda &> /dev/null; then
    echo "conda not found. Please install conda first."
    exit 1
fi

# Create and activate conda environment with SRA toolkit
echo "Setting up conda environment with SRA toolkit..."
conda create -n sra_env -c bioconda sra-tools -y || true
eval "$(conda shell.bash hook)"
conda activate sra_env

# Function to download SRA data directly from AWS S3
download_sra() {
    local ACCESSION=$1
    local SAMPLE_NAME=$2
    echo "Downloading $ACCESSION (Sample: $SAMPLE_NAME)..."
    
    cd "${DATA_DIR}"
    
    # Skip if the file already exists
    if [ -f "${ACCESSION}_1_${SAMPLE_NAME}.fastq.gz" ]; then
        echo "Sample $ACCESSION already downloaded, skipping..."
        return 0
    fi
    
    # Construct the S3 URL
    local S3_URL="https://sra-pub-run-odp.s3.amazonaws.com/sra/${ACCESSION}/${ACCESSION}"
    
    echo "Downloading from: $S3_URL"
    wget -c "$S3_URL" -O "${ACCESSION}.sra" || {
        echo "Failed to download $ACCESSION"
        return 1
    }
    
    # Convert SRA to FASTQ using fastq-dump
    echo "Converting SRA to FASTQ..."
    fastq-dump --split-files --gzip "${ACCESSION}.sra" || {
        echo "Failed to convert $ACCESSION to FASTQ"
        return 1
    }
    
    # Clean up SRA file
    rm "${ACCESSION}.sra"
    
    # Rename files to include sample name
    for f in ${ACCESSION}*.fastq.gz; do
        mv "$f" "${f/.fastq.gz/_${SAMPLE_NAME}.fastq.gz}"
    done
    
    cd "${SCRIPT_DIR}"
}

# Define samples and their names
declare -A SAMPLE_INFO=(
    ["SRR2554344"]="liver_rep1"
    ["SRR2554345"]="liver_rep2"
    ["SRR2554362"]="muscle_rep1"
    ["SRR2554363"]="muscle_rep2"
    ["SRR2554364"]="brain_rep1"
    ["SRR2554365"]="brain_rep2"
    ["SRR2554366"]="kidney_rep1"
    ["SRR2554367"]="kidney_rep2"
)

# Download all remaining datasets
echo "Starting downloads of remaining samples..."
FAILED_DOWNLOADS=()

for acc in "${!SAMPLE_INFO[@]}"; do
    if download_sra "$acc" "${SAMPLE_INFO[$acc]}"; then
        echo "Successfully processed $acc"
    else
        echo "Failed to process $acc"
        FAILED_DOWNLOADS+=("$acc")
    fi
done

# Report results
echo "Downloads complete! Data is in ${DATA_DIR}"
if [ ${#FAILED_DOWNLOADS[@]} -eq 0 ]; then
    echo "All samples were downloaded successfully!"
else
    echo "The following samples failed to download:"
    printf '%s\n' "${FAILED_DOWNLOADS[@]}"
fi

# Create sample information file
echo "Creating sample information file..."
cat > "${SCRIPT_DIR}/sample_info.txt" << EOL
# Sample information
declare -A samples=(
    ["SRR2554344"]="liver_rep1"
    ["SRR2554345"]="liver_rep2"
    ["SRR2554362"]="muscle_rep1"
    ["SRR2554363"]="muscle_rep2"
    ["SRR2554364"]="brain_rep1"
    ["SRR2554365"]="brain_rep2"
    ["SRR2554366"]="kidney_rep1"
    ["SRR2554367"]="kidney_rep2"
)
EOL

echo "Setup complete! You can now run the analysis pipeline." 