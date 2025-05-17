#!/bin/bash

# Exit on error, but with cleanup
set -e
trap 'cleanup_on_error $?' ERR

# Function to clean up on error
cleanup_on_error() {
    local exit_code=$1
    echo "Error occurred with exit code: ${exit_code}"
    # Remove partial downloads
    rm -f "${GENOME_DIR}"/*.gz.tmp
    exit "${exit_code}"
}

# Function to check available disk space (in GB)
check_disk_space() {
    local dir=$1
    local required_space=$2
    local available_space=$(df -BG "${dir}" | awk 'NR==2 {print $4}' | sed 's/G//')
    if [ "${available_space}" -lt "${required_space}" ]; then
        echo "Error: Not enough disk space. Required: ${required_space}GB, Available: ${available_space}GB"
        exit 1
    fi
}

# Function to check available memory (in GB)
check_memory() {
    local required_mem=$1
    local available_mem=$(free -g | awk 'NR==2 {print $7}')
    if [ "${available_mem}" -lt "${required_mem}" ]; then
        echo "Error: Not enough memory. Required: ${required_mem}GB, Available: ${available_mem}GB"
        exit 1
    fi
}

# Function to verify downloaded file
verify_download() {
    local file=$1
    if [ ! -f "${file}" ]; then
        echo "Error: Failed to download ${file}"
        return 1
    fi
    # Check if file is not empty
    if [ ! -s "${file}" ]; then
        echo "Error: Downloaded file ${file} is empty"
        rm -f "${file}"
        return 1
    fi
    return 0
}

# Set up base directories
BASE_DIR="$(pwd)"
GENOME_DIR="${BASE_DIR}/chicken_shear_force/genome"
mkdir -p "${GENOME_DIR}"

# URLs for Ensembl data (GCA_000002315.5)
GENOME_URL="https://ftp.ensembl.org/pub/release-110/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz"
GTF_URL="https://ftp.ensembl.org/pub/release-110/gtf/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.110.gtf.gz"

# Check for required disk space (4GB for genome, 2GB for GTF, 10GB for index)
echo "Checking available disk space..."
check_disk_space "${GENOME_DIR}" 16

# Check for required memory (minimum 16GB for STAR indexing)
echo "Checking available memory..."
check_memory 16

# Download and prepare genome if not already present
if [ ! -f "${GENOME_DIR}/chicken_genome.fa" ]; then
    echo "Downloading reference genome..."
    wget -c ${GENOME_URL} -O "${GENOME_DIR}/chicken_genome.fa.gz.tmp"
    verify_download "${GENOME_DIR}/chicken_genome.fa.gz.tmp"
    mv "${GENOME_DIR}/chicken_genome.fa.gz.tmp" "${GENOME_DIR}/chicken_genome.fa.gz"
    gunzip -f "${GENOME_DIR}/chicken_genome.fa.gz"
else
    echo "Genome file already exists, skipping download..."
fi

# Download and prepare GTF if not already present
if [ ! -f "${GENOME_DIR}/chicken_genes.gtf" ]; then
    echo "Downloading gene annotations..."
    wget -c ${GTF_URL} -O "${GENOME_DIR}/chicken_genes.gtf.gz.tmp"
    verify_download "${GENOME_DIR}/chicken_genes.gtf.gz.tmp"
    mv "${GENOME_DIR}/chicken_genes.gtf.gz.tmp" "${GENOME_DIR}/chicken_genes.gtf.gz"
    gunzip -f "${GENOME_DIR}/chicken_genes.gtf.gz"
else
    echo "GTF file already exists, skipping download..."
fi

# Create conda environment for STAR if it doesn't exist
if ! conda env list | grep -q "^align_env "; then
    echo "Creating conda environment for STAR..."
    conda create -n align_env -c bioconda star -y
fi

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate align_env

# Index the genome with STAR if index doesn't exist
if [ ! -f "${GENOME_DIR}/star_index/SA" ]; then
    echo "Indexing genome with STAR..."
    mkdir -p "${GENOME_DIR}/star_index"
    
    # Run STAR indexing with appropriate memory management
    STAR --runMode genomeGenerate \
         --genomeDir "${GENOME_DIR}/star_index" \
         --genomeFastaFiles "${GENOME_DIR}/chicken_genome.fa" \
         --sjdbGTFfile "${GENOME_DIR}/chicken_genes.gtf" \
         --runThreadN $(nproc --all) \
         --sjdbOverhang 99 \
         --limitGenomeGenerateRAM 40000000000
else
    echo "STAR index already exists, skipping indexing..."
fi

echo "Reference genome preparation complete!"
echo "Genome files are in: ${GENOME_DIR}"
echo "STAR index is in: ${GENOME_DIR}/star_index" 