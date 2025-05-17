#!/bin/bash

#SBATCH --job-name=download_ref
#SBATCH --output=download_ref_%j.out
#SBATCH --error=download_ref_%j.err
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=fast

# Exit on error, but with cleanup
set -e
trap 'cleanup_on_error $?' ERR

# Function to clean up on error
cleanup_on_error() {
    local exit_code=$1
    echo "Error occurred with exit code: ${exit_code}"
    if [ "${exit_code}" -eq 137 ]; then
        echo "ERROR: Job killed due to out of memory (OOM)"
        echo "Please resubmit with more memory"
    fi
    # Remove partial downloads and temp files
    rm -f "${GENOME_DIR}"/*.gz.tmp
    rm -f "${GENOME_DIR}"/Genome*
    # Save status for resuming
    echo "${exit_code}" > "${GENOME_DIR}/last_error.txt"
    exit "${exit_code}"
}

# Function to verify downloaded file
verify_download() {
    local file=$1
    if [ ! -f "${file}" ]; then
        echo "Error: Failed to download ${file}"
        return 1
    fi
    if [ ! -s "${file}" ]; then
        echo "Error: Downloaded file ${file} is empty"
        rm -f "${file}"
        return 1
    fi
    return 0
}

# Function to check file size
check_file_size() {
    local file=$1
    local size=$(stat -f%z "$file" 2>/dev/null || stat -c%s "$file")
    echo "File size of $file: $size bytes"
    if [ "$size" -lt 1000000 ]; then  # Moins de 1MB
        echo "Warning: File $file seems too small"
        return 1
    fi
    return 0
}

# Set up base directories
WORK_DIR="${HOME}/ondemand/data/sys/dashboard/batch_connect/sys/jupyter/core/chicken_rnaseq"
GENOME_DIR="${WORK_DIR}/genome"
mkdir -p "${GENOME_DIR}"

# Create status directory for checkpoints
STATUS_DIR="${GENOME_DIR}/status"
mkdir -p "${STATUS_DIR}"

# Load required modules
module purge
module load star

# URLs for Ensembl data
GENOME_URL="https://ftp.ensembl.org/pub/release-110/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz"
GTF_URL="https://ftp.ensembl.org/pub/release-110/gtf/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.110.gtf.gz"

# Download and prepare genome if not already done
if [ ! -f "${STATUS_DIR}/genome_done" ]; then
    echo "=== Téléchargement du génome de référence ==="
    cd "${GENOME_DIR}"
    if [ ! -f "chicken_genome.fa" ]; then
        echo "Téléchargement du génome..."
        wget -c ${GENOME_URL} -O "chicken_genome.fa.gz.tmp"
        verify_download "chicken_genome.fa.gz.tmp"
        mv "chicken_genome.fa.gz.tmp" "chicken_genome.fa.gz"
        gunzip -f "chicken_genome.fa.gz"
        check_file_size "chicken_genome.fa"
    fi
    touch "${STATUS_DIR}/genome_done"
else
    echo "Le génome est déjà téléchargé"
fi

# Download and prepare GTF if not already done
if [ ! -f "${STATUS_DIR}/gtf_done" ]; then
    echo "=== Téléchargement des annotations GTF ==="
    if [ ! -f "chicken_genes.gtf" ]; then
        echo "Téléchargement des annotations..."
        wget -c ${GTF_URL} -O "chicken_genes.gtf.gz.tmp"
        verify_download "chicken_genes.gtf.gz.tmp"
        mv "chicken_genes.gtf.gz.tmp" "chicken_genes.gtf.gz"
        gunzip -f "chicken_genes.gtf.gz"
        check_file_size "chicken_genes.gtf"
    fi
    touch "${STATUS_DIR}/gtf_done"
else
    echo "Les annotations sont déjà téléchargées"
fi

# Calculate available memory in GB
MEM_BYTES=$(free -b | awk '/^Mem:/{print $2}')
MEM_GB=$((MEM_BYTES/1024/1024/1024))
STAR_MEM=$((MEM_GB * 90 / 100))  # Use 90% of available memory

# Index the genome with STAR if not already done
if [ ! -f "${STATUS_DIR}/star_done" ]; then
    echo "=== Création de l'index STAR ==="
    if [ ! -f "star_index/SA" ]; then
        echo "Indexation du génome..."
        echo "Mémoire disponible: ${MEM_GB}GB, Allocation pour STAR: ${STAR_MEM}GB"
        mkdir -p "star_index"
        
        # Use more memory and threads for STAR
        STAR --runMode genomeGenerate \
             --genomeDir "star_index" \
             --genomeFastaFiles "chicken_genome.fa" \
             --sjdbGTFfile "chicken_genes.gtf" \
             --runThreadN $SLURM_CPUS_PER_TASK \
             --sjdbOverhang 99 \
             --limitGenomeGenerateRAM $((STAR_MEM * 1024**3))

        # Verify index creation
        if [ -f "star_index/SA" ]; then
            touch "${STATUS_DIR}/star_done"
            echo "Index STAR créé avec succès"
        else
            echo "Error: STAR index creation failed"
            exit 1
        fi
    fi
else
    echo "L'index STAR existe déjà"
fi

echo "=== Préparation du génome de référence terminée ! ==="
echo "Fichiers générés dans : ${GENOME_DIR}"
ls -lh "${GENOME_DIR}"
echo "Index STAR dans : ${GENOME_DIR}/star_index"

# Create completion flag
touch "${GENOME_DIR}/reference_preparation_complete" 