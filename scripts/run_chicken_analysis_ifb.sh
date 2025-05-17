#!/bin/bash

#SBATCH --job-name=chicken_rnaseq
#SBATCH --output=chicken_rnaseq_%j.out
#SBATCH --error=chicken_rnaseq_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=fast

# Exit on error and print commands
set -e
set -x

# Base directories
BASE_DIR="$HOME/ondemand/data/sys/dashboard/batch_connect/sys/jupyter/core/chicken_rnaseq"
GENOME_DIR="${BASE_DIR}/genome"
DATA_DIR="${BASE_DIR}/data"
RESULTS_DIR="${BASE_DIR}/results"
TMP_DIR="/tmp/${USER}_chicken_rnaseq_${SLURM_JOB_ID}"

# Create directory structure
echo "Creating directory structure..."
mkdir -p ${DATA_DIR}/{raw,processed} \
    ${GENOME_DIR} \
    ${RESULTS_DIR}/{qc,alignment,counts,expression_atlas} \
    "$TMP_DIR"

# Cleanup function
cleanup() {
    echo "Cleaning up temporary files..."
    rm -rf "$TMP_DIR"
}
trap cleanup EXIT

# Load required modules
module purge
module load sra-tools
module load fastqc
module load trimmomatic
module load star
module load samtools
module load subread  # for featureCounts
module load multiqc

# Verify modules are loaded correctly by checking for their main executables
declare -A module_commands=(
    ["sra-tools"]="fastq-dump"
    ["fastqc"]="fastqc"
    ["trimmomatic"]="trimmomatic"
    ["star"]="STAR"
    ["samtools"]="samtools"
    ["subread"]="featureCounts"
    ["multiqc"]="multiqc"
)

echo "Verifying module availability..."
for module in "${!module_commands[@]}"; do
    cmd="${module_commands[$module]}"
    if ! command -v "$cmd" &> /dev/null; then
        echo "Error: Command '$cmd' from module '$module' is not available. Please check module loading."
        exit 1
    fi
done
echo "All required modules are available."

# Use conda R instead of module
export PATH="/shared/ifbstor1/software/miniconda/envs/r-4.2/bin:$PATH"

# URLs for Ensembl data (latest version)
GENOME_URL="https://ftp.ensembl.org/pub/release-110/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz"
GTF_URL="https://ftp.ensembl.org/pub/release-110/gtf/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.110.gtf.gz"

# Download and prepare reference if not already done
if [ ! -f "${GENOME_DIR}/reference_preparation_complete" ]; then
    echo "Downloading and preparing reference genome..."
    cd ${GENOME_DIR}
    
    # Download files
    wget -c ${GENOME_URL} -O chicken_genome.fa.gz
    wget -c ${GTF_URL} -O chicken_genes.gtf.gz
    
    # Decompress
    gunzip -f chicken_genome.fa.gz
    gunzip -f chicken_genes.gtf.gz
    
    # Index with STAR
    mkdir -p star_index
    STAR --runMode genomeGenerate \
         --genomeDir star_index \
         --genomeFastaFiles chicken_genome.fa \
         --sjdbGTFfile chicken_genes.gtf \
         --runThreadN $SLURM_CPUS_PER_TASK \
         --sjdbOverhang 99
    
    touch "${GENOME_DIR}/reference_preparation_complete"
fi

# Function to process a single RNA-seq sample
process_sample() {
    local SAMPLE_ID=$1
    local TISSUE_TYPE=$2
    
    echo "Processing sample: $SAMPLE_ID ($TISSUE_TYPE)"
    
    # Quality control on raw data
    echo "Running FastQC on raw data..."
    mkdir -p ${RESULTS_DIR}/qc/raw/${TISSUE_TYPE}
    fastqc -t $SLURM_CPUS_PER_TASK -o ${RESULTS_DIR}/qc/raw/${TISSUE_TYPE} \
        ${DATA_DIR}/raw/${SAMPLE_ID}_1_${TISSUE_TYPE}.fastq.gz
    
    # Trim reads
    echo "Trimming reads..."
    cd ${DATA_DIR}/processed
    trimmomatic SE -threads $SLURM_CPUS_PER_TASK \
        ${DATA_DIR}/raw/${SAMPLE_ID}_1_${TISSUE_TYPE}.fastq.gz \
        ${TISSUE_TYPE}_trimmed.fastq.gz \
        ILLUMINACLIP:/shared/apps/trimmomatic/0.39/adapters/TruSeq3-SE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    
    # QC on trimmed reads
    echo "Running FastQC on trimmed data..."
    mkdir -p ${RESULTS_DIR}/qc/trimmed/${TISSUE_TYPE}
    fastqc -t $SLURM_CPUS_PER_TASK -o ${RESULTS_DIR}/qc/trimmed/${TISSUE_TYPE} \
        ${TISSUE_TYPE}_trimmed.fastq.gz
    
    # Align reads with STAR
    echo "Aligning reads..."
    cd ${RESULTS_DIR}/alignment
    STAR --genomeDir ${GENOME_DIR}/star_index \
         --readFilesIn ${DATA_DIR}/processed/${TISSUE_TYPE}_trimmed.fastq.gz \
         --readFilesCommand zcat \
         --outFileNamePrefix ${TISSUE_TYPE}_ \
         --outSAMtype BAM SortedByCoordinate \
         --runThreadN $SLURM_CPUS_PER_TASK \
         --outBAMsortingThreadN $SLURM_CPUS_PER_TASK
    
    # Index BAM
    samtools index ${TISSUE_TYPE}_Aligned.sortedByCoord.out.bam
    
    # Count features
    echo "Counting features..."
    mkdir -p ${RESULTS_DIR}/counts
    featureCounts -T $SLURM_CPUS_PER_TASK \
        -a ${GENOME_DIR}/chicken_genes.gtf \
        -o ${RESULTS_DIR}/counts/${TISSUE_TYPE}_counts.txt \
        ${TISSUE_TYPE}_Aligned.sortedByCoord.out.bam
}

# Process available samples
echo "Processing available samples..."
cd ${DATA_DIR}/raw

# Check if we only need to create the expression matrix
ALL_COUNTS_EXIST=true
for sample_id in "${!samples[@]}"; do
    tissue_desc="${samples[$sample_id]}"
    if [ ! -f "${RESULTS_DIR}/counts/${tissue_desc}_counts.txt" ]; then
        ALL_COUNTS_EXIST=false
        break
    fi
done

if [ "$ALL_COUNTS_EXIST" = true ]; then
    echo "All count files exist. Skipping to expression matrix creation..."
else
    # Check if raw data files exist
    if [ ! "$(ls -A ${DATA_DIR}/raw/*.fastq.gz 2>/dev/null)" ]; then
        echo "Error: No .fastq.gz files found in ${DATA_DIR}/raw/"
        echo "Please transfer the raw data files first using scp from your local machine"
        exit 1
    fi

    # Define samples to process
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

    # Process each sample
    for sample_id in "${!samples[@]}"; do
        tissue_desc="${samples[$sample_id]}"
        echo "Processing sample: $sample_id ($tissue_desc)"
        process_sample "$sample_id" "$tissue_desc"
    done

    # Generate MultiQC report
    echo "Generating MultiQC report..."
    cd ${RESULTS_DIR}/qc
    multiqc .
fi

# Create expression matrix
echo "Creating expression matrix..."
cd ${RESULTS_DIR}/counts

# Verify we're in the correct directory
if [ ! -d "${RESULTS_DIR}/counts" ]; then
    echo "Error: Counts directory not found at ${RESULTS_DIR}/counts"
    exit 1
fi

# Get list of count files and verify count
count_files=($(ls *_counts.txt 2>/dev/null))
expected_count=${#samples[@]}  # Use the number of samples defined earlier
if [ ${#count_files[@]} -eq 0 ]; then
    echo "Error: No count files found in ${RESULTS_DIR}/counts"
    exit 1
elif [ ${#count_files[@]} -ne $expected_count ]; then
    echo "Warning: Found ${#count_files[@]} count files, expected $expected_count"
    echo "Continuing with available files..."
fi

# Extract gene IDs from first file
echo "Extracting gene IDs..."
awk 'NR>1 {print $1}' "${count_files[0]}" > genes.txt

# Create header line with sample names
echo -n "Gene_ID" > expression_matrix.txt
for f in "${count_files[@]}"; do
    sample_name=$(basename "$f" _counts.txt)
    echo -n -e "\t$sample_name" >> expression_matrix.txt
done
echo "" >> expression_matrix.txt

# Combine counts from all files
echo "Combining count data..."
awk 'NR>1 {printf "%s", $1; for(i=7; i<=NF; i++) printf "\t%s", $i; printf "\n"}' "${count_files[0]}" > temp_counts.txt

for f in "${count_files[@]:1}"; do
    awk 'NR>1 {for(i=7; i<=NF; i++) printf "\t%s", $i; printf "\n"}' "$f" > temp.txt
    paste temp_counts.txt temp.txt > temp_combined.txt
    mv temp_combined.txt temp_counts.txt
done

# Create final matrix
cat temp_counts.txt >> expression_matrix.txt

# Clean up temporary files
rm -f temp_counts.txt temp.txt genes.txt

# Verify the matrix was created and has expected format
if [ -f "expression_matrix.txt" ]; then
    echo "Expression matrix created successfully"
    
    # Check matrix format
    header_count=$(head -n1 expression_matrix.txt | awk '{print NF}')
    if [ $header_count -lt 2 ]; then
        echo "Error: Matrix header appears malformed (less than 2 columns)"
        exit 1
    fi
    
    # Report dimensions
    total_lines=$(wc -l < expression_matrix.txt)
    genes=$((total_lines - 1))  # Subtract header line
    samples=$((header_count - 1))  # Subtract Gene_ID column
    
    echo "Matrix dimensions:"
    echo "- $genes genes"
    echo "- $samples samples"
    
    # Verify first column header
    if ! head -n1 expression_matrix.txt | grep -q "^Gene_ID"; then
        echo "Warning: First column header is not 'Gene_ID' as expected"
    fi
    
    echo "Analysis complete! Results are in: ${RESULTS_DIR}"
    echo "Expression matrix is in: ${RESULTS_DIR}/counts/expression_matrix.txt"
else
    echo "Error: Failed to create expression matrix"
    exit 1
fi 