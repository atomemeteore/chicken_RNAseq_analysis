#!/bin/bash

# Exit on error and print commands
set -e
set -x

# System-specific settings
N_THREADS=6  # Using 75% of available threads for stability
MAX_MEM="2G"  # Memory usage for sorting
TMP_DIR="/tmp/chicken_rnaseq_analysis"  # Temporary directory for processing

# Create directory structure
echo "Creating directory structure..."
mkdir -p {data/{raw,processed},genome,results/{qc,alignment,counts,expression_atlas}} "$TMP_DIR"

# Cleanup function
cleanup() {
    echo "Cleaning up temporary files..."
    rm -rf "$TMP_DIR"
}
trap cleanup EXIT

# Activate conda environment
source activate chicken_rnaseq_env

# Download chicken reference genome and annotation
echo "Downloading reference genome and annotation..."
cd genome
wget --no-check-certificate https://ftp.ensembl.org/pub/release-89/fasta/gallus_gallus/dna/Gallus_gallus.GalGal5.dna.toplevel.fa.gz
wget --no-check-certificate https://ftp.ensembl.org/pub/release-89/gtf/gallus_gallus/Gallus_gallus.GalGal5.89.gtf.gz

# Decompress reference files
gunzip -f Gallus_gallus.GalGal5.dna.toplevel.fa.gz
gunzip -f Gallus_gallus.GalGal5.89.gtf.gz

# Index reference genome with HISAT2
echo "Indexing reference genome..."
hisat2-build -p $N_THREADS Gallus_gallus.GalGal5.dna.toplevel.fa chicken_index

# Function to process a single RNA-seq sample
process_sample() {
    local SAMPLE_ID=$1
    local TISSUE_TYPE=$2
    local R1="$SAMPLE_ID"_1.fastq.gz
    local R2="$SAMPLE_ID"_2.fastq.gz
    
    echo "Processing sample: $SAMPLE_ID ($TISSUE_TYPE)"
    
    # Quality control
    cd ../data/raw
    mkdir -p ../../results/qc/raw/${TISSUE_TYPE}
    fastqc -t $N_THREADS -o ../../results/qc/raw/${TISSUE_TYPE} ${R1} ${R2}
    
    # Trim adapters and low quality bases
    cd ../processed
    trimmomatic PE -threads $N_THREADS \
        ../raw/${R1} ../raw/${R2} \
        ${TISSUE_TYPE}_1_trimmed_paired.fastq.gz ${TISSUE_TYPE}_1_trimmed_unpaired.fastq.gz \
        ${TISSUE_TYPE}_2_trimmed_paired.fastq.gz ${TISSUE_TYPE}_2_trimmed_unpaired.fastq.gz \
        ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    
    # QC on trimmed reads
    mkdir -p ../../results/qc/trimmed/${TISSUE_TYPE}
    fastqc -t $N_THREADS -o ../../results/qc/trimmed/${TISSUE_TYPE} \
        ${TISSUE_TYPE}_1_trimmed_paired.fastq.gz ${TISSUE_TYPE}_2_trimmed_paired.fastq.gz
    
    # Align reads
    echo "Aligning reads for ${TISSUE_TYPE}..."
    hisat2 -p $N_THREADS --max-intronlen 50000 \
        -x ../../genome/chicken_index \
        -1 ${TISSUE_TYPE}_1_trimmed_paired.fastq.gz \
        -2 ${TISSUE_TYPE}_2_trimmed_paired.fastq.gz \
        -S ../../results/alignment/${TISSUE_TYPE}.sam \
        2> ../../results/alignment/${TISSUE_TYPE}_alignment_stats.txt
    
    # Convert SAM to BAM and sort
    cd ../../results/alignment
    samtools view -@ $N_THREADS -bS ${TISSUE_TYPE}.sam > "$TMP_DIR/${TISSUE_TYPE}.bam"
    samtools sort -@ $N_THREADS \
        -m $MAX_MEM \
        -T "$TMP_DIR/sort_tmp" \
        "$TMP_DIR/${TISSUE_TYPE}.bam" \
        -o ${TISSUE_TYPE}.sorted.bam
    
    # Index BAM
    samtools index ${TISSUE_TYPE}.sorted.bam
    
    # Remove intermediate files
    rm -f ${TISSUE_TYPE}.sam "$TMP_DIR/${TISSUE_TYPE}.bam"
    
    # Count features
    echo "Counting features for ${TISSUE_TYPE}..."
    cd ..
    featureCounts -T $N_THREADS -p \
        -a ../genome/Gallus_gallus.GalGal5.89.gtf \
        -o counts/${TISSUE_TYPE}_counts.txt \
        alignment/${TISSUE_TYPE}.sorted.bam
}

# Process all samples from samples.txt
echo "Processing all samples..."
while read -r line; do
    # Skip comments and empty lines
    [[ $line =~ ^#.*$ ]] && continue
    [[ -z "$line" ]] && continue
    
    # Extract accession and tissue type
    accession=$(echo $line | cut -f1)
    tissue_type=$(echo $line | cut -f2)
    process_sample $accession $tissue_type
done < ../data/samples.txt

# Generate MultiQC report
echo "Generating MultiQC report..."
cd results/qc
multiqc .

# Create expression atlas
echo "Creating expression atlas..."
cd ../counts
Rscript ../../scripts/create_expression_atlas.R

echo "Analysis complete! Check the results in the 'results' directory." 