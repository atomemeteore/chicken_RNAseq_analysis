#!/bin/bash

# Exit on error and print commands
set -e
set -x

# System-specific settings
N_THREADS=6  # Using 75% of available threads for stability
MAX_MEM="2G"  # Memory usage for sorting
TMP_DIR="/tmp/chicken_shear_force_analysis"  # Temporary directory for processing

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
source activate chicken_shear_force_env

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
    local CONDITION=$2
    local R1="$SAMPLE_ID"_1.fastq.gz
    local R2="$SAMPLE_ID"_2.fastq.gz
    
    echo "Processing sample: $SAMPLE_ID ($CONDITION)"
    
    # Quality control
    cd ../data/raw
    mkdir -p ../../results/qc/raw/${CONDITION}
    fastqc -t $N_THREADS -o ../../results/qc/raw/${CONDITION} ${R1} ${R2}
    
    # Trim adapters and low quality bases
    cd ../processed
    trimmomatic PE -threads $N_THREADS \
        ../raw/${R1} ../raw/${R2} \
        ${CONDITION}_1_trimmed_paired.fastq.gz ${CONDITION}_1_trimmed_unpaired.fastq.gz \
        ${CONDITION}_2_trimmed_paired.fastq.gz ${CONDITION}_2_trimmed_unpaired.fastq.gz \
        ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    
    # QC on trimmed reads
    mkdir -p ../../results/qc/trimmed/${CONDITION}
    fastqc -t $N_THREADS -o ../../results/qc/trimmed/${CONDITION} \
        ${CONDITION}_1_trimmed_paired.fastq.gz ${CONDITION}_2_trimmed_paired.fastq.gz
    
    # Align reads
    echo "Aligning reads for ${CONDITION}..."
    hisat2 -p $N_THREADS --max-intronlen 50000 \
        -x ../../genome/chicken_index \
        -1 ${CONDITION}_1_trimmed_paired.fastq.gz \
        -2 ${CONDITION}_2_trimmed_paired.fastq.gz \
        -S ../../results/alignment/${CONDITION}.sam \
        2> ../../results/alignment/${CONDITION}_alignment_stats.txt
    
    # Convert SAM to BAM and sort
    cd ../../results/alignment
    samtools view -@ $N_THREADS -bS ${CONDITION}.sam > "$TMP_DIR/${CONDITION}.bam"
    samtools sort -@ $N_THREADS \
        -m $MAX_MEM \
        -T "$TMP_DIR/sort_tmp" \
        "$TMP_DIR/${CONDITION}.bam" \
        -o ${CONDITION}.sorted.bam
    
    # Index BAM
    samtools index ${CONDITION}.sorted.bam
    
    # Remove intermediate files
    rm -f ${CONDITION}.sam "$TMP_DIR/${CONDITION}.bam"
    
    # Count features
    echo "Counting features for ${CONDITION}..."
    cd ..
    featureCounts -T $N_THREADS -p \
        -a ../genome/Gallus_gallus.GalGal5.89.gtf \
        -o counts/${CONDITION}_counts.txt \
        alignment/${CONDITION}.sorted.bam
}

# Process all samples from samples.txt
echo "Processing all samples..."
while read -r line; do
    # Skip comments and empty lines
    [[ $line =~ ^#.*$ ]] && continue
    [[ -z "$line" ]] && continue
    
    # Extract accession and condition
    accession=$(echo $line | cut -f1)
    condition=$(echo $line | cut -f2)
    process_sample $accession $condition
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