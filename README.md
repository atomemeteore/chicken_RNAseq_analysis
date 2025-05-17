# Chicken Muscle Shear Force RNA-seq Analysis

## Overview
This repository contains an RNA-seq analysis pipeline that examines gene expression patterns in chicken breast muscles under different shear force conditions. The data appears to be from Piórkowska, et al. (2016) in Animal Genetics, comparing transcriptional profiles between high and low shear force groups, though this needs to be confirmed as the original article is not publicly accessible.

## Data Source
- **Original Study**: The data is believed to be from Piórkowska K, et al. (2016). Genome-wide RNA-Seq analysis of breast muscles of two broiler chicken groups differing in shear force. Animal Genetics, 47(1):68-80. However, this needs to be confirmed as the article is not publicly available.
- **BioProject**: PRJNA297364
- **Sample Groups**: 
  - High shear force (4 replicates):
    - SRR2554364 (replicate 1)
    - SRR2554365 (replicate 2)
    - SRR2554366 (replicate 3)
    - SRR2554367 (replicate 4)
  - Low shear force (4 replicates):
    - SRR2554362 (replicate 1)
    - SRR2554363 (replicate 2)
    - SRR2554345 (replicate 3)
    - SRR2554344 (replicate 4)
- **Experimental Design**: Four biological replicates per condition
- **Sequencing Platform**: Illumina HiScanSQ
- **Library Type**: RNA-seq, single-end reads
- **Data Processing**: RNA-seq analysis performed on IFB cluster
- **Important Note**: While the BioProject and SRA accession numbers are confirmed, the association with the cited publication needs verification due to limited access to the original article. Users are encouraged to verify the data source independently.

## Directory Structure
```
.
├── data/               # Raw sequencing data and metadata
├── genome/            # Reference genome and annotation files
├── scripts/           # Analysis and plotting scripts
├── results/           # Analysis output files
├── plots/             # Generated figures and visualizations
└── requirements.txt   # Python dependencies
```

## Analysis Pipeline
1. Quality Control and Preprocessing
   - Raw data QC using FastQC (v0.11.9)
     * Quality score distribution analysis
     * Sequence duplication assessment
     * Adapter contamination detection
     * Base composition analysis
   - Adapter trimming with Trimmomatic (v0.39)
     * Removal of Illumina adapters
     * Quality filtering (leading/trailing bases)
     * Sliding window quality trimming
     * Minimum length filtering
   - Post-trimming QC with FastQC and MultiQC
     * Verification of adapter removal
     * Quality metrics visualization
     * Sample comparison reports

2. Alignment and Quantification
   - Reference genome alignment using HISAT2 (v2.2.1)
     * Splice-aware alignment
     * Multi-threaded processing
     * Optimized for RNA-seq data
   - SAMtools (v1.15) for BAM processing
     * SAM to BAM conversion
     * BAM sorting and indexing
     * Alignment statistics
   - Expression quantification with featureCounts (v2.0.1)
     * Gene-level counting
     * Multi-mapping handling
     * Strand-specific counting support

3. Differential Expression Analysis
   - DESeq2 (v1.32.0) for statistical analysis
     * Normalization for sequencing depth
     * Dispersion estimation
     * Differential expression testing
     * Log2 fold change shrinkage
   - Statistical thresholds:
     * Adjusted p-value < 0.05
     * |Log2 fold change| > 1

4. Visualization (R v4.1.0)
   - Correlation plots using pheatmap
     * Sample correlation heatmaps
     * Hierarchical clustering
   - MA plots with ggplot2
     * Log fold change vs mean expression
     * Significance highlighting
   - PCA plots using ggplot2
     * Dimension reduction
     * Sample clustering visualization
   - Custom expression pattern plots
     * Gene-specific visualizations
     * Group comparisons

## Results
The analysis reveals expression patterns between high and low shear force groups in chicken breast muscle. Key findings include:

### 1. Sample Correlation Analysis
The correlation heatmap shows the pairwise correlation coefficients between samples from different shear force conditions:

![Sample Correlation Heatmap](images/sample_correlation_heatmap_v2.png)

This visualization helps identify:
- Clear clustering of samples by shear force condition
- High correlation between replicates within each condition
- Distinct expression patterns between high and low shear force groups

### 2. Principal Component Analysis
The PCA plot demonstrates the separation between high and low shear force groups in reduced dimensional space:

![PCA Plot](images/pca_plot_v2.png)

Key observations:
- Clear segregation between conditions along principal components
- Tight clustering of biological replicates
- Explained variance for each principal component

### 3. Differential Expression Analysis
MA plots showing the relationship between mean expression and log fold change:

![MA Plot](images/MA_plot_v2.png)

This plot reveals:
- Highlights differentially expressed genes between conditions
- Shows the distribution of up and down-regulated genes
- Indicates statistical significance thresholds

### 4. Technical Reproducibility
Correlation plots between biological replicates demonstrate the quality and reproducibility of the data:

#### High Shear Force Replicates
![High Shear Force Correlations](images/high_shear_force_replicate_correlations_v2.png)

The correlations show:
- Pairwise comparisons between all high shear force replicates
- Strong correlation coefficients (R² > 0.95)
- Consistent expression patterns across replicates

#### Low Shear Force Replicates
![Low Shear Force Correlations](images/low_shear_force_replicate_correlations_v2.png)

The correlations demonstrate:
- Pairwise comparisons between all low shear force replicates
- High technical reproducibility (R² > 0.95)
- Consistent expression patterns across biological replicates

### Key Findings Summary
- Clear separation between high and low shear force groups in both PCA and correlation analysis
- High reproducibility between biological replicates (R² > 0.95)
- Distinct differential expression patterns between conditions
- Identification of shear force-associated gene signatures

## Tools and Dependencies

### Core Analysis Tools

- **SRA Tools**
  - Purpose: Download and handle SRA format sequencing data
  - Features:
    * Direct download from NCBI SRA
    * FASTQ format conversion
    * Compression handling

- **FastQC** (v0.11.9)
  - Purpose: Quality control of raw sequencing data
  - Features:
    * Per base quality scores
    * Sequence duplication levels
    * Adapter content detection
    * Base composition analysis
  - Usage: Initial and post-trimming QC

- **Trimmomatic** (v0.39)
  - Purpose: Adapter and quality trimming
  - Features: 
    * Illumina adapter removal
    * Sliding window trimming
    * Minimum length filtering
  - Parameters:
    * ILLUMINACLIP: TruSeq3-SE.fa:2:30:10
    * LEADING: 3
    * TRAILING: 3
    * SLIDINGWINDOW: 4:15
    * MINLEN: 36

- **STAR** (v2.7+)
  - Purpose: RNA-seq read alignment
  - Features:
    * Splice-aware alignment
    * Multi-threaded processing
    * Memory-efficient operation
  - Parameters:
    * Genome indexing with sjdbOverhang 99
    * BAM output sorted by coordinate
    * Multi-threading support

- **SAMtools** (v1.15)
  - Purpose: SAM/BAM manipulation
  - Features:
    * Format conversion
    * BAM sorting and indexing
    * Alignment statistics
    * BAM file compression

- **Subread package/featureCounts** (v2.0.1)
  - Purpose: Read counting and quantification
  - Features:
    * Gene-level counting
    * Multi-threading support
    * GTF/GFF format support
    * Strand-specific counting
  - Parameters:
    * Multi-threaded processing
    * Gene-level summarization

- **MultiQC**
  - Purpose: Aggregate QC reports
  - Features:
    * Combines reports from multiple tools
    * Interactive visualizations
    * Sample comparison
    * Quality metric summaries

### Statistical Analysis and Visualization
- **R** (v4.2)
  - Purpose: Statistical analysis and visualization
  - Environment: Local installation (not on IFB cluster)
  - Usage: Analysis performed locally using `create_expression_atlas.R` after transferring count data from IFB cluster
  - Key packages:
    * **DESeq2**: Differential expression analysis
    * **edgeR**: RNA-seq analysis
    * **ggplot2**: Data visualization
    * **pheatmap**: Heatmap generation
    * **RColorBrewer**: Color palettes
    * **gridExtra**: Plot layouts

Note: The RNA-seq data processing (quality control, alignment, and counting) is performed on the IFB cluster, but the final statistical analysis and visualization using R is done locally. This requires transferring the count data from the IFB cluster to your local machine using:
```bash
# From your local machine
scp USERNAME@core.cluster.france-bioinformatique.fr:~/workspace/project_name/results/counts/expression_matrix.txt ./
```

### Resource Management
- **SLURM** Workload Manager
  - Purpose: Job scheduling and resource allocation
  - Features:
    * Job queuing and management
    * Resource monitoring
    * Multi-node support
  - Parameters:
    * CPUs per task: 8
    * Memory: 64GB
    * Runtime: 24 hours
    * Partition: fast

### Version Control and Documentation
- **Git** (v2.34+)
  - Purpose: Version control and collaboration
  - Features:
    * Change tracking
    * Branch management
    * Collaborative development

All tools are configured with parameters optimized for chicken RNA-seq analysis on the IFB cluster. For specific version requirements and compatibility, see `requirements.txt`.

## Requirements
See requirements.txt for Python package dependencies. Main requirements include:
- Python 3.6+
- NumPy
- Pandas
- Matplotlib
- Seaborn
- DESeq2

## Usage

The analysis is performed in two parts:

1. Main RNA-seq Processing (On IFB Cluster)
   - Quality control
   - Read alignment
   - Feature counting
   - See `IFB_GUIDE.md` or `IFB_GUIDE_UNIVERSAL.md` for detailed instructions

2. Statistical Analysis (Local)
   - Transfer count data from IFB:
   ```bash
   # From your local machine
   scp USERNAME@core.cluster.france-bioinformatique.fr:~/workspace/project_name/results/counts/expression_matrix.txt ./
   ```
   - Run R analysis:
   ```bash
   Rscript scripts/create_expression_atlas.R
   ```

Note: The main data processing is performed on the IFB cluster due to computational requirements, while the final statistical analysis and visualization are done locally using R.

### IFB Cluster Analysis
This analysis can be run on the IFB (Institut Français de Bioinformatique) cluster for better performance and resource management. Two versions of the IFB guide are available:

1. `IFB_GUIDE.md` - A personalized guide with specific paths and settings for the original setup
2. `IFB_GUIDE_UNIVERSAL.md` - A universal guide that can be adapted for any user

The analysis uses the following IFB-specific scripts:
- `scripts/download_reference_ifb.sh` - Downloads and prepares the reference genome on the IFB cluster, with SLURM job configuration
- `scripts/run_chicken_analysis_ifb.sh` - Runs the complete RNA-seq analysis pipeline on the IFB cluster, with SLURM job configuration

These scripts have been specifically configured for the IFB cluster environment, including:
- SLURM job scheduling parameters
- IFB-specific module loading
- Proper resource allocation (CPU, memory, time)
- IFB-specific paths and environment setup
- Checkpoint system for long-running jobs

#### Prerequisites
- An IFB account (request at https://my.cluster.france-bioinformatique.fr/)
- Basic command line knowledge
- SSH client installed on your local machine
- Your IFB username and password

#### 1. Cluster Connection and Setup
1. Connect to the cluster:
```bash
ssh USERNAME@core.cluster.france-bioinformatique.fr
```

2. Start a tmux session (for disconnect protection):
```bash
tmux new -s rnaseq
# To reattach if disconnected: tmux attach -t rnaseq
```

#### 2. Data Transfer
1. On your local machine, organize and compress your data:
```bash
# Create a project directory and organize files
mkdir -p chicken_rnaseq/{data,genome,scripts,status}

# Create an archive
tar czf chicken_rnaseq.tar.gz chicken_rnaseq/
```

2. Transfer files to the cluster:
```bash
# Transfer the archive
scp chicken_rnaseq.tar.gz USERNAME@core.cluster.france-bioinformatique.fr:~/workspace/

# For large raw data files, use compression:
scp -C data/*.fastq.gz USERNAME@core.cluster.france-bioinformatique.fr:~/workspace/chicken_rnaseq/data/
```

#### 3. Cluster Setup
1. Set up your working directory:
```bash
# Go to your workspace
cd ~/workspace
mkdir -p chicken_rnaseq/{data,genome,scripts,status,results}
cd chicken_rnaseq

# Extract your files
tar xzf ../chicken_rnaseq.tar.gz

# Make scripts executable
chmod +x scripts/*.sh
```

#### 4. Running Analysis
1. Download and index reference genome:
```bash
# Submit the genome preparation job
sbatch scripts/download_reference_ifb.sh

# Monitor the job
watch -n 30 'squeue -u $USER'
tail -f download_ref_*.out
```

2. Launch the RNA-seq analysis:
```bash
# Check if genome preparation is complete
if [ -f "genome/reference_preparation_complete" ]; then
    # Submit the analysis job
    sbatch scripts/run_chicken_analysis_ifb.sh
    
    # Monitor the job
    watch -n 30 'squeue -u $USER'
    tail -f chicken_rnaseq_*.out
else
    echo "Genome preparation not complete"
fi
```

3. Monitor progress in the results directory:
```bash
# Check the progress of different steps
ls -l results/{qc,alignment,counts}/

# View MultiQC report when available
firefox results/qc/multiqc_report.html  # or your preferred browser
```

## Citation
If you use this analysis pipeline, please cite both this repository and the original data source:

```
Original Data:
Piórkowska K, et al. (2016). Genome-wide RNA-Seq analysis of breast muscles of two broiler 
chicken groups differing in shear force. Animal Genetics, 47(1):68-80.

Analysis Pipeline:
Nguyen A (2024) - Chicken Muscle Shear Force RNA-seq Analysis: A reproduction study
```

## Acknowledgments
This project was developed with the assistance of:
- [Cursor](https://cursor.sh/) - The AI-powered code editor that provided intelligent coding assistance throughout the project
- Claude (Anthropic) - The AI model powering Cursor that helped with code development, documentation, and analysis pipeline design

## Contact
For questions about this reproduction analysis, please open a GitHub issue or contact Alexis NGUYEN (alexisnguyen97@yahoo.fr)