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
   - Raw data QC
   - Adapter trimming
   - Quality filtering

2. Alignment and Quantification
   - Reference genome alignment
   - Expression quantification
   - Generation of expression matrix

3. Differential Expression Analysis
   - Comparison between high and low shear force groups
   - Statistical analysis
   - Visualization

4. Visualization
   - Correlation plots
   - MA plots
   - PCA plots
   - Shear force-specific expression patterns

## Requirements
See requirements.txt for Python package dependencies. Main requirements include:
- Python 3.6+
- NumPy
- Pandas
- Matplotlib
- Seaborn
- DESeq2

## Usage

### Local Analysis
1. Clone the repository:
```bash
git clone https://github.com/username/chicken_shear_force.git
cd chicken_shear_force
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Run the analysis pipeline:
```bash
python scripts/plot_expression_matrix.py
```

### IFB Cluster Analysis
This analysis can be run on the IFB (Institut Français de Bioinformatique) cluster for better performance and resource management. For detailed instructions, see `IFB_GUIDE_UNIVERSAL.md`.

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
sbatch scripts/prepare_genome.sh
```

2. Launch the analysis:
```bash
# After genome preparation is complete
if [ -f "genome/preparation_complete" ]; then
    sbatch scripts/run_analysis.sh
else
    echo "Genome preparation not complete"
fi
```

#### 5. Monitoring Jobs
- Check job status: `squeue -u $USER`
- View job details: `scontrol show job JOB_ID`
- Monitor logs: `tail -f job_name_*.out`
- Check resource usage: `sstat -j JOB_ID`

#### 6. Retrieving Results
From your local machine:
```bash
cd /path/to/local/project
scp -r USERNAME@core.cluster.france-bioinformatique.fr:~/workspace/chicken_rnaseq/results/ ./
```

#### Important Notes
- Resource Management:
  ```bash
  #SBATCH --cpus-per-task=8     # Number of CPUs
  #SBATCH --mem=32G             # Memory
  #SBATCH --time=24:00:00       # Time limit
  ```
- Directory Structure:
  - `genome/`: Reference files and indices
  - `data/`: Raw data
  - `scripts/`: Analysis scripts
  - `results/`: Output files
  - `status/`: Checkpoint files
- Essential tmux shortcuts:
  - `Ctrl-b d`: Detach session
  - `Ctrl-b "`: Split horizontally
  - `Ctrl-b %`: Split vertically
  - `Ctrl-b arrows`: Navigate between panes
  - `Ctrl-b c`: Create new window
  - `Ctrl-b n`: Next window
  - `Ctrl-b p`: Previous window

#### Best Practices
- Use checkpoints in long-running jobs
- Monitor resource usage
- Keep logs organized
- Use appropriate job time limits
- Clean up unnecessary files regularly

For more detailed instructions and troubleshooting, please refer to `IFB_GUIDE_UNIVERSAL.md`.

## Results
The analysis reveals expression patterns between high and low shear force groups in chicken breast muscle. Key findings include:

### 1. Sample Correlation Analysis
The correlation heatmap (sample_correlation_heatmap.pdf) shows the pairwise correlation coefficients between samples from different shear force conditions. This visualization helps identify:
- Clear clustering of samples by shear force condition
- High correlation between replicates within each condition
- Distinct expression patterns between high and low shear force groups

### 2. Principal Component Analysis
The PCA plot (pca_plot.pdf) demonstrates the separation between high and low shear force groups in reduced dimensional space:
- Clear segregation between conditions along principal components
- Tight clustering of biological replicates
- Explained variance for each principal component

### 3. Differential Expression Analysis
MA plots (MA_plot.pdf) showing the relationship between mean expression and log fold change:
- Highlights differentially expressed genes between conditions
- Shows the distribution of up and down-regulated genes
- Indicates statistical significance thresholds

### 4. Technical Reproducibility
Correlation plots between biological replicates demonstrate the quality and reproducibility of the data:

#### High Shear Force Replicates
The high_shear_force_replicate_correlations.pdf shows:
- Pairwise comparisons between all high shear force replicates
- Strong correlation coefficients (R² > 0.95)
- Consistent expression patterns across replicates

#### Low Shear Force Replicates
The low_shear_force_replicate_correlations.pdf demonstrates:
- Pairwise comparisons between all low shear force replicates
- High technical reproducibility (R² > 0.95)
- Consistent expression patterns across biological replicates

Key findings from the analysis:
- Clear separation between high and low shear force groups in both PCA and correlation analysis
- High reproducibility between biological replicates (R² > 0.95)
- Distinct differential expression patterns between conditions
- Identification of shear force-associated gene signatures

## Contributing
Please feel free to submit issues and pull requests.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Citation
If you use this analysis pipeline, please cite both this repository and the original data source:

```
Original Data:
Piórkowska K, et al. (2016). Genome-wide RNA-Seq analysis of breast muscles of two broiler 
chicken groups differing in shear force. Animal Genetics, 47(1):68-80.

Analysis Pipeline:
Nguyen A (2024) - Chicken Muscle Shear Force RNA-seq Analysis: A reproduction study
```

## Contact
For questions about this reproduction analysis, please open a GitHub issue or contact Alexis NGUYEN (alexisnguyen97@yahoo.fr) 