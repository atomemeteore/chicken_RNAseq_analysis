# IFB Cluster Guide for RNA-seq Analysis

## Overview
This guide provides instructions for running RNA-seq analyses on the IFB (Institut Français de Bioinformatique) cluster. It covers cluster connection, data transfer, analysis execution, and results retrieval.

## Prerequisites
- An IFB account (request at https://my.cluster.france-bioinformatique.fr/)
- Basic command line knowledge
- SSH client installed on your local machine
- Your IFB username and password

## Important Note About Directory Paths
The paths used in this guide are examples and will need to be adapted based on your IFB account organization. Common base directories on IFB include:
- `~/workspace/` - General workspace directory
- `~/work/` - Work directory with larger storage allocation
- `~/projects/` - Project-specific directory
- `~/ondemand/data/sys/dashboard/batch_connect/sys/jupyter/core/` - Jupyter-related workspace

Before running any scripts, make sure to:
1. Check your IFB account structure and available directories
2. Modify all paths in the scripts to match your account organization
3. Ensure you have sufficient storage space in your chosen directory

Example of path customization in scripts:
```bash
# Original path
WORK_DIR="~/workspace/project_name"

# Modify according to your IFB organization, e.g.:
WORK_DIR="~/work/project_name"  # If using work directory
# or
WORK_DIR="~/projects/project_name"  # If using projects directory
# or
WORK_DIR="~/ondemand/data/sys/dashboard/batch_connect/sys/jupyter/core/project_name"  # If using Jupyter workspace
```

## 1. Cluster Connection and Setup

1. Connect to the cluster:
```bash
ssh USERNAME@core.cluster.france-bioinformatique.fr
```
Replace `USERNAME` with your IFB username.

2. Start a tmux session (for disconnect protection):
```bash
# Create a new tmux session
tmux new -s rnaseq

# If disconnected, reattach to your session with:
tmux attach -t rnaseq
```

## 2. Data Transfer

1. On your local machine, organize and compress your data:
```bash
# Create a project directory
mkdir -p project_name/{data,genome,scripts,status}

# Organize your files into appropriate directories
# - Put raw data in data/
# - Put reference files in genome/
# - Put analysis scripts in scripts/

# Create an archive
tar czf project_data.tar.gz project_name/
```

2. Transfer files to the cluster:
```bash
# Transfer the archive
scp project_data.tar.gz USERNAME@core.cluster.france-bioinformatique.fr:~/workspace/

# For large raw data files, use compression:
scp -C data/*.fastq.gz USERNAME@core.cluster.france-bioinformatique.fr:~/workspace/project_name/data/
```

## 3. Cluster Setup

1. Set up your working directory:
```bash
# Go to your workspace
cd ~/workspace

# Create project directory structure
mkdir -p project_name/{data,genome,scripts,status,results}
cd project_name

# Extract your files (if using archive)
tar xzf ../project_data.tar.gz

# Make scripts executable
chmod +x scripts/*.sh
```

2. Verify your setup:
```bash
# Check directory structure
ls -R

# Check available space
df -h .
```

3. Update paths in your scripts:
```bash
# Set working directory
WORK_DIR="$PWD"
# Update paths in your scripts accordingly
```

## 4. Running Analysis

1. Submit genome preparation job (if needed):
```bash
# Submit with appropriate resources
sbatch scripts/prepare_genome.sh

# Monitor progress
watch -n 30 'squeue -u $USER'  # Press Ctrl-C to exit
tail -f prepare_genome_*.out   # Press Ctrl-C to exit
```

2. Submit analysis job:
```bash
# Check genome preparation completion (if applicable)
if [ -f "genome/preparation_complete" ]; then
    sbatch scripts/run_analysis.sh
else
    echo "Genome preparation not complete"
fi
```

## 5. Monitoring Jobs

1. Check job status:
```bash
# List your jobs
squeue -u $USER

# View job details
scontrol show job JOB_ID

# Check resource usage
sstat -j JOB_ID
```

2. Monitor logs:
```bash
# In a separate tmux pane (Ctrl-b then ")
tail -f job_name_*.out
tail -f job_name_*.err
```

## 6. Retrieving Results

From your local machine:
```bash
cd /path/to/local/project
scp -r USERNAME@core.cluster.france-bioinformatique.fr:~/workspace/project_name/results/ ./
```

## Important Notes

### Resource Management
- Request appropriate resources in your SLURM scripts:
```bash
#SBATCH --cpus-per-task=8     # Number of CPUs
#SBATCH --mem=32G             # Memory
#SBATCH --time=24:00:00       # Time limit
```

### Directory Structure
- Maintain a clear directory structure:
  - `genome/`: Reference files and indices
  - `data/`: Raw data
  - `scripts/`: Analysis scripts
  - `results/`: Output files
  - `status/`: Checkpoint files

### tmux Commands
- Essential tmux shortcuts:
  - `Ctrl-b d`: Detach session
  - `Ctrl-b "`: Split horizontally
  - `Ctrl-b %`: Split vertically
  - `Ctrl-b arrows`: Navigate between panes
  - `Ctrl-b c`: Create new window
  - `Ctrl-b n`: Next window
  - `Ctrl-b p`: Previous window

### Disconnection Handling
1. If disconnected:
   - Reconnect to the cluster
   - Reattach to tmux: `tmux attach -t rnaseq`
   - Your jobs continue running regardless of connection status

### Best Practices
- Use checkpoints in long-running jobs
- Monitor resource usage
- Keep logs organized
- Use appropriate job time limits
- Clean up unnecessary files regularly

## Useful Resources

- IFB Documentation: https://ifb-elixirfr.gitlab.io/cluster/doc/
- Cluster Status: https://www.france-bioinformatique.fr/cluster-status/
- Support: https://support.cluster.france-bioinformatique.fr/
- tmux Guide: https://tmuxcheatsheet.com/
- SLURM Documentation: https://slurm.schedmd.com/documentation.html

## Troubleshooting

### Common Issues
1. Connection Problems
   - Check your internet connection
   - Verify your IFB credentials
   - Ensure you're using the correct hostname

2. Job Failures
   - Check error logs
   - Verify resource requests
   - Ensure input files exist
   - Check disk space

3. Data Transfer Issues
   - Use `-C` flag for compression
   - Try smaller chunks
   - Check disk space on both ends

### Getting Help
- Contact IFB support
- Check cluster documentation
- Review job logs
- Use the IFB user forum 