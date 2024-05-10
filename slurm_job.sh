#!/bin/bash
#SBATCH --job-name=nextflow_pbmc     # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --output=next_flow_DEG/logs/output.%J.%x.txt
#SBATCH --error=next_flow_DEG/logs/error.%J.%x.txt
#SBATCH --mail-user=mmsaad100@gmail.com    # Where to send mail	
### Time to execute
#SBATCH --time=96:00:00

#SBATCH --mem-per-cpu=15G

### OpenMP threads
#SBATCH --cpus-per-task=3

# Load necessary modules (if any)
module load scRNA/1.0.4  # Nextflow requires Java
module load R/4.2.3  # If Nextflow is installed as a module

# Run Nextflow
//home/iu627335/nextflow run main.nf --object "/data/pbmc.Rds" --cond1 "g1" --cond2 "g2" --annotation "RNA_snn_res.1" --cond_colname "group"

