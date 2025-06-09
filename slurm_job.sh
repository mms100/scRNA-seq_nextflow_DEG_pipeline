#!/bin/bash
#SBATCH --job-name=nextflow_pbmc     # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --output=Path/to/logs/output.%J.%x.txt
#SBATCH --error=Path/to/logs/error.%J.%x.txt
### Time to execute
#SBATCH --time=96:00:00

### amount of memory 
#SBATCH --mem-per-cpu=100G

### amount of cores
#SBATCH --cpus-per-task=10

# Load necessary modules on HPC
module load scRNA/1.0.4  
module load R/4.2.3  


# Run pipeline with batch
//path/to/nextflow_excutor run path/to/main.nf  --object "/path/to/seurat_object.rds" --cond1 "g1" --cond2 "g2" --annotation "RNA_snn_res.0.8" --cond_colname "groups"  --batch_colname "letter.idents" --output_1 "Pbmc_batch"


# Run pipeline without batch
#//path/to/nextflow_excutor run path/to/main.nf  --object "/path/to/seurat_object.rds" --cond1 "g1" --cond2 "g2" --annotation "RNA_snn_res.0.8" --cond_colname "groups"  --output_1 "Pbmc_no_batch"