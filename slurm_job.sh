#!/bin/bash
#SBATCH --job-name=nextflow_pbmc     # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --output=Path/to/scRNA-seq_nextflow_DEG_pipeline/log/output.%J.%x.txt
#SBATCH --error=Path/to/scRNA-seq_nextflow_DEG_pipeline/log/error.%J.%x.txt
### Time to execute
#SBATCH --time=96:00:00

### amount of memory 
#SBATCH --mem-per-cpu=100G

### amount of cores
#SBATCH --cpus-per-task=10

# Load necessary modules if you have on your HPC
module load R/4.3.2


# Run pipeline with batch
//path/to/nextflow_excutor run path/to/scRNA-seq_nextflow_DEG_pipeline/main.nf  --object "/path/to/seurat_object.rds" --cond1 "g1" --cond2 "g2" --annotation "RNA_snn_res.0.8" --cond_colname "stage"  --batch_colname "library_name" 


# Run pipeline without batch
#//path/to/nextflow_excutor run path/to/scRNA-seq_nextflow_DEG_pipeline/main.nf  --object "/path/to/seurat_object.rds" --cond1 "g1" --cond2 "g2" --annotation "RNA_snn_res.0.8" --cond_colname "stage"  
