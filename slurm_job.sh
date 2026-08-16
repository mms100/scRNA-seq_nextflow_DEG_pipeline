#!/bin/bash
#SBATCH --job-name=DEG_pipeline     # Job name
#SBATCH --mail-user=email    # Where to send mail
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --output=/scRNA-seq_nextflow_DEG_pipeline/log/output.%J.%x.txt
### Time to execute
#SBATCH --time=96:00:00

### amount of memory
#SBATCH --mem-per-cpu=8G

### amount of cores
#SBATCH --cpus-per-task=2

# Load necessary modules if you have on your HPC
module load scRNA/1.0.4
module load R/4.4.1


# Run pipeline without batch
//path/to/executor/nextflow \
    run /scRNA-seq_nextflow_DEG_pipeline/main.nf \
    --results_dir "/scRNA-seq_nextflow_DEG_pipeline/output_WO_batch/"  \
    --object "/scRNA-seq_nextflow_DEG_pipeline/pbmc_object.Rds" \
    --cond1 "g1" \
    --cond2 "g2" \
    --annotation "letter.idents" \
    --batch_colname "NULL"   \
    --cond_colname "groups" \
    --species "human" \
    --cell_to_filter "Cell_type1, Cell_type_2"
    
# Run pipeline with batch
#//path/to/executor/nextflow \
#    run /scRNA-seq_nextflow_DEG_pipeline/main.nf \
#    --results_dir "/scRNA-seq_nextflow_DEG_pipeline/output_W_batch/"  \
#    --object "/scRNA-seq_nextflow_DEG_pipeline/pbmc_object.Rds" \
#    --cond1 "g1" \
#    --cond2 "g2" \
#    --annotation "letter.idents" \
#    --batch_colname "RNA_snn_res.1"   \
#    --cond_colname "groups" \
#    --species "human" \
#    --cell_to_filter "Cell_type1, Cell_type_2"
