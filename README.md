# Differential gene expresssion nextflow pipeline for seurat V4
A repository establishing a nextflow pipeline for applying MAST DEG analysis (Finak et al., 2015) 


# Parameters list

//path/to/nextflow_excutor run path/to/main.nf  --object "/path/to/seurat_object.rds" --cond1 "g1" --cond2 "g2" --annotation "RNA_snn_res.0.8" --cond_colname "stage"  --batch_colname "library_name" --output_1 "batch_on"

**#parameters list**

main.nf = the pipeline protocol

--object = path for the seurat object

--cond1 = the first condition in the pariwise comparison

--cond2 = the second condition in the pariwise comparison

--annotation = the celltype column in the seurat metadata 

--cond_colname = the column name that contains conditions 

--batch_colname = the column name that contains batch info

--output_1 = name of the output repository 

**Note:** no sapce or special character is allowed in any of the previously listed metadata ( cond1, cond2, ... etc.)
**Example:** 

**instead of ** --annotation= annotation 1 (B cells 1, T cells 2, ... etc.) 

**use**         --annotation= annotation_1 (B_cells_1, T_cells_2, ... etc.) 

# Running using slurm executor

sbatch slurm_job.sh 

# prerequisite 

**nextflow version** 24.04.2.5914

## Reproducing the R/4.3.2 Environment

This project uses [renv](https://rstudio.github.io/renv/) for reproducible R environments.

To recreate the environment:

1. Install renv (if not already installed):

   install.packages("renv")


2. In the project directory, run:
   
   setwd("Path/to/scRNA-seq_nextflow_DEG_pipeline")
   
   renv::restore()

This will install all packages as specified in `renv.lock`.

####



<img width="467" alt="image" src="https://github.com/user-attachments/assets/5d250dc7-2849-4340-a1fe-f0325891685b" />


Finak, G., McDavid, A., Yajima, M., Deng, J., Gersuk, V., Shalek, A. K., Slichter, C. K., Miller, H. W., McElrath, M. J., Prlic, M., Linsley, P. S., & Gottardo, R. (2015). MAST: A flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data. Genome Biology, 16(1), 278. https://doi.org/10.1186/s13059-015-0844-5

Butler, A., Hoffman, P., Smibert, P., Papalexi, E., & Satija, R. (2018). Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nature Biotechnology, 36(5), 411–420. https://doi.org/10.1038/nbt.4096


