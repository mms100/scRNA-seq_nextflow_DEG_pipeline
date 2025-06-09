# Differential gene expresssion nextflow pipeline
A repository establishing a nextflow pipeline for applying MAST DEG analysis (Finak et al., 2015)


using pbmc_small object provided from seurat package (Butler et al., 2018)

# Parameters list

//path/to/nextflow_excutor run path/to/main.nf  --object "/path/to/seurat_object.rds" --cond1 "g1" --cond2 "g2" --annotation "RNA_snn_res.0.8" --cond_colname "groups"  --batch_colname "letter.idents" --output_1 "Pbmc_batch"

**#parameters list**

main.nf = the pipeline protocol

--object = path for the seurat object

--cond1 = the first condition in the pariwise comparison

--cond2 = the second condition in the pariwise comparison

--annotation = the celltype column in the seurat metadat

--cond_colname = the column that contains conditions 

--batch_colname = the column that contains batch info.

--output_1 = name of the output repository 

# Running using slurm executor

sbatch slurm_job.sh 

# prerequisite packages to be installed in R studio
## Reproducing the R Environment

This project uses [renv](https://rstudio.github.io/renv/) for reproducible R environments.

To recreate the environment:

1. Install renv (if not already installed):

   install.packages("renv")


2. In the project directory, run:

   renv::restore()

This will install all packages as specified in `renv.lock`.

####



<img width="467" alt="image" src="https://github.com/user-attachments/assets/5d250dc7-2849-4340-a1fe-f0325891685b" />


Finak, G., McDavid, A., Yajima, M., Deng, J., Gersuk, V., Shalek, A. K., Slichter, C. K., Miller, H. W., McElrath, M. J., Prlic, M., Linsley, P. S., & Gottardo, R. (2015). MAST: A flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data. Genome Biology, 16(1), 278. https://doi.org/10.1186/s13059-015-0844-5

Butler, A., Hoffman, P., Smibert, P., Papalexi, E., & Satija, R. (2018). Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nature Biotechnology, 36(5), 411–420. https://doi.org/10.1038/nbt.4096


