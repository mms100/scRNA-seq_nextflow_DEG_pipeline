# Differential gene expresssion nextflow pipeline for seurat V4
A repository establishing a nextflow pipeline for applying MAST DEG analysis (Finak et al., 2015) 

input: Seurat object with Raw counts

outputs: 

1- csv files of Differentailly expressed genes (DEGs) per each cell type

2- csv files of DEGs after removal of mitochondrial, hemoglobin, Immunoglobulins/Plasma genes

3- Volcano plots

4- barplot for cells with highest DEGs number

5- barplot ranking genes using t_stat = -log10(pval) x sign(lfc)

# Parameters list

//scRNA-seq_nextflow_DEG_pipeline/nextflow \
    run /scRNA-seq_nextflow_DEG_pipeline/main.nf \
    --results_dir "/scRNA-seq_nextflow_DEG_pipeline/output_WO_batch/"  \
    --object "/scRNA-seq_nextflow_DEG_pipeline/pbmc_object.Rds" \
    --cond1 "g1" \
    --cond2 "g2" \
    --annotation "letter.idents" \
    --batch_colname "NULL"   \
    --cond_colname "groups" \
    --species "human"
    

**#parameters list**

main.nf = the pipeline protocol

--results_dir = path to output

--object = path for the seurat object

--cond1 = the first condition in the pariwise comparison

--cond2 = the second condition in the pariwise comparison

--annotation = the celltype column in the seurat metadata 

--batch_colname = the column name that contains batch info

--cond_colname = the column name that contains conditions 

--species = "human" / "mouse"



**Note:** no sapce or special character is allowed in cond_colname ( cond1, cond2)
**Example:** 

<pre>
instead of --cond_colname= "cond_1" 
use        --cond_colname= "cond1" 
</pre>

# Running using slurm executor

sbatch slurm_job.sh 

# prerequisite 

**nextflow version** 24.04.2.5914

## Reproducing the R/4.4.1 Environment

# following packages are needed
Bioconductor Version: V(3.20)
library(tidyverse) V(2.0.0)
library(scales) V(1.4.0)
library(SingleCellExperiment) V(1.28.1)
library(Seurat) V(4.4.0)
library(SeuratObject) V(4.1.4)
library(MAST) V(1.32.0)
library(EnhancedVolcano) V(1.24.0)
library(optparse) V(1.8.2)
library(archive) V(1.1.13)

# note  please install V4 of seurat using the following:
remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

####



<img width="467" alt="image" src="https://github.com/user-attachments/assets/5d250dc7-2849-4340-a1fe-f0325891685b" />


Finak, G., McDavid, A., Yajima, M., Deng, J., Gersuk, V., Shalek, A. K., Slichter, C. K., Miller, H. W., McElrath, M. J., Prlic, M., Linsley, P. S., & Gottardo, R. (2015). MAST: A flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data. Genome Biology, 16(1), 278. https://doi.org/10.1186/s13059-015-0844-5

Butler, A., Hoffman, P., Smibert, P., Papalexi, E., & Satija, R. (2018). Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nature Biotechnology, 36(5), 411–420. https://doi.org/10.1038/nbt.4096


