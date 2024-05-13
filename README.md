A repository establishing a nextflow pipeline for applying MAST DEG analysis (Finak et al., 2015)


using pbmc_small object provided from seurat package (Butler et al., 2018)


Running using slurm executor

//path/to/nextflow_excutor run main.nf --object "/path/to/seurat_object.rds" --cond1 "IRF8KO" --cond2 "WT" --annotation "RNA_snn_res.0.1" --cond_colname "stage"

###parameters list
main.nf = the pipeline protocol

--object = path for the seurat object

--cond1 = the first condition in the pariwise comparison

--cond2 = the second condition in the pariwise comparison

--annotation = the celltype column in the seurat metadat

--cond_colname = the column that contains conditions 
####

<img width="670" alt="image" src="https://github.com/mms100/scRNA-seq_nextflow_DEG_pipeline/assets/60142059/aeab1978-17bb-4352-a75d-29329a875c31">


Finak, G., McDavid, A., Yajima, M., Deng, J., Gersuk, V., Shalek, A. K., Slichter, C. K., Miller, H. W., McElrath, M. J., Prlic, M., Linsley, P. S., & Gottardo, R. (2015). MAST: A flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data. Genome Biology, 16(1), 278. https://doi.org/10.1186/s13059-015-0844-5

Butler, A., Hoffman, P., Smibert, P., Papalexi, E., & Satija, R. (2018). Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nature Biotechnology, 36(5), 411–420. https://doi.org/10.1038/nbt.4096


