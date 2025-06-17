#!/usr/bin/env Rscript

# Load required packages
# Attempt to load required packages and source an R script
tryCatch({
  suppressPackageStartupMessages({
    library(ggplot2)
    library(limma)
    library(reshape2)
    library(data.table)
    library(knitr)
    library(stringr)
    library(lme4)
    library(rsvd)
    library(RColorBrewer)
    library(MAST)
    library(dplyr)
    library(SeuratData)
    library(patchwork)
    library(SingleCellExperiment)
    library(Seurat)
    library(SeuratObject)
    library(optparse)
  })
}, error = function(e) {
  # Handle the error (e.g., log or print an error message)
  cat("Error during package loading or script sourcing:", conditionMessage(e), "\n")
  # You can also log the error to a file or a monitoring system
})


suppressPackageStartupMessages(library(archive))

save_object <- function(object, file_name, file_format="zstd"){

  stopifnot(file_format %in% c("zstd", "lz4", "gzip", "bzip2", "xz", "nocomp"))

  file_name.tmp <- paste0(file_name, ".tmp")

  if(file_format %in% "nocomp"){
    saveRDS(object = object, file = file_name.tmp, compress = FALSE)
    file.rename(file_name.tmp, file_name)
    return(invisible(NULL))
  }

  if(file_format %in% c("zstd", "lz4")){
    con <- archive::file_write(file = file_name.tmp, filter = file_format)
    open(con)
    saveRDS(object = object, file = con)
    close(con)
    file.rename(file_name.tmp, file_name)
  }else{
    saveRDS(object = object, file = file_name.tmp, compress = file_format)
    file.rename(file_name.tmp, file_name)
  }
}

load_object <- function(file_name){
  con <- archive::file_read(file = file_name)
  res <- readRDS(file = con)
  close(con)
  return(res)
}


# use_tools: if True, use the tools directory to load the object
#            if False, use the save_dir
seutools_partition <- function(scrna, partition, save_dir, allinone=FALSE, use_tools=FALSE){
  #scrna@tools, store into partition if allinone is FALSE
  out <- NULL
  if(allinone == FALSE){
    #assertthat::assert_that(scrna@tools$allinone == FALSE)
    if(!(partition %in% names(scrna@tools))){
      stop("partition not found in scrna@tools")
    }
    if(endsWith(scrna@tools[[partition]], "Rds")){
      if(use_tools == TRUE){
        out <- load_object(scrna@tools[[partition]])
      }else{
        out <- load_object(file.path(save_dir, "partition", glue::glue("{partition}.Rds")))
      }
    }else{
      stop(glue::glue("The {partition}.Rds is not existing!"))
    }
  }else{
    out <- scrna@tools[[partition]]
  }

  return(out)
}



seu_assay <- function(scrna, assay, save_dir, allinone=FALSE, use_tools=FALSE){
  #scrna@tools, store into assay if allinone is FALSE
  out <- NULL
  if(allinone == FALSE){
    #assertthat::assert_that(scrna@tools$allinone == FALSE)
    if(!((assay %in% names(scrna@tools$assay_info)) | (assay %in% names(scrna@assays)))){
      stop("assay not found in scrna@tools$assay_info or scrna@assay")
    }
    if(endsWith(scrna@tools$assay_info[[assay]], "Rds")){
      if(use_tools == TRUE){
        assay_data <- load_object(scrna@tools$assay_info[[assay]])
        assertthat::assert_that(all(colnames(scrna) %in% colnames(assay_data$assay)))
        scrna[[assay]] <- subset(assay_data$assay, cells=colnames(scrna))
        rm(assay_data)
      }else{
        assay_data <- load_object(file.path(save_dir, "assays", glue::glue("{assay}.Rds")))
        assertthat::assert_that(all(colnames(scrna) %in% colnames(assay_data$assay)))
        scrna[[assay]] <- subset(assay_data$assay, cells=colnames(scrna))
        rm(assay_data)
      }
    }else{
      stop(glue::glue("The {assay}.Rds is not existing!"))
    }
  }
  gc()
  return(scrna)
}


# Define command-line option

args <- commandArgs(trailingOnly = TRUE)

# Extract variable values from arguments
object <- args[which(args == "--object") + 1] 
outdir <- args[which(args == "--outdir") + 1] 
cond1 <- args[which(args == "--cond1") + 1] 
cond2 <- args[which(args == "--cond2") + 1] 
anno <- args[which(args == "--annotation") + 1]
cond_colname <- args[which(args == "--cond_colname") + 1]
batch_colname <- args[which(args == "--batch_colname") + 1]
if (is.null(batch_colname) || batch_colname == "" || batch_colname == "NULL") {
  use_batch <- FALSE
} else {
  use_batch <- TRUE
}
# Import the Seurat object
#Seurat_object <- readRDS(object)

Seurat_object <- load_object(object)
print("object_imported")

# Normalize and process Seurat object
counts <- Seurat_object@assays$RNA@counts
meta_data <- Seurat_object@meta.data

Seurat_object <- CreateSeuratObject(counts = counts, meta.data = meta_data)
Seurat_object <- NormalizeData(Seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)

#change column name of the condations from whatever it is called  to stage
column_index <- which(colnames(Seurat_object@meta.data) == cond_colname)
colnames(Seurat_object@meta.data)[column_index] <- "stage"

#change column name of the batch from whatever it is called  to "name"
if (use_batch && batch_colname %in% colnames(Seurat_object@meta.data)) {
  column_index <- which(colnames(Seurat_object@meta.data) == batch_colname)
  colnames(Seurat_object@meta.data)[column_index] <- "name"
}
#change column name of the cell_type from whatever it is called  to "annotation"
column_index <- which(colnames(Seurat_object@meta.data) == anno)
colnames(Seurat_object@meta.data)[column_index] <- "annotation"


print("object_Normalized")

# Conditions to include
include_conditions <- c(cond1, cond2)

# Subset to include only rows with the desired conditions
Seurat_object_D0WT <- subset(Seurat_object, stage %in% include_conditions)

#remove cells that has no cells in one of the conditions
table_cells <- as.data.frame.matrix(table(Seurat_object_D0WT$annotation, Seurat_object_D0WT$stage))
df_no_zeros <- table_cells[apply(table_cells, 1, function(row) all(row > 1)), ]
no_zero_cells <- rownames(df_no_zeros)
Seurat_object_D0WT <- SetIdent(Seurat_object_D0WT, value = "annotation")
Seurat_object_D0WT <- subset(Seurat_object_D0WT, idents = no_zero_cells)
print("object_cleaned_from_zero_cells")

# Drop unused levels
Seurat_object_D0WT$stage <- factor(Seurat_object_D0WT$stage) #remeber to remove this in real data
Seurat_object_D0WT$stage <- droplevels(Seurat_object_D0WT$stage)
print("object_Subseted")

# Split into sub-objects
list_of_subpops_D0vsWT <- SplitObject(Seurat_object_D0WT, split.by = "annotation")
print("object_Splitted")

# Convert Seurat objects to SingleCellExperiment
for (i in 1:length(list_of_subpops_D0vsWT)) {
  list_of_subpops_D0vsWT[[i]] <- as.SingleCellExperiment(list_of_subpops_D0vsWT[[i]])
}
print("object_Singlecellexperimet_converted")





#first function
find_de_MAST_D0vsWT <- function(adata_) {
  sca <- SceToSingleCellAssay(adata_, class = "SingleCellAssay")
  sca <- sca[freq(sca) > 0.1, ]
  cdr2 <- colSums(assay(sca) > 0)
  colData(sca)$ngeneson <- scale(cdr2)
  label <- factor(colData(sca)$stage)
  label <- relevel(label, cond2)
  colData(sca)$label <- label
  
  if (use_batch && "name" %in% colnames(colData(sca))) {
    replicate <- factor(colData(sca)$name)
    colData(sca)$replicate <- replicate
    zlmCond <- zlm(formula = ~ngeneson + label + (1 | replicate),
                   sca = sca,
                   method = 'glmer',
                   ebayes = FALSE)
  } else {
    zlmCond <- zlm(formula = ~ngeneson + label,
                   sca = sca,
                   method = 'glm',
                   ebayes = FALSE)
  }
  
  summaryCond <- summary(zlmCond, doLRT = paste0('label', cond1))
  summaryDt <- summaryCond$datatable
  result <- merge(
    summaryDt[contrast == paste0('label', cond1) & component == 'H', .(primerid, `Pr(>Chisq)`)],
    summaryDt[contrast == paste0('label', cond1) & component == 'logFC', .(primerid, coef)],
    by = 'primerid'
  )
  result[, coef := result[, coef] / log(2)]
  result[, FDR := p.adjust(`Pr(>Chisq)`, 'fdr')]
  result <- stats::na.omit(as.data.frame(result))
  return(result)
}

list_of_restlts_D0vsWT <- list()



#first comparison

for(i in 1:length(list_of_subpops_D0vsWT)){
  tryCatch({
    print(paste0("calculating_MAST_", cond1, "vs" , cond2, "_for_", names(list_of_subpops_D0vsWT)[[i]]))
    list_of_restlts_D0vsWT[[i]] <-find_de_MAST_D0vsWT(list_of_subpops_D0vsWT[[i]])
    names(list_of_restlts_D0vsWT)[[i]] <- names(list_of_subpops_D0vsWT)[[i]]
    colnames(list_of_restlts_D0vsWT[[i]]) <- c("GeneID", "pval", "lfc", "FDR")
    list_of_restlts_D0vsWT[[i]]  <- list_of_restlts_D0vsWT[[i]]  %>%
      mutate(t_stat = (-log10(pval) * sign(lfc)))
    write.csv(list_of_restlts_D0vsWT[[i]], file = paste0( cond1 ,"vs", cond2, "_", names(list_of_restlts_D0vsWT)[[i]], ".csv"))
  }, error = function(e){
    # Handle the error (e.g., print an error message)
    cat("Error in iteration", i, ":", conditionMessage(e), "\n")
  })
}


