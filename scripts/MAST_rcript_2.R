#!/usr/bin/env Rscript

# Load required packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(limma)
  library(reshape2)
  library(data.table)
  library(knitr)
  library(stringr)
  library(NMF)
  library(rsvd)
  library(RColorBrewer)
  library(MAST)
  library(dplyr)
  library(SeuratData)
  library(patchwork)
  library(SingleCellExperiment)
  library(Seurat)
  library(SeuratObject)
  library(SeuratWrappers)
  library(optparse)  # For command-line argument parsing
  source("//data/Schneider_lab/biopsy_MPN_CML/CML/data/cloned_pipeline/R/save_load_helper.R")

})

# Define command-line option

args <- commandArgs(trailingOnly = TRUE)

# Extract variable values from arguments
object <- args[which(args == "--object") + 1] 
outdir <- args[which(args == "--outdir") + 1] 
cond1 <- args[which(args == "--cond1") + 1] 
cond2 <- args[which(args == "--cond2") + 1] 
anno <- args[which(args == "--annotation") + 1] 
# Import the Seurat object
CML_object <- load_object(object)
print("object_imported")

# Normalize and process Seurat object
counts <- CML_object@assays$RNA@layers$counts
meta_data <- CML_object@meta.data

CML_object <- CreateSeuratObject(counts = counts, meta.data = meta_data)
CML_object <- NormalizeData(CML_object, normalization.method = "LogNormalize", scale.factor = 10000)

print("object_Normalized")

# Conditions to include
include_conditions <- c(cond1, cond2)

# Subset to include only rows with the desired conditions
CML_object_D0WT <- subset(CML_object, stage %in% include_conditions)

# Drop unused levels
CML_object_D0WT$stage <- factor(CML_object_D0WT$stage) #remeber to remove this in real data
CML_object_D0WT$stage <- droplevels(CML_object_D0WT$stage)
print("object_Subseted")

# Split into sub-objects
list_of_subpops_D0vsWT <- SplitObject(CML_object_D0WT, split.by = "integrated_snn_res.0.1")
print("object_Splitted")

# Convert Seurat objects to SingleCellExperiment
for (i in 1:2) {
  list_of_subpops_D0vsWT[[i]] <- as.SingleCellExperiment(list_of_subpops_D0vsWT[[i]])
}
print("object_Singlecellexperimet_converted")





#first function
find_de_MAST_D0vsWT <- function(adata_){
  # create a MAST object
  sca <- SceToSingleCellAssay(adata_, class = "SingleCellAssay")
  print("Dimensions before subsetting:")
  print(dim(sca))
  print("")
  # keep genes that are expressed in more than 10% of all cells
  #sca <- sca[freq(sca)>0.1,]
  print("Dimensions after subsetting:")
  print(dim(sca))
  print("")
  # add a column to the data which contains scaled number of genes that are expressed in each cell
  cdr2 <- colSums(assay(sca)>0)
  colData(sca)$ngeneson <- scale(cdr2)
  # store the columns that we are interested in as factors
  label <- factor(colData(sca)$stage)
  # set the reference level
  label <- relevel(label, cond2)
  colData(sca)$label <- label
  # define and fit the model
  zlmCond <- zlm(formula = ~ngeneson + label ,
                 sca=sca, 
                 method='glm', 
                 ebayes=F) # to speed up calculations
  
  # perform likelihood-ratio test for the condition that we are interested in    
  summaryCond <- summary(zlmCond, doLRT= paste0('label', cond1))
  # get the table with log-fold changes and p-values
  summaryDt <- summaryCond$datatable
  result <- merge(summaryDt[contrast== paste0('label', cond1) & component=='H',.(primerid, `Pr(>Chisq)`)], # p-values
                  summaryDt[contrast== paste0('label', cond1) & component=='logFC', .(primerid, coef)],
                  by='primerid') # logFC coefficients
  # MAST uses natural logarithm so we convert the coefficients to log2 base to be comparable to edgeR
  result[,coef:=result[,coef]/log(2)]
  # do multiple testing correction
  result[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  #result = result[result$FDR<0.05,, drop=F]
  result <- stats::na.omit(as.data.frame(result))
  return(result)
}


list_of_restlts_D0vsWT <- list()



#first comparison

for(i in 1:2){
  tryCatch({
    print(paste0("calculating_MAST_D0vsWT_for_", names(list_of_subpops_D0vsWT)[[i]]))
    list_of_restlts_D0vsWT[[i]] <-find_de_MAST_D0vsWT(list_of_subpops_D0vsWT[[i]])
    names(list_of_restlts_D0vsWT)[[i]] <- names(list_of_subpops_D0vsWT)[[i]]
    colnames(list_of_restlts_D0vsWT[[i]]) <- c("GeneID", "pval", "lfc", "FDR")
    list_of_restlts_D0vsWT[[i]]  <- list_of_restlts_D0vsWT[[i]]  %>%
      mutate(t_stat = (-log10(pval) * sign(lfc)))
    write.csv(list_of_restlts_D0vsWT[[i]], file = paste0(outdir, "/", cond1 ,"vs", cond2, "_", paste0(names(list_of_restlts_D0vsWT)[[i]]), ".csv"))
  }, error = function(e){
    # Handle the error (e.g., print an error message)
    cat("Error in iteration", i, ":", conditionMessage(e), "\n")
  })
}


