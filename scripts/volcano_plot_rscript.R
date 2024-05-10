
library(EnhancedVolcano)
library(dplyr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(SeuratData)
library(patchwork)
library(SingleCellExperiment)
library(pheatmap)
library(SeuratObject) 
library(SeuratWrappers)
library(tidyverse)
library(cowplot)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)
library(readr)
library(optparse)  # For command-line argument parsing
library(scales)

#first capture the arguments from the shell
args <- commandArgs(trailingOnly = TRUE)


inputdir <- args[which(args == "--input_dir_1") + 1] # Get the value after '--var2'
outdir <- args[which(args == "--outdir_3") + 1] # Get the value after '--var2'
cond1 <- args[which(args == "--cond1") + 1] # Get the value after '--var2'
cond2 <- args[which(args == "--cond2") + 1] # Get the value after '--var2'




#first comparison
#read the csv files
csv_files <- list.files(path = inputdir, pattern = "\\.csv$", full.names = TRUE)


#generate a vector with cluster names

names_subcluster <- strsplit(csv_files, "/")
names_subcluster_2 <- list()
split_strings <- list()
names_subcluster_3 <- character()

print("str_splitted")
str(names_subcluster)

for(i in 1:length(names_subcluster)){
  names_subcluster_2[i] <- strsplit(names_subcluster[[i]][length(names_subcluster[[1]])], ".csv")
  split_strings[i] <- strsplit(names_subcluster_2[[i]], "_")
  names_subcluster_3[[i]] <- paste(split_strings[[i]][-1], collapse = "_")
}

print("str_captured")

filtered_list_CML_D90vsWT <- list()

#import the CSV files
filtered_list_CML_D90vsWT <- lapply(csv_files, read.csv)


for(i in 1:length(names_subcluster_3)){
  names(filtered_list_CML_D90vsWT)[[i]] <- names_subcluster_3[[i]]
}

#include the top 20 lfc genes in each cell type

#order genes
tlist <- list()
tlist_up <- list()
tlist_down <- list()
mlist <- list()

for(i in 1:length(filtered_list_CML_D90vsWT)){
  tlist[[i]] <- filtered_list_CML_D90vsWT[[i]] 
  tlist[[i]] <- tlist[[i]][order(tlist[[i]]$lfc,decreasing = T),]
  tlist_up[[i]] <- tlist[[i]] %>% filter(FDR <= 0.05 ) %>% filter(lfc >= 0.25 )
  tlist_down[[i]] <- tlist[[i]] %>% filter(FDR <= 0.05 ) %>% filter(lfc <=-0.25)
  mlist[[i]] <- rbind(head(tlist_up[[i]],n = 20),tail(tlist_down[[i]],n = 20))
  names(tlist)[[i]] <- names(filtered_list_CML_D90vsWT)[[i]]
  names(mlist)[[i]] <- names(filtered_list_CML_D90vsWT)[[i]]
}

#plot volcano plot
plist <- list()

for(i in 1:length(filtered_list_CML_D90vsWT)){
  list_of_comparisons <- as.data.frame(table(sign(filtered_list_CML_D90vsWT[[i]]$lfc)))
  plist[[i]] <- EnhancedVolcano(tlist[[i]],
                                lab = tlist[[i]]$GeneID,
                                selectLab = mlist[[i]]$GeneID,
                                axisLabSize = 8,
                                caption = NULL,
                                titleLabSize = 10,
                                subtitleLabSize = 8,
                                x = 'lfc',
                                pCutoff = 0.05,
                                y = 'FDR',
                                pointSize = 1,
				FCcutoff = 0.25,
                                labSize = 3,
                                ylab = bquote(~-Log[10] ~ italic(FDR)),
                                title = paste0(names(tlist)[[i]]),
                                subtitle =paste0( cond1,"vs", cond2 )   ,
                                max.overlaps=Inf,
                                drawConnectors = TRUE,
                                arrowheads = FALSE,
                                maxoverlapsConnectors = Inf)
}



for(i in 1: length(filtered_list_CML_D90vsWT)){
	tryCatch({
	       	pdf(file = paste0("volcano",  cond1 ,"vs", cond2, "_", paste0(names(tlist)[[i]]), ".pdf"))
		print(plist[[i]])
		dev.off()
	}, error = function(e){
    # Handle the error (e.g., print an error message)
		cat("Error in iteration", i, ":", conditionMessage(e), "\n")
	})
}
                

