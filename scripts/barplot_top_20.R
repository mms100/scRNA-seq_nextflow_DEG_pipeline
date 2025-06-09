
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
library(RColorBrewer)
library(data.table)
library(readr)
library(optparse)  # For command-line argument parsing
library(scales)

#first capture the arguments from the shell
args <- commandArgs(trailingOnly = TRUE)


inputdir <- args[which(args == "--input_dir_1") + 1] # Get the value after '--var2'
outdir <- args[which(args == "--outdir_4") + 1] # Get the value after '--var2'
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

filtered_list_seurat_object_COND1vsCOND2 <- list()

#import the CSV files
filtered_list_seurat_object_COND1vsCOND2 <- lapply(csv_files, read.csv)


for(i in 1:length(names_subcluster_3)){
  names(filtered_list_seurat_object_COND1vsCOND2)[[i]] <- names_subcluster_3[[i]]
}



#start to prepare the plot


mlist <- filtered_list_seurat_object_COND1vsCOND2
tlist <- list()
plist <- list()


for(i in 1:length(mlist)){
  mlist[[i]] <- mlist[[i]][!is.infinite(mlist[[i]]$t_stat), ]
  mlist[[i]] <- mlist[[i]][order(mlist[[i]]$t_stat,decreasing = T),]
  top_n <- min(20, nrow(mlist[[i]]))  # Ensure it doesn't exceed the number of rows
  bottom_n <- min(20, nrow(mlist[[i]])) 
  top_entries <- head(mlist[[i]], n = top_n)# Select the top 
  bottom_entries <- tail(mlist[[i]], n = bottom_n)# Select the bottom entries
  tlist[[i]] <- rbind(head(mlist[[i]],n = 20),tail(mlist[[i]],n = 20))
  color_vector <- c(rep("red2", top_n), rep("slateblue4", bottom_n))
  tlist[[i]][, 2] <- factor(tlist[[i]][, 2], levels = rev(tlist[[i]][, 2]))

}

names(tlist) <- names(mlist)



plist <- list()
# Set up color vectors based on the number of entries



# Create the bar plot with the appropriate color vector

for(i in 1:length(tlist)){
  pdf(file = paste0("barplot_top_20", cond1 ,"vs", cond2, "_", paste0(names(tlist)[[i]]), ".pdf"))
  barplot(
    tlist[[i]]$t_stat,
    names.arg = tlist[[i]][, 2],
    horiz = TRUE,
    las = 2,
    col = color_vector,
    border = NA,
    main = names(tlist)[i],
    xlab = "t_stat",
    xlim = c(max(tlist[[i]]$t_stat), min(tlist[[i]]$t_stat )),
    cex.names = 0.5  # Replace lower_limit and upper_limit with your desired limits
  )
  dev.off()
  
}


