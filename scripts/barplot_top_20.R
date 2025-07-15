# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyverse)
  library(ggplot2)  
  library(scales)     
  library(optparse)   
})


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
  tlist[[i]] <- rbind(head(mlist[[i]],n = 20),tail(mlist[[i]],n = 20))
  
}
names(tlist) <- names(mlist)
#loop to generate list of tables for plotting
symmetric_limits <- function (x) 
{
  max <- max(abs(x))
  c(-max, max)
}

for(i in 1:length(tlist)){
  tlist[[i]]$DEG <- c(rep("UP", 20), rep("DOWN", 20))
  
  # Reverse factor levels so positive values are on top
  tlist[[i]]$GeneID <- factor(tlist[[i]]$GeneID, levels = rev(tlist[[i]]$GeneID))
  
  plist[[i]] <- ggplot(tlist[[i]], aes(x = t_stat, y = GeneID, fill = DEG)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    scale_fill_manual(values = c("slateblue4", "red2")) +
    scale_x_continuous(limits = symmetric_limits) +
    theme(legend.position = "bottom", aspect.ratio = 3) +
    ggtitle(names(mlist)[[i]]) +
    geom_vline(xintercept = 0, linetype = "dashed")
}

pdf(file = paste0("barplot_top_20_", cond1 ,"vs", cond2, ".pdf"))
for (i in 1:length(plist)) {
  print(plist[[i]])
  
}
dev.off()

