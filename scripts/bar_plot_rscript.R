#prepare the plotting datafame####

#load needed libraries####

library(dplyr)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(SeuratData)
library(patchwork)
library(SingleCellExperiment)
library(pheatmap)
library(SeuratObject) 
library(SeuratWrappers)
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


args <- commandArgs(trailingOnly = TRUE)


inputdir <- args[which(args == "--input_dir_1") + 1] # Get the value after '--var2'
outdir <- args[which(args == "--outdir_2") + 1] # Get the value after '--var2'


#in my PC
filtered_list_seurat_object_COND1vsCOND2 <- list()



#prepare the plotting datafame####

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

#import the CSV files
filtered_list_seurat_object_COND1vsCOND2 <- lapply(csv_files, read.csv)


for(i in 1:length(names_subcluster_3)){
  names(filtered_list_seurat_object_COND1vsCOND2)[[i]] <- names_subcluster_3[[i]]
}

#filter to include only significant DEG

for(i in 1:length(names_subcluster_3)){
  names(filtered_list_seurat_object_COND1vsCOND2)[[i]] <- names_subcluster_3[[i]]
  filtered_list_seurat_object_COND1vsCOND2[[i]] <- filtered_list_seurat_object_COND1vsCOND2[[i]] %>% filter(FDR <= 0.05) %>% filter(lfc <= -0.25 | lfc >= 0.25)
}

filtered_list_seurat_object_COND1vsCOND2 <- Filter(function(df) nrow(df) > 0, filtered_list_seurat_object_COND1vsCOND2)



comparisons <- c("Data_seurat_object_COND1vsCOND2")
list_of_comparisons <- list()


for(i in 1){
  list_of_comparisons[[1]]<- data.frame(matrix(nrow = 2, ncol = length(filtered_list_seurat_object_COND1vsCOND2)))
  names(list_of_comparisons)[[i]] <- comparisons[i]
}

colnames(list_of_comparisons[[1]]) <- names(filtered_list_seurat_object_COND1vsCOND2)



data_1<- data.frame(matrix(nrow = 2, ncol = length(filtered_list_seurat_object_COND1vsCOND2)))


for(i in 1:length(filtered_list_seurat_object_COND1vsCOND2)) {
  #first comparison
  list_of_comparisons[[1]] <- as.data.frame(table(sign(filtered_list_seurat_object_COND1vsCOND2[[i]]$lfc)))
  data_1[,i] <- list_of_comparisons[[1]]$Freq
  colnames(data_1) <- names(filtered_list_seurat_object_COND1vsCOND2)
  
}



#multiplie downregulated values by -1
data_1[1, ] <- data_1[1, ] * -1


data_1$regulation <- c("Down", "UP")

dat_pivoted_1 <- pivot_longer(data_1, cols = 1:ncol(data_1)-1, names_to = "cluster", values_to = "Freq")


#pivot table
dat_pivoted_1 <- dat_pivoted_1[order(dat_pivoted_1$Freq, decreasing = T), ]
dat_pivoted_1$order <- rep(1:nrow(dat_pivoted_1))
exp_cluster <- dat_pivoted_1$cluster
exp_cluster <- exp_cluster[1:length(exp_cluster)]
dat_pivoted_1 <- dat_pivoted_1[dat_pivoted_1$cluster %in% exp_cluster,]


## find the order
list_of_pivoted <- list(dat_pivoted_1)
list_of_tmp <- list()
the_order_list <- list()
for(i in 1){
  list_of_tmp[[i]] <- list_of_pivoted[[i]]%>% 
    filter(regulation == "UP") %>% 
    arrange(Freq)
  the_order_list[[i]] <- list_of_tmp[[i]]$cluster
}


#plotting one by one to check
upper_limit <- max(list_of_pivoted[[1]]$Freq)
lower_limit <- min(list_of_pivoted[[1]]$Freq)
pdf(file= "barplot.pdf")

list_of_pivoted[[1]]%>% 
  ggplot(aes(x = cluster, y = Freq, group = regulation, fill = regulation)) +
  geom_bar(stat = "identity", width = 0.75) +
  coord_flip() +
  scale_x_discrete(limits = the_order_list[[1]]) +
  # another trick!
  scale_y_continuous(
    breaks = pretty_breaks(n = 50),  # Adjust 'n' for desired number of breaks
    labels = abs
  )+
  labs(x = "cluster", y = "Number of regulated genes with threshold (FDR=0.05, lfc=+/-0.25) ", title = comparisons[1]) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill =  "white")
  ) + theme_bw()+
  scale_fill_manual(values=c("#1777B9", "#B93217"),
                    name="",
                    breaks=c("Down", "UP"),
                    labels=c("Down", "UP"))+
  theme(axis.text.y = element_text(face= "bold", size= 10))+
  theme(axis.text.x = element_text(  size= 7, angle = 45))


dev.off()


