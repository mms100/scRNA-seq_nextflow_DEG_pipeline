
#load needed libraries####
# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)       
  library(tidyverse)  
  library(scales)     
  library(optparse)   
})

args <- commandArgs(trailingOnly = TRUE)


inputdir <- args[which(args == "--input_dir_1") + 1] # Get the value after '--var2'
outdir <- args[which(args == "--outdir_2") + 1] # Get the value after '--var2'
organ <- args[which(args == "--params.output_1") + 1] # Get the value after '--var2'
cond1 <- args[which(args == "--cond1") + 1] 
cond2 <- args[which(args == "--cond2") + 1] 

#in my PC
filtered_list <- list()



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
filtered_list <- lapply(csv_files, read.csv)


for(i in 1:length(names_subcluster_3)){
  names(filtered_list)[[i]] <- names_subcluster_3[[i]]
}

#filter to include only significant DEG

for(i in 1:length(names_subcluster_3)){
  names(filtered_list)[[i]] <- names_subcluster_3[[i]]
  filtered_list[[i]] <- filtered_list[[i]] %>% filter(FDR <= 0.05) %>% filter(lfc <= -0.25 | lfc >= 0.25)
}

filtered_list <- Filter(function(df) nrow(df) > 0, filtered_list)



comparisons <- c("Data_CML_D0vsWT")
list_of_comparisons <- list()


for(i in 1){
  list_of_comparisons[[1]]<- data.frame(matrix(nrow = 2, ncol = length(filtered_list)))
  names(list_of_comparisons)[[i]] <- comparisons[i]
}

colnames(list_of_comparisons[[1]]) <- names(filtered_list)



#data_1<- data.frame(matrix(nrow = 2, ncol = length(filtered_list)))
data_1 <- data.frame(matrix(0, nrow = 2, ncol = length(filtered_list)))



for(i in 1:length(filtered_list)) {
  #first comparison
  list_of_comparisons[[1]] <- as.data.frame(table(sign(filtered_list[[i]]$lfc)))
  list_of_comparisons[[1]]$Var1 <- as.numeric(as.character(list_of_comparisons[[1]]$Var1))
  if(!(-1 %in% list_of_comparisons[[1]]$Var1)) {
    list_of_comparisons[[1]] <- rbind(list_of_comparisons[[1]], data.frame(Var1 = -1, Freq = 0))
  }
  if(!(1 %in% list_of_comparisons[[1]]$Var1)) {
    list_of_comparisons[[1]] <- rbind(list_of_comparisons[[1]], data.frame(Var1 = 1, Freq = 0))
  }
  list_of_comparisons[[1]] <- list_of_comparisons[[1]][order(list_of_comparisons[[1]]$Var1), ]
  list_of_comparisons[[1]]$Product <- list_of_comparisons[[1]]$Var1 * list_of_comparisons[[1]]$Freq
  data_1[,i] <- list_of_comparisons[[1]]$Product
  colnames(data_1) <- names(filtered_list)
  }



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
    breaks = pretty_breaks(n = 20), 
    labels = abs
  )+
  labs(x = "cluster", y = "Number of regulated genes with threshold (FDR=0.05, lfc=+/-0.25) ", title = paste0(organ, '_', cond1,'_', cond2)) +
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
  theme(axis.text.x = element_text(  size= 10, angle = 45, vjust = 1, hjust = 1))


dev.off()


