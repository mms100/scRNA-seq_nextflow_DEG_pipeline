#!/usr/bin/env Rscript

# ---- Load Required Packages ----
suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
})

# ---- Setup Arguments ----
option_list <- list(
  make_option(c("--input_dir"), type="character", default=NULL, help="Input directory with CSVs"),
  make_option(c("--outdir"), type="character", default=NULL, help="Output directory for filtered CSVs"),
  make_option(c("--species"), type="character", default="human", help="Species: 'human' or 'mouse'")
)
opt <- parse_args(OptionParser(option_list=option_list))

# ---- Validate Species ----
species <- tolower(opt$species)
if (!species %in% c("human", "mouse")) {
  stop("Species must be either 'human' or 'mouse'")
}

message(paste(">>> Filtering DEG tables for species:", species, "..."))

csv_files <- list.files(path = opt$input_dir, pattern = "\\.csv$", full.names = TRUE)

for (f in csv_files) {
  df <- read.csv(f)
  all_genes <- df$GeneID
  
  # ---- Define Regex Patterns based on Species ----
  if (species == "human") {
    # Patterns for Human (ALL CAPS)
    # MT- = Mito; RPS/RPL = Ribosomal; HB = Hemoglobin; IG = Immunoglobulins/Plasma
    exclude_pattern <- "^(MT-|RPS|RPL|HBA|HBB|HBD|HBG|HBE|HBZ|IG[KLC][A-Z]|IGH[A-Z])"
    
    # Filter genes
    filtered_genes <- all_genes[!grepl(exclude_pattern, all_genes)]
    
  } else if (species == "mouse") {
    # Patterns for Mouse (Sentence Case)
    # mt- = Mito; Rps/Rpl = Ribosomal; Hb = Hemoglobin; Gm = LncRNA/Pseudogenes; Ig = Immunoglobulins
    exclude_pattern <- "^(mt-|Rps|Rpl|Hba|Hbb|Hbq|Gm|Ig[klc][a-z]|Igh[a-z])"
    
    # Filter genes
    filtered_genes <- all_genes[!grepl(exclude_pattern, all_genes)]
    
    # Remove Rik suffixes & specific car/hemoglobin-associated artifacts
    filtered_genes <- filtered_genes[!grepl("Rik$|Rik2$", filtered_genes)]
    filtered_genes <- filtered_genes[!grepl("^Car[1-2]$", filtered_genes)]
  }
  
  # Apply filter and write output
  df_filtered <- df[df$GeneID %in% filtered_genes, ]
  
  # Rename the file to avoid collision
  new_name <- paste0('filtered_', basename(f))
  write.csv(df_filtered, file = file.path(opt$outdir, new_name), row.names = FALSE)
}

message(">>> Filtering completed successfully ✅")