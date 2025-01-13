#### -------- INSTALLING AND LOADING PACKAGES ------- ####
list.of.packages <- c("BiocGenerics","tximport","S4Vectors", "DESeq2", "biomaRt","data.table",
                      "ggplot2", "ggsignif", "ggpubr", "sva", "devtools", "org.Hs.eg.db", 
                      "org.Mm.eg.db", "limma","stringr","KEGGREST","ggrepel", "openxlsx", 'splitstackshape',
                      "fgsea","clusterProfiler","pheatmap","ggpubr","cowplot",'dplyr','janitor',
                      "RColorBrewer",'AnnotationDbi', 'tidyverse','pheatmap', 'dendextend')

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
if(length(new.packages)) BiocManager::install(new.packages)

# Packages loading
invisible(lapply(list.of.packages, library, character.only = TRUE))


###### ----- SETTING WORK DIRECTORY -----#####
setwd("~/sciebo - Kleemann, Kilian (kleemann@uni-bonn.de)@uni-bonn.sciebo.de/Ivana_IFN_KI_experiments/BMDM_IFNa4KI_PE_241217_P2024-266-RNA")

#Create Directories
dir.create('results')
dir.create('plots')
dir.create('plots/PCA')
dir.create('plots/barplots')
dir.create('plots/heatmaps')
dir.create('plots/donut/')
dir.create('plots/volcano/')

##### Import Data and metadata #####


# Define the directory containing FeatureCounts files
dir_path <- "featureCounts_out/"  # Replace with your directory path

# List all .txt files in the directory
file_list <- list.files(dir_path, pattern = "\\.txt$", full.names = TRUE)

# Initialize an empty list to store data frames
count_list <- list()

# Loop through each file and extract counts
for (file in file_list) {
  # Read the file
  count_data <- fread(file, skip = 1)  # Skip the first line (header)
  
  # Extract relevant columns (Geneid and counts for the sample)
  sample_name <- sub("\\.txt$", "", basename(file))  # Remove .txt extension
  count_list[[sample_name]] <- count_data[, .(Geneid, Counts = V7)]  # Adjust V7 if counts are in a different column
}

# Merge all files into a single data frame
count_matrix <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), count_list)

# Rename columns for clarity
colnames(count_matrix)[2:ncol(count_matrix)] <- names(count_list)



