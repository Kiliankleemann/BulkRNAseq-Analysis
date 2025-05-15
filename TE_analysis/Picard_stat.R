#### -------- INSTALLING AND LOADING PACKAGES ------- ####
list.of.packages <- c("BiocGenerics","tximport","S4Vectors", "DESeq2", "biomaRt",
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
setwd("/Users/kiliankleemann/sciebo - Kleemann, Kilian (kleemann@uni-bonn.de)@uni-bonn.sciebo.de/Ivana_IFN_KI_experiments/Rayk_Aref_Ivana_MEF_DAXX_October_2024")

#Import metadata
sample_data <- read.xlsx("metadata.xlsx", sheet = 1)


###### ----- PICARD ALIGN STATS -----#####
# Set the directory containing the .txt files
file_path <- "Picard_metrics"

# List all .txt files in the directory
file_list <- list.files(path = file_path, pattern = "\\.txt$", full.names = TRUE)

# Create an empty list to store dataframes
data_picard <- data.frame()
for (i in 1:6){
    df <- read.delim(file_list[i], header=T, comment.char="#")
    df <- df[-c(2:200),]
    data_picard <- rbind(df, data_picard)
  }
  
file_list_name <- gsub("Picard_metrics/","",file_list)
file_list_name <- gsub("_output_rna_metrics.txt","",file_list_name)

#Merge with metadata
data_picard$Sample_ID <- file_list_name
data_picard_merge <- left_join(data_picard,sample_data,'Sample_ID')
data_picard_merge$PF_BASES <- as.numeric(data_picard_merge$PF_BASES)
data_picard_merge$PF_ALIGNED_BASES <- as.numeric(data_picard_merge$PF_ALIGNED_BASES)

data_picard_merge <- data_picard_merge %>% column_to_rownames('Sample_ID')
data_picard_merge <- data_picard_merge[,-c(3,16,28:35,37)]

data_picard_merge <- data_picard_merge %>% rownames_to_column('Sample_ID')

#make Barplots for all conditions
#Per sample
dir.create(paste0('plots/barplots/',"Picard_stats_per_sample"))

lineWidth=0.5
pointSize =30
for (i in picard_stats) {
  column_sym <- sym(i)
  barplot <- ggplot(data_picard_merge, aes(x = Sample_ID, y = !!column_sym)) + 
    geom_bar(position = 'dodge', stat = 'summary', fun = mean, width = 0.7, colour="black") +
    expand_limits(x = 0, y = 0) +
    theme(panel.spacing = unit(1, "lines")) +
    labs(x = NULL, y = c("Count")) +
    ggtitle(str_wrap(paste0(i), width = 30))+
    theme(text = element_text(size = pointSize, colour = "black"),
      rect = element_blank(),
      line = element_line(linewidth = lineWidth, colour = "black"),
      plot.title  = element_text(color="black", size=pointSize),
      axis.title  = element_text(size = pointSize, colour = "black"),
      axis.text.x  = element_text(size = pointSize , colour = "black", angle = 45, hjust = 1),
      axis.text.y  = element_text(size = pointSize , colour = "black"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = pointSize , colour = "black"),
      # legend.key.height = unit(0.1, "cm"),
      # legend.key.width = unit(0.2, "cm"),
      axis.line = element_line(linewidth = lineWidth, colour = "black"),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) #+ scale_y_continuous(expand = c(0, 0, .05, 0))
  
  pdf(file = paste0('plots/barplots/',"Picard_stats_per_sample",'/',i,'.pdf'), pointsize = 10, width = 10, height= 10)
  print(barplot)
  dev.off()
}



#Per condition
dir.create(paste0('plots/barplots/',"Picard_stats_per_condition"))
data_picard_merge$Condition <- data_picard_merge$Genotype

for (i in picard_stats) {
  formula <- as.formula(paste(i, "~ Condition"))
  comparisons <- compare_means(
    data = data_picard_merge,
    formula = formula,
    method = "t.test",
    p.adjust.method = "BH") %>% filter(p<0.05)
  
  comparison_list_sign <- comparisons %>% mutate(comparison_list = map2(group1, group2,c)) %>% pull(comparison_list)
  
  column_sym <- sym(i)
  
  barplot <- barplot <- ggplot(data_picard_merge, aes(x = Condition, y = !!column_sym, fill = Condition)) + 
    geom_bar(position = 'dodge', stat = 'summary', fun = mean, width = 0.7, colour="black") +
    geom_errorbar(stat = 'summary', position = 'dodge', width = 0.4) +
    geom_point(aes(x = Condition), position = position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width=0.9)) +
    expand_limits(x = 0, y = 0) +
    theme(panel.spacing = unit(1, "lines")) +
    labs(x = NULL, y = i) +
    ggtitle(str_wrap(paste0('Stats Per Condition'), width = 30))+
    theme(text = element_text(size = pointSize, colour = "black"),
          rect = element_blank(),
          line = element_line(linewidth = lineWidth, colour = "black"),
          plot.title  = element_text(color="black", size=pointSize*0.5),
          axis.title  = element_text(size = pointSize, colour = "black"),
          axis.text.x  = element_text(size = pointSize , colour = "black", angle = 45, hjust = 1),
          axis.text.y  = element_text(size = pointSize , colour = "black"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = pointSize , colour = "black"),
          # legend.key.height = unit(0.1, "cm"),
          # legend.key.width = unit(0.2, "cm"),
          axis.line = element_line(linewidth = lineWidth, colour = "black"),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
    stat_compare_means(comparisons = comparison_list_sign, method = 't.test', label = "p.signif" ) +
    scale_fill_brewer(palette = "Paired") + scale_y_continuous(expand = c(0, 0, .05, 0))
  
  pdf(file = paste0('plots/barplots/',"Picard_stats_per_condition",'/',i,'.pdf'), pointsize = 10, width = 5, height= 10)
  print(barplot)
  dev.off()
}


