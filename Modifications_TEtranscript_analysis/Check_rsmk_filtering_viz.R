#### -------- INSTALLING AND LOADING PACKAGES ------- ####
list.of.packages <- c("BiocGenerics","tximport","S4Vectors", "DESeq2", "biomaRt",'data.table',
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
setwd <- ("/home/kilian/sciebo - Kleemann, Kilian (kleemann@uni-bonn.de)@uni-bonn.sciebo.de/Immune_priming_TE/Repeatmasker_modification")

#Creating folders
dir.create('results')
dir.create('plots')


#GTF file to be investigated
file_path <- "/media/kilian/OS/GTF_files_TEtranscript/GRCm38_rmsk_TE_transcript_overlap.gtf"
file_prefix <- "GRCm38_rmsk_TE_transcript_overlap"

# Function to read GTF file
gtf <- fread(file_path, header = FALSE)
colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

# Function to count occurrences of each class ID
count_class_occurrences <- function(gtf) {
  gtf[, class_id := gsub(".*class_id \"(.*?)\".*", "\\1", attribute)]
  class_occurrences <- gtf[, .N, by = class_id]
  return(class_occurrences)
}
# Function to calculate percentage of each class ID
calculate_percentage <- function(class_occurrences) {
  total_records <- sum(class_occurrences$N)
  class_occurrences <- class_occurrences %>%
    mutate(percentage = N / total_records * 100)
  return(class_occurrences)
}
# Main function
gtf <- read_gtf(file_path)
cat("Total number of records in GTF file:", nrow(gtf), "\n\n")
  
class_occurrences <- count_class_occurrences(gtf)
class_percentage <- calculate_percentage(class_occurrences)
  
cat("Summary of Class IDs and Percentage:\n")
print(class_percentage)
  
# Export table as Excel file
write.xlsx(class_percentage, "class_percentage.xlsx")
cat("\nTable exported as 'class_percentage.xlsx'\n")



#Making piechart for all elements 
hsize <- 2
pointSize = 20
lineWidth = 0.5

#Setting up colors: 
palette1 <- brewer.pal(9, "Paired")  # Example palette 1 with 9 colors
palette2 <- brewer.pal(8, "Dark2")  # Example palette 2 with 11 colors

# Combine the two palettes
combined_palette <- c(palette1, palette2)


donut_data <- class_percentage %>% 
  mutate(perc = `N` / sum(`N`)) %>% 
  arrange(perc) %>%
  mutate(labels = paste( sprintf("%.2f%%", N / total_records * 100))) %>%
  arrange(desc(class_id)) %>% ## arrange in the order of the legend
  mutate(text_y = cumsum(N) - N/2) ### calculate where to place the text labels

# Basic Donut of TE elements
piechart <- ggplot(donut_data, aes(x=hsize , y=N, fill=class_id)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + 
  geom_label(aes(label = paste(labels), y = text_y),
             colour = "white", fontface = "bold", size = pointSize *0.1, nudge_x = 1,
             show.legend = FALSE) +
  labs(title = paste('TEs overlapping with transcripts in mm10')) +
  guides(fill=guide_legend(title=paste0('TE Class'))) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        plot.title = element_text(color="black", hjust = -1, size=pointSize, face = "italic"),
        legend.title = element_text(size = pointSize*0.5 , colour = "black"),
        legend.text = element_text(size = pointSize*0.5 , colour = "black"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  scale_fill_manual(values = combined_palette) 

dir.create(paste0("plots/donut"))
dir.create(paste0("plots/donut/",file_prefix))

ggsave(paste0('plots/donut/', file_prefix,'/','TE_Class.pdf'),
       plot = piechart,
       device = NULL,
       path = NULL,
       #scale = 1,
       width = 15,
       height = 15,
       units = c("cm"),
       dpi = 600)



