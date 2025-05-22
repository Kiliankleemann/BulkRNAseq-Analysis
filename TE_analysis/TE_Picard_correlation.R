#### -------- INSTALLING AND LOADING PACKAGES ------- ####
list.of.packages <- c("BiocGenerics","tximport","S4Vectors", "DESeq2", "biomaRt","data.table","readr",
                      "ggplot2", "ggsignif", "ggpubr", "sva", "devtools", "org.Hs.eg.db", "ggnewscale","VennDiagram",
                      "org.Mm.eg.db", "limma","stringr","KEGGREST","ggrepel", "openxlsx", 'splitstackshape',
                      "fgsea","clusterProfiler","pheatmap","ggpubr","cowplot",'dplyr','janitor',"purrr",
                      "RColorBrewer",'AnnotationDbi', 'tidyverse','pheatmap', 'dendextend',"factoextra")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
if(length(new.packages)) BiocManager::install(new.packages)

# Packages loading
invisible(lapply(list.of.packages, library, character.only = TRUE))


###### ----- SETTING WORK DIRECTORY -----#####
setwd("/home/kilian/sciebo - Kleemann, Kilian (kleemann@uni-bonn.de)@uni-bonn.sciebo.de/Aref_DNA_damage_Oct_2023/Kilian_results/")

#Create Directories
dir.create('results/TE_picard_correlation')
dir.create('plots/TE_picard_correlation/')


###### ----- FILE PREPARATION -----#####
# Load metadata file.
sample_data <- read.xlsx("metadata.xlsx", sheet = 1)


# Load TE expression data (TE_family x sample)
TEtranscript_multi_counts <- read.table(paste0('TETranscripts_multi/TEtranscripts_out.cntTable')) %>% row_to_names(row_number = 1)  
TEtranscript_multi_counts2 <- TEtranscript_multi_counts[,-1]
rownames(TEtranscript_multi_counts2) <- TEtranscript_multi_counts[,1]
TEtranscript_multi_counts <- TEtranscript_multi_counts2
colClean1 <- function(TEtranscript_multi_counts){colnames(TEtranscript_multi_counts) <- gsub("Aligned.sortedByCoord.out.bam.T", "", colnames(TEtranscript_multi_counts)); TEtranscript_multi_counts } 
colClean2 <- function(TEtranscript_multi_counts){colnames(TEtranscript_multi_counts) <- gsub("Aligned.sortedByCoord.out.bam.C", "", colnames(TEtranscript_multi_counts)); TEtranscript_multi_counts } 
colClean3 <- function(TEtranscript_multi_counts){colnames(TEtranscript_multi_counts) <- gsub(".*?/", "", colnames(TEtranscript_multi_counts)); TEtranscript_multi_counts } 
TEtranscript_multi_counts <- colClean1(TEtranscript_multi_counts)
TEtranscript_multi_counts <- colClean2(TEtranscript_multi_counts)
TEtranscript_multi_counts <- colClean3(TEtranscript_multi_counts)
TEtranscript_multi_counts <- mutate_all(TEtranscript_multi_counts, function(x) as.numeric(as.character(x)))
TEtranscript_multi_counts_reordered <- TEtranscript_multi_counts %>% select(paste(sample_data$Sample_ID))

TEtranscript_multi_counts_reordered_TE = TEtranscript_multi_counts_reordered %>% rownames_to_column("gene") %>%
  filter(grepl(':',gene)) %>%
  column_to_rownames("gene") %>% t() %>%
  as.data.frame()



# Load Picard metrics (sample x metrics)
file_path <- "Picard_metrics"
file_list <- list.files(path = file_path, pattern = "\\.txt$", full.names = TRUE)
data_picard <- data.frame()
for (i in 1:9){
  df <- read.delim(file_list[i], header=T, comment.char="#")
  df <- df[-c(2:200),]
  df$Sample_ID <- paste0(file_list[i])
  data_picard <- rbind(df, data_picard)
}
#Formatting
data_picard$Sample_ID <- gsub("Picard_metrics/","",data_picard$Sample_ID)
data_picard$Sample_ID <- gsub("_output_rna_metrics.txt","",data_picard$Sample_ID)
#Merge with metadata
data_picard_merge <- left_join(data_picard,sample_data,'Sample_ID')
data_picard_merge$PF_BASES <- as.numeric(data_picard_merge$PF_BASES)
data_picard_merge$PF_ALIGNED_BASES <- as.numeric(data_picard_merge$PF_ALIGNED_BASES)
data_picard_merge <- data_picard_merge %>% column_to_rownames('Sample_ID')
data_picard_merge <- data_picard_merge[,-c(3,8:10,16,23,28:34)]




# Ensure samples match
common_samples <- intersect(rownames(TEtranscript_multi_counts_reordered_TE), rownames(data_picard_merge))
te_expr <- TEtranscript_multi_counts_reordered_TE[common_samples, ]
picard <- data_picard_merge[common_samples, ]

# Perform correlation analysis
results <- expand.grid(
  TE_family = colnames(te_expr),
  Picard_metric = colnames(picard),
  stringsAsFactors = FALSE) %>%
  mutate(Spearman_rho = map2_dbl(TE_family, Picard_metric, ~ cor(te_expr[[.x]], picard[[.y]], method = "spearman")),
    p_value = map2_dbl(TE_family, Picard_metric, ~ cor.test(te_expr[[.x]], picard[[.y]], method = "spearman")$p.value))

# FDR correction
results <- results %>% mutate(FDR = p.adjust(p_value, method = "BH"))

# Save results
write.xlsx(results, "results/TE_picard_correlation/TE_Picard_correlation_results.xlsx")


#Plotting significant result
results_padj05 <- results %>% filter(FDR < 0.1)
heatmap_data <- results_padj05 %>%
  select(TE_family, Picard_metric, Spearman_rho) %>%
  pivot_wider(names_from = Picard_metric, values_from = Spearman_rho) %>%
  column_to_rownames("TE_family")

heatmap_data[is.na(heatmap_data)] <- 0  # or use a small neutral value like 1e-6


phm <- pheatmap( mat = heatmap_data,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Spearman Correlation: TE Expression vs Picard Metrics",
  fontsize_row = 8,
  fontsize_col = 8,
  cellheight = 8,
  cellwidth = 8,)

ggsave(paste0("plots/TE_picard_correlation/TE_Picard_padj05.pdf"),
       phm,
       width = 15,
       height = 40,
       dpi = 300)

#Export Heatmap for correlation per Family
heatmap_data_family <- heatmap_data %>% rownames_to_column('gene')
heatmap_data_family <- cSplit(heatmap_data_family, "gene", ":")

families_unique <- c("LTR"   ,    "DNA"  ,     "SINE"   ,     "LINE")
Pval_cutoff <- "padj0.1"
dir.create(paste0('plots/heatmaps/','TE_Picard_correlation' ,"_",Pval_cutoff))
heatmap_dir <- paste0('plots/heatmaps/','TE_Picard_correlation' ,"_",Pval_cutoff,"/")

for (TE in families_unique) {
  heatmap_data_family_final <- heatmap_data_family %>% 
    filter(gene_3 %in% TE) 
  heatmap_data_family_final$gene <- paste0(heatmap_data_family_final$gene_1,":", heatmap_data_family_final$gene_2,":", heatmap_data_family_final$gene_3)
  heatmap_data_family_final <- heatmap_data_family_final %>% 
    select(!c("gene_1","gene_2", "gene_3")) %>%
    column_to_rownames('gene')
  height_heatmap <- as.numeric(paste0(nrow(heatmap_data_family_final)))
  phm_Fam <- pheatmap(heatmap_data_family_final,
                      color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                      breaks = seq(-2, 2, by = 0.1),
                      kmeans_k = NA,
                      cluster_rows = T,
                      cutree_row = 1,
                      cluster_cols = F,
                      #cutree_cols = 4,
                      #gaps_col = gaps,              
                      legend = TRUE,
                      show_rownames = T,
                      main = paste(TE),
                      #cellwidth = 25,
                      cellheight = 9,
                      treeheight_col = 0,
                      treeheight_row = 50,
                      border_color = 'NA',
                      fontsize = 8,
                      scale = 'row')
  ggsave(paste0(heatmap_dir, TE,'_labelled.pdf'),
         phm_Fam,
         height = height_heatmap,
         width = 20,
         units = 'cm',
         limitsize = FALSE,
         dpi = 300)
}


#Plotting only for intronic reads
# Filter significant TE families for PCT_INTRONIC_BASES
sig_intronic <- results %>%
  filter(Picard_metric == "PCT_INTRONIC_BASES", !is.na(FDR), FDR < 0.05)

sig_te_list <- sig_intronic$TE_family 

dir.create("plots/TE_picard_correlation/PCT_INTRONIC_BASES/")
pointSize <- 20
for (te in sig_te_list) {
  df_plot <- data.frame(
    TE_expr = te_expr[[te]],
    PCT_INTRONIC_BASES = picard$PCT_INTRONIC_BASES,
    Sample_ID = rownames(picard))
  p <- ggplot(df_plot, aes(x = TE_expr, y = PCT_INTRONIC_BASES)) +
    geom_point(color = "black", size = 2) +
    geom_smooth(data = df_plot, method = "gam", se=TRUE, color="grey", formula = y ~ x) +
    stat_cor(data = df_plot, color = 'black', label.x.npc = "left", label.y.npc = "top") +
    labs(title = paste("Correlation: ", te, "vs. PCT_INTRONIC_BASES"),
      x = paste("Expression of", te),
      y = "PCT_INTRONIC_BASES") +
    theme(text = element_text(size = pointSize, colour = "black"),
          rect = element_blank(),
          line = element_line(size = lineWidth, colour = "black"),
          plot.title  = element_text(color="black", size=pointSize),
          axis.title.x = element_text(size = pointSize , colour = "black"),
          axis.title.y = element_blank(),
          axis.text.x  = element_text(size = pointSize , colour = "black", vjust=0.1),
          axis.text.y  = element_text(size = pointSize , colour = "black"),
          axis.ticks = element_line(size = lineWidth, colour = "black"),
          axis.line = element_line(size = lineWidth, colour = "black"),
          legend.position = "right",
          legend.title = element_text(size = pointSize , colour = "black"),
          legend.text = element_text(size = pointSize , colour = "black"),
          legend.key.height = unit(1, "cm"),
          legend.key.width = unit(0.5, "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.margin=unit(.5, "lines"))
  
  ggsave(filename = paste0("plots/TE_picard_correlation/PCT_INTRONIC_BASES/", gsub("[:/]", "_", te), ".pdf"),
         plot = p, width = 6, height = 4)
}




### Plotting overlapp between DETEs and correlated TEs
#Data directories:
dir_1 <- 'results/TE_picard_correlation/TE_Picard_correlation_results.xlsx'
dir_2 <- 'results/TEtranscript_multi_all_samples/DEGene_statistics_padj05.xlsx'

#Names for datasets
gene_set_comparison_name <- "DETEs_vs_Picard_Intron_correlated_TEs"
name_1 <- "TE_Picard_intron_correlation_padj05"
name_2 <- "DEGene_statistics_padj05"

#Data 1
results_1 <- read.xlsx(dir_1)
results_1_genes <- results_1 %>% filter(Picard_metric == "PCT_INTRONIC_BASES", !is.na(FDR), FDR < 0.05) %>% pull(TE_family)

#Data 2
results_2 <- read.xlsx(dir_2)
results_2_genes <- results_2 %>% filter(grepl(':',gene)) %>% pull(gene)

# #Calculate overlap significance - Chi squared test
# nrow(data_result1)
# nrow(overlap_genes)
# x = 60      #Overlap 
# m = 785     #data_result1
# n = 11878   #Total number of instances tested i.e. total RNA editing sites
# k = 728     #data_result2
# hyper_overlap_stat <- phyper(x, m, n, k, lower.tail = FALSE)
# hyper_overlap_stat

# Make venn diagram for overlap
# Chart
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
dev.off()
venn <- venn.diagram(
  x = list(results_1_genes, results_2_genes),
  category.names = c(paste0(name_1),paste0(name_2)),
  #filename = paste0('plots/venn/',gene_set_comparison_name,'.pdf'),
  filename = NULL,
  output=TRUE,
  # Output features
  imagetype = "png",
  height = 900, 
  width = 900, 
  resolution = 600,
  compression = "lzw",
  #Circles
  lwd = 1,
  col=c("red","purple" ), #'orchid3'), 
  fill = c(alpha("red" ,0.5),alpha("purple",0.5)), #alpha('orchid3',0.5)),#
  # Numbers
  cex = 3,
  #fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 0.5,
  cat.default.pos = "outer",
  cat.pos = c(0.2,0.2),
  cat.dist = c(0.3,0.3),
  #cat.fontface = "bold",
  cat.col = c("red","purple"), #"orchid3"),
  cat.fontfamily = "sans")
grid.draw(venn)

dev.off()
pdf(file=paste0('plots/venn/',gene_set_comparison_name,'.pdf'))
grid.draw(venn)
dev.off()


###### ----- DONUT PLOTS ----- #####
dir.create(paste0('plots/donut/',gene_set_comparison_name))
donut_dir <- paste0('plots/donut/',gene_set_comparison_name,"/")

#Barplot showing percentage of families in overlap
overlap_genes <-  Reduce(intersect,list(results_1_genes, results_2_genes))
overlap_genes <- overlap_genes %>% as.data.frame()
overlap_genes <- cSplit(overlap_genes, ".", ":")

# Piechart of overlap
donut_data <- overlap_genes

hsize <- 2

donut_data <- donut_data %>% 
  group_by(._3) %>% # Variable to be transformed
  count() %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc)) %>%
  arrange(desc(._3)) %>% ## arrange in the order of the legend
  mutate(text_y = cumsum(n) - n/2) ### calculate where to place the text labels

# Basic Donut of TE elements
piechart <- ggplot(donut_data, aes(x=hsize , y= n, fill=._3)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + 
  geom_label(aes(label = paste(labels), y = text_y),
             colour = "white", fontface = "bold", size = pointSize *0.1, nudge_x = 1,
             show.legend = FALSE) +
  labs(title = paste(paste())) +
  guides(fill=guide_legend(title=paste0(''))) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        plot.title = element_text(color="black", hjust = , size=pointSize *0.7),
        legend.title = element_text(size = pointSize*0.5 , colour = "black"),
        legend.text = element_text(size = pointSize*0.5 , colour = "black"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())+
  #scale_fill_manual(values = c('blue3','red3'))
  scale_fill_brewer(palette = "Paired") 
piechart

ggsave(paste0(donut_dir,gene_set_comparison_name, '.pdf'),
       plot = piechart,
       device = NULL,
       path = NULL,
       #scale = 1,
       width = 15,
       height = 15,
       units = c("cm"),
       dpi = 600)

