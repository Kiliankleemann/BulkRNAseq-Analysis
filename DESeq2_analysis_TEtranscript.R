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
setwd("//")
dir.create('results')
dir.create('plots')
dir.create('plots/PCA')
dir.create('plots/barplots')
dir.create('plots/heatmaps')



###### ----- FILE PREPARATION -----#####
# Create a metadata map using your metadata file.
sample_data <- read.xlsx("metadata.xlsx", sheet = 1)
outliers = c('')
sample_data <- sample_data %>%            
  filter(!Sample_ID %in% outliers)

#Import TEtranscript counts
TEtranscript_multi_counts <- read.table(paste0('TETranscripts_multi/TEtranscripts_out.cntTable')) %>% row_to_names(row_number = 1)  
TEtranscript_multi_counts <- TEtranscript_multi_counts %>% column_to_rownames('gene/TE')
TEtranscript_multi_counts2 <- TEtranscript_multi_counts[,-1]
rownames(TEtranscript_multi_counts2) <- TEtranscript_multi_counts[,1]
TEtranscript_multi_counts <- TEtranscript_multi_counts2

colClean1 <- function(TEtranscript_multi_counts){ colnames(TEtranscript_multi_counts) <- gsub("Aligned.sortedByCoord.out.bam.T", "", colnames(TEtranscript_multi_counts)); TEtranscript_multi_counts } 
colClean2 <- function(TEtranscript_multi_counts){ colnames(TEtranscript_multi_counts) <- gsub("Aligned.sortedByCoord.out.bam.C", "", colnames(TEtranscript_multi_counts)); TEtranscript_multi_counts } 
colClean3 <- function(TEtranscript_multi_counts){ colnames(TEtranscript_multi_counts) <- gsub(".*?multi/", "", colnames(TEtranscript_multi_counts)); TEtranscript_multi_counts } 

TEtranscript_multi_counts <- colClean1(TEtranscript_multi_counts)
TEtranscript_multi_counts <- colClean2(TEtranscript_multi_counts)
TEtranscript_multi_counts <- colClean3(TEtranscript_multi_counts)

#Make numeric
TEtranscript_multi_counts <- mutate_all(TEtranscript_multi_counts, function(x) as.numeric(as.character(x)))


#Arrange count_data columns as ordered in sample_data
TEtranscript_multi_counts_filtered <- TEtranscript_multi_counts %>% select(-contains(outliers))
TEtranscript_multi_counts_reordered <- TEtranscript_multi_counts_filtered %>% select(paste(sample_data$Sample_ID))


#Analyzis parameters
file_prefix = 'TEtranscript_multi_counts_KI_vs_WT'
experiment = TEtranscript_multi_counts_reordered
experiment
sample_data <- sample_data
sample_data


###### ----- DESEQ2 ANALYSIS -----#####
#Designs 
#dds <- DESeqDataSetFromTximport(txi, colData = experiment, design = ~  Genotype) 
dds <- DESeqDataSetFromMatrix(countData = experiment, colData = sample_data, design = ~ Genotype) 

#prefiltering on minimum of 5 reads (NOT REQUIRED AS DESEQ2 OPTIMIZES AND FILTERS AUTOMATICALLY)
dds <- estimateSizeFactors(dds)
# idx <- rowSums(counts(dds, normalized = TRUE)) >= 5
# dds <- dds[idx,]

#Setting the reference level (control group to compare against)
dds@colData@listData$Genotype <- relevel(dds@colData@listData$Genotype, ref = "WT")


## Run DESeq analysis to gather differential expression results
#Run DESeq (LRT)
dds_run <- DESeq(dds, test = 'LRT', reduced = ~ 1)
# Run DESeq (Wald)
#dds_run <- DESeq(dds, betaPrior = F)

# View names of estimated effects
resultsNames(dds_run)

# Create table of effects showing log2fold, standard error, test stats, and p-vals
# LRT
dds_result <- results(dds_run, pAdjustMethod = 'BH', independentFiltering = T, cooksCutoff = T)

# Wald
# dds_result <- results(dds_run, contrast = c('Condition_1','LPS', 'ctrl'), independentFiltering = T, pAdjustMethod = 'BH')
# dds_result <- lfcShrink(dds_run, contrast = c('Condtition_1', 'LPS', 'ctrl'), res = dds_result, type = 'normal')

#dds_result
summary(dds_result, alpha = 0.05)


## View Total number of normalized counts per sample
dds_run <- estimateSizeFactors(dds_run)
DS_norm_counts <- counts(dds_run, normalized = TRUE)
colSums(DS_norm_counts)

## Plot dispersion estimates
#plotDispEsts(dds_run)

#plotCounts(dds_run, 'Fus', intgroup = 'condition')

# Transform counts for data visualization
vst <- vst(dds_run, blind=TRUE)

###### ----- PCA -----#####
# Add nametags
z <- plotPCA(vst, intgroup=c('Genotype'),ntop = 200)

theme_PCA <- theme(aspect.ratio = 1, 
                   panel.background = element_blank(),
                   panel.border=element_rect(fill=NA),
                   panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   strip.background=element_blank(),
                   axis.text.x=element_text(colour="black"),
                   axis.text.y=element_text(colour="black"),
                   axis.ticks=element_line(colour="black"),
                   legend.key=element_blank(),
                   plot.margin=unit(c(1,1,1,1),"line"))

pdf(file = paste0('plots/PCA/', file_prefix,'', '.pdf'), pointsize = 10)
ggplot(z$data,aes(x=z$data$PC1, y=z$data$PC2, )) +
  geom_point(aes(color = Genotype)) +
  theme_PCA +
  labs(title = 'PCA (Top 200 variable genes)',   x=paste(z$labels$x), y=paste(z$labels$y))
dev.off()


pdf(file = paste0('plots/PCA/', file_prefix,'_labelled', '.pdf'), pointsize = 10)
ggplot(z$data,aes(x=z$data$PC1, y=z$data$PC2, label = z$data$name)) +
  geom_point(aes(color = Genotype)) +
  theme_PCA + labs(title = 'PCA (Top 200 variable genes)',   x=paste(z$labels$x), y=paste(z$labels$y)) +
  geom_text_repel(aes(label = sample_data$Sample_ID),max.overlaps = Inf )
dev.off()


###### ----- SIG DATA & EXPORT -----#####
#unwanted_genes = paste(c('^Gm', '^mt-', '^Vmn', '^Rpl', '^Rps', '^Olfr','Rik'), collapse = '|')

# Build significant gene table and extract list of sorted DE genes.
sig_res <- dds_result %>%
  data.frame() %>%
  rownames_to_column(var='gene') %>%    # makes a 'gene' column using column1
  as_tibble %>%
  filter(pvalue < 0.05) %>%
  arrange(pvalue) #%>% 
#filter(!str_detect(gene, unwanted_genes))
sorted_DEGenes <- sig_res$gene

# Export significant gene count data
name_list <- c('gene', (paste0(dds_run$Genotype, "_", dds_run$Sample_ID)))
name_list

sig_export <- DS_norm_counts %>%
  data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  as_tibble() %>%
  `colnames<-`(name_list) %>%
  semi_join(sig_res) %>%
  dplyr::slice(match(sorted_DEGenes, gene)) #%>%
#filter(!str_detect(gene, unwanted_genes)) 

#Manual Reordering of Columns (if necessary)
reordered_index <- c(
  grep("WT", names(sig_export), ignore.case = T),
  grep("KI", names(sig_export), ignore.case = T))

sig_export <- sig_export %>%
  select(c(1, reordered_index))

# Column Names
check_cols <- c(colnames(sig_export))
check_cols


# Normalized counts
check_counts <- DS_norm_counts %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  dplyr::select(c(1, reordered_index)) %>%
  as_tibble %>% 
  `colnames<-`(check_cols) #%>%
#filter(!str_detect(gene, unwanted_genes))


# Results
check_res <- dds_result %>%
  as.data.frame() %>%
  rownames_to_column('gene') #%>%
#filter(!str_detect(gene, unwanted_genes))

#TPM File
# TPM_name_list <- c((paste0(dds_run$position, "_", dds_run$clinical_2, "_",dds_run$apoe, "_", dds_run$sex)))
# TPM_name_list  
# TPM  <- txicounts %>%
#   as.data.frame() %>%
#   `colnames<-`(TPM_name_list) %>%
#   rownames_to_column('gene')


# Output files
#Statistics
write.xlsx(sig_res, file = paste0('results/', file_prefix, 'DEGene_statistics_pval05.xlsx'), overwrite = T)
write.xlsx(check_res, file = paste0('results/', file_prefix, 'statistics.xlsx'), overwrite = T)

#Counts
write.xlsx(sig_export, file = paste0('results/', file_prefix, 'DEGene_counts_pval05.xlsx'), overwrite = T)
write.xlsx(check_counts, file = paste0('results/', file_prefix, 'DS_counts_check.xlsx'), overwrite = T)
# write.xlsx(TPM,file = paste0('results/', file_prefix, 'TPM.xlsx'), overwrite = T)


#TE results
sig_res_TE <- sig_res %>% filter(grepl(':',gene))
sig_counts_TE <- sig_export %>% filter(grepl(':',gene))


###### ----- BARPLOTS -----#####
goi <- sig_res_TE %>% head(10) %>% pull('gene')
goi_data <-  check_counts %>% filter(gene %in% goi)
goi_data2 <- rename(goi_data, 
                    'WT-' = contains('WT'),
                    'KI-' = contains('KI')
) 

dir.create(paste0('plots/barplots/',file_prefix))


lineWidth = 1
pointSize = 40
for (i in goi) { 
  goi_data3 <- goi_data2 %>% filter(gene == i) %>% pivot_longer(
    cols = -(1:1),
    values_to = c("DS_Counts"),
    names_to = c("condition", "replicate"),
    names_sep  = "-")
  
  goi_data3$condition <- factor(goi_data3$condition, levels=unique(goi_data3$condition))
  
  comparisons <- compare_means(
    data = goi_data3,
    formula = DS_Counts ~ condition,
    method = "t.test",
    p.adjust.method = "BH") %>% filter(p<0.05)
  
  comparison_list_sign <- comparisons %>% mutate(comparison_list = map2(group1, group2,c)) %>% pull(comparison_list)
  
  barplot <- ggplot(goi_data3, aes(x = condition, y = DS_Counts, fill = condition)) + 
    geom_bar(position = 'dodge', stat = 'summary', fun = mean, width = 0.7, colour="black") +
    geom_errorbar(stat = 'summary', position = 'dodge', width = 0.4) +
    geom_point(aes(x = condition), position = position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width=0.9)) +
    expand_limits(x = 0, y = 0) +
    theme(panel.spacing = unit(1, "lines")) +
    labs(x = NULL, y = c("Normalized Counts ±SE")) +
    ggtitle(paste0(i))+
    theme(text = element_text(size = pointSize, colour = "black"),
      rect = element_blank(),
      line = element_line(size = lineWidth, colour = "black"),
      plot.title  = element_text(color="black", size=40, face="bold.italic"),
      axis.title  = element_text(size = pointSize, colour = "black"),
      axis.text.x  = element_text(size = pointSize , colour = "black", angle = 45, hjust = 1),
      axis.text.y  = element_text(size = pointSize , colour = "black"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = pointSize , colour = "black"),
      axis.line = element_line(size = lineWidth, colour = "black"),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
    stat_compare_means(comparisons = comparison_list_sign, method = 't.test', label = "p.signif" ) +
    scale_fill_brewer(palette = "Paired") + scale_y_continuous(expand = c(0, 0, .05, 0))
  
  pdf(file = paste0('plots/barplots/',file_prefix,'/', i,'.pdf'), pointsize = 10, width = 5, height= 10)
  print(barplot)
  dev.off()
}



###### ----- HEATMAP -----#####
# Define a normalization function to calculate Z scores.
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_final <- sig_counts_TE %>% column_to_rownames('gene')
data_subset_50 <- data_final %>% head(n = 50)
data_subset_200 <- data_final %>% head(n = 200)
data_subset_400 <- data_final %>% head(n = 400)

# Normalize your data according to Z score using cal_z_score function.
data_norm <- t(apply(data_final, 1, cal_z_score))
data_subset_200 <- t(apply(data_subset_200, 1, cal_z_score))
data_subset_400 <- t(apply(data_subset_400, 1, cal_z_score))
#data_grouped_norm <- t(apply(data_grouped, 1, cal_z_score))

#Heatmap setup
gaps = c(5)
clusters = 2
dir.create(paste0('plots/heatmaps/',file_prefix))

# Create heatmap
pdf(file = paste0('plots/heatmaps/', file_prefix, '/','ALL_DEGs_pval05.pdf'), pointsize = 10)
phm_full <- pheatmap(data_norm,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 4,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = F,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     #cellwidth = 25,
                     scale = 'row'
                     #cellheight = 3,
)
dev.off()


pdf(file = paste0('plots/heatmaps/', file_prefix, '/','ALL_DEGs_labelled_pval05.pdf'), pointsize = 10, height = 200)
phm_full <- pheatmap(data_norm,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 4,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 10,
                     #cellwidth = 25,
                     scale = 'row',
                     cellheight = 8
)
dev.off()

pdf(file = paste0('plots/heatmaps/', file_prefix, '/','Top_50.pdf'), pointsize = 10, height = 10)
phm_50 <- pheatmap(data_subset_50,
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                   breaks = seq(-2, 2, by = 0.1),
                   kmeans_k = NA,
                   cluster_rows = T,
                   cutree_row = clusters,
                   cluster_cols = F,
                   #cutree_cols = 4,
                   gaps_col = gaps,              
                   legend = TRUE,
                   show_rownames = T,
                   #cellwidth = 25,
                   cellheight = 9,
                   treeheight_col = 0,
                   treeheight_row = 50,
                   border_color = 'NA',
                   fontsize = 8,
                   scale = 'row')
dev.off()

pdf(file = paste0('plots/heatmaps/', file_prefix, '/','Top_200.pdf'), pointsize = 10, height = 35)
phm_200 <- pheatmap(data_subset_200,
                    color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                    breaks = seq(-2, 2, by = 0.1),
                    kmeans_k = NA,
                    cluster_rows = T,
                    cutree_row = clusters,
                    cluster_cols = F,
                    #cutree_cols = 4,
                    gaps_col = gaps,              
                    legend = TRUE,
                    show_rownames = T,
                    #cellwidth = 25,
                    cellheight = 9,
                    treeheight_col = 0,
                    treeheight_row = 50,
                    border_color = 'NA',
                    fontsize = 8,
                    scale = 'row')
dev.off()

pdf(file = paste0('plots/heatmaps/', file_prefix, '/','Top_400.pdf'), pointsize = 10, height = 60)
phm_full <- pheatmap(data_subset_400,
                     color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                     breaks = seq(-2, 2, by = 0.1),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = gaps,             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     #cellwidth = 25,
                     scale = 'row',
                     cellheight = 9)
dev.off()

