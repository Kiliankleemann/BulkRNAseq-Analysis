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
setwd("/home/kilian/sciebo - Kleemann, Kilian (kleemann@uni-bonn.de)@uni-bonn.sciebo.de/Immune_priming_TE/3prime_vs_fulllength/3prime_Samhd1/")

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

#sample_data <- sample_data %>% filter(Genotype_2 == 'WT')

#Import TElocal locus sepcific information 
locus_df <- read.table(paste0('/media/kilian/OS/GTF_files_TEtranscript/mm10_rmsk_TE.gtf.locInd.locations')) %>% row_to_names(row_number = 1)
locus_df <- cSplit(locus_df, 'chromsome:start-stop:strand' , ":")
locus_df <- cSplit(locus_df, 'chromsome:start-stop:strand_2' , "-")

locus_df <- locus_df %>% `colnames<-`(c('TE','Chromosome','Strand','Start','Stop'))
locus_df$Length = locus_df$Stop - locus_df$Start

#Import TElocal counts
sample_files <- list.files(".//TElocal_multi", full.names = T, pattern=NULL, all.files=FALSE)
sample_files
experiment_ctrl <- experiment #%>% filter(Condition_1 == 'ctrl' & Condition_3 == 'no_guide_ctrl')

TElocal_uniq_counts_all<- data.frame(row.names=1:3751066)
for (sample in experiment_ctrl$Sample_ID) {
  TElocal_uniq_counts <- read.table(paste0('TElocal_multi/',sample,'.cntTable')) %>% row_to_names(row_number = 1)
  TElocal_uniq_counts_all <-  cbind(TElocal_uniq_counts_all,TElocal_uniq_counts)
}

TElocal_uniq_counts_all2 <- TElocal_uniq_counts_all[,-1]
rownames(TElocal_uniq_counts_all2) <- TElocal_uniq_counts_all[,1]

TElocal_uniq_counts_all3 <- TElocal_uniq_counts_all2 %>% dplyr::select(!contains('gene/TE'))

TElocal_uniq_counts_all3 <- mutate_all(TElocal_uniq_counts_all3, function(x) as.numeric(as.character(x)))

colClean1 <- function(TEtranscript_multi_counts){colnames(TElocal_uniq_counts_all3) <- gsub("_1Aligned.sortedByCoord.out.bam", "", colnames(TElocal_uniq_counts_all3)); TElocal_uniq_counts_all3 } 
colClean3 <- function(TEtranscript_multi_counts){colnames(TElocal_uniq_counts_all3) <- gsub(".*?multi/", "", colnames(TElocal_uniq_counts_all3)); TElocal_uniq_counts_all3 } 

TElocal_uniq_counts_all3 <- colClean1(TElocal_uniq_counts_all3)
TElocal_uniq_counts_all3 <- colClean3(TElocal_uniq_counts_all3)

#Arrange count_data columns as ordered in sample_data
# TEtranscript_multi_counts_filtered <- TEtranscript_multi_counts %>% select(-contains(outliers))
# TEtranscript_multi_counts_reordered <- TEtranscript_multi_counts_filtered %>% select(paste(sample_data$Sample_ID))


#Analyzis parameters
file_prefix = 'TELocal_multi_WT_vs_Samhd1KO'
experiment = TElocal_uniq_counts_all3
experiment



###### ----- DESEQ2 ANALYSIS -----#####
#Designs 
#dds <- DESeqDataSetFromTximport(txi, colData = experiment, design = ~  Genotype) 
dds <- DESeqDataSetFromMatrix(countData = experiment, colData = sample_data, design = ~ Genotype_2) 

#prefiltering on minimum of 5 reads (NOT REQUIRED AS DESEQ2 OPTIMIZES AND FILTERS AUTOMATICALLY)
dds <- estimateSizeFactors(dds)
# idx <- rowSums(counts(dds, normalized = TRUE)) >= 5
# dds <- dds[idx,]

#Setting the reference level (control group to compare against)
#dds@colData@listData$Condition_1 <- relevel(dds@colData@listData$Condition_1, ref = "ctrl")


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
z <- plotPCA(vst, intgroup=c('Genotype_2'),ntop = 200)

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
  geom_point(aes(color = group)) +
  theme_PCA +
  labs(title = 'PCA (Top 200 variable genes)',   x=paste(z$labels$x), y=paste(z$labels$y))
dev.off()


pdf(file = paste0('plots/PCA/', file_prefix,'_labelled', '.pdf'), pointsize = 10)
ggplot(z$data,aes(x=z$data$PC1, y=z$data$PC2, label = z$data$name)) +
  geom_point(aes(color = Condition_1)) +
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
name_list <- c('gene', (paste0(dds_run$Genotype_2, "_", dds_run$Sample_ID)))
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
# reordered_index <- c(
#   grep("WT", names(sig_export), ignore.case = T),
#   grep("KI", names(sig_export), ignore.case = T))
# 
# sig_export <- sig_export %>%
#   select(c(1, reordered_index))

# Column Names
check_cols <- c(colnames(sig_export))
check_cols


# Normalized counts
check_counts <- DS_norm_counts %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  #dplyr::select(c(1, reordered_index)) %>%
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
dir.create(paste0('results/',file_prefix))
#Statistics
write.xlsx(sig_res, file = paste0('results/', file_prefix, '/DEGene_statistics_pval05.xlsx'), overwrite = T)
write.xlsx(check_res, file = paste0('results/', file_prefix, '/statistics.xlsx'), overwrite = T)

#Counts
write.xlsx(sig_export, file = paste0('results/', file_prefix, '/DEGene_counts_pval05.xlsx'), overwrite = T)
write.xlsx(check_counts, file = paste0('results/', file_prefix, '/DS_counts_check.xlsx'), overwrite = T)
# write.xlsx(TPM,file = paste0('results/', file_prefix, 'TPM.xlsx'), overwrite = T)

#TE results
check_res_TE <- check_res %>% filter(grepl(':',gene))
sig_counts_TE <- sig_export %>% filter(grepl(':',gene))

check_res_TE <- cSplit(check_res_TE,'gene',':')
sig_res_TE <- sig_res_TE[,c(1:7)]

sig_counts_TE <- cSplit(sig_counts_TE,'gene',':')
#sig_counts_TE <- sig_counts_TE %>% select(-c('gene_2','gene_3', 'gene_4'))

#Ctrl samples
TElocal_uniq_counts_all4 <- TElocal_uniq_counts_all3 %>% rownames_to_column('gene') %>% filter(grepl(':',gene)) %>% column_to_rownames('gene')
TElocal_uniq_counts_all4$average <- rowMeans(TElocal_uniq_counts_all4)

TElocal_uniq_counts_all4 <- TElocal_uniq_counts_all4 %>% rownames_to_column('gene') 

TElocal_uniq_counts_all4 <- cSplit(TElocal_uniq_counts_all4,'gene',':')

TElocal_uniq_counts_all4_top <- TElocal_uniq_counts_all4 %>% filter(average > 50)

TElocal_uniq_counts_all4_top <- TElocal_uniq_counts_all4_top %>% rename(TE = gene_1)
TElocal_uniq_counts_all4_top_locus <- left_join(TElocal_uniq_counts_all4_top, locus_df, by = 'TE')
TElocal_mutli_counts_LINE1_top_locus <- TElocal_uniq_counts_all4_top_locus %>% filter(gene_4 == 'LINE' & Length > 5000)

#Merging with Loocus informaiton dataframe
# colnames_TE <- colnames(sig_res_TE)[c(1:6)]
# sig_res_TE <- sig_res_TE %>% `colnames<-`(c( colnames_TE,'TE'))

# check_res_TE <- check_res_TE %>% rename(TE = gene_1)
# check_res_TE_locus <- left_join(check_res_TE, locus_df, by = 'TE')



#Filter sig_res_TE with mean expression > 50
check_res_TE_locus_filter <- check_res_TE_locus %>% filter(gene_4 == 'LINE') %>% filter(baseMean > 10) %>% arrange(desc(baseMean)) %>% filter(Length > 6000)

dir.create(paste0('results/', file_prefix))
write.xlsx(TElocal_mutli_counts_LINE1_top_locus, file = paste0('results/', file_prefix, '/TOP_LINE_TE_filter_3CONTROL.xlsx'), overwrite = T)



###### ----- BARPLOTS -----#####
goi <- sig_res_TE %>% head(10) %>% pull('gene')
goi_data <-  sig_counts_TE %>% filter(gene %in% goi)
goi_data2 <- rename(goi_data, 
                    'ctrl-' = contains('ctrl'),
                    'KI-' = contains('KI')
) 


dir.create(paste0('plots/barplots/', file_prefix))

lineWidth = 1
pointSize = 40
for (i in goi) { 
  goi_data3 <- goi_data2 %>% filter(gene == i) %>% pivot_longer(
    cols = -(1:1),
    values_to = c("DS_Counts"),
    names_to = c("condition", "replicate"),
    names_sep  = "-")
  
  goi_data3$condition <- factor(goi_data3$condition, levels=unique(goi_data3$condition))
  
  pvalue_data <- round(sig_res_TE %>% filter(gene == i) %>% pull(padj),digits = 10)
  
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
    labs(title = paste0(i),
         subtitle = paste0("Padj = ", pvalue_data),
         x = NULL, 
         y = c("Normalized Counts ±SE")) +
    theme(text = element_text(size = pointSize, colour = "black"),
          rect = element_blank(),
          line = element_line(size = lineWidth, colour = "black"),
          plot.title  = element_text(color="black", size=pointSize*0.5 ,hjust= -0.5, face = "italic"),
          plot.subtitle= element_text(size=pointSize *0.5, hjust= -1, face="italic", color="black"),
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
    #stat_compare_means(comparisons = comparison_list_sign, method = 't.test', label = "p.signif" ) +
    scale_fill_brewer(palette = "Paired") + scale_y_continuous(expand = c(0, 0, .05, 0))
  
  pdf(file = paste0('plots/barplots/',file_prefix,'/',i,'.pdf'), pointsize = 10, width = 5, height= 10)
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
data_subset_50 <-  t(apply(data_subset_50, 1, cal_z_score))
data_subset_200 <- t(apply(data_subset_200, 1, cal_z_score))
data_subset_400 <- t(apply(data_subset_400, 1, cal_z_score))
#data_grouped_norm <- t(apply(data_grouped, 1, cal_z_score))

#Heatmap setup
gaps = c(4)
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

##### ----- DONUT PLOTS ----- #####
dir.create(paste0('plots/donut/'))
dir.create(paste0('plots/donut/', file_prefix))

#Piechart showing families 
sig_res_TE <- sig_res %>% filter(grepl(':',gene))
sig_res_TE <- cSplit(sig_res_TE, "gene", ":")
sig_res_TE$UP_or_DOWN[sig_res_TE$log2FoldChange>0] <- 'UP'
sig_res_TE$UP_or_DOWN[sig_res_TE$log2FoldChange<0] <- 'DOWN'

sig_res_TE_UP <- sig_res_TE %>% filter(log2FoldChange > 0)
sig_res_TE_DOWN <- sig_res_TE %>% filter(log2FoldChange < 0)

donut_prefix <- 'sig_res_TE_DOWN'
donut_data <- sig_res_TE_DOWN

hsize <- 2

donut_data <- donut_data %>% 
  group_by(gene_3) %>% # Variable to be transformed
  count() %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc)) %>%
  arrange(desc(gene_3)) %>% ## arrange in the order of the legend
  mutate(text_y = cumsum(n) - n/2) ### calculate where to place the text labels

# Basic Donut of TE elements
piechart <- ggplot(donut_data, aes(x=hsize , y= n, fill=gene_3)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) + 
  geom_label(aes(label = paste(labels), y = text_y),
             colour = "white", fontface = "bold", size = pointSize *0.1, nudge_x = 1,
             show.legend = FALSE) +
  # geom_label_repel(aes(label = labels, y = text_y),
  #                  nudge_x = 0.6, nudge_y = 0.6,
  #                  size = 5, show.legend = F) +
  labs(title = paste('DE TEs Pval05')) +
  guides(fill=guide_legend(title=paste0('Downregulated TEs'))) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        plot.title = element_text(color="black", hjust = -1, size=pointSize, face = "italic"),
        legend.title = element_text(size = pointSize*0.5 , colour = "black"),
        legend.text = element_text(size = pointSize*0.5 , colour = "black"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())+
  #scale_fill_manual(values = c('dodgerblue3','firebrick'))+
  scale_fill_brewer(palette = "Blues") 
#xlim(c(0.2, hsize + 0.5))
piechart

ggsave(paste0('plots/donut/',file_prefix,'/', donut_prefix,'TE_families', '.pdf'),
       plot = piechart,
       device = NULL,
       path = NULL,
       #scale = 1,
       width = 15,
       height = 15,
       units = c("cm"),
       dpi = 600)

#Barplots of overall TE results
bar_prefix <- 'sig_res_TE'
bar_data <- sig_res_TE
bar_data <- bar_data %>% 
  group_by(UP_or_DOWN) %>% # Variable to be transformed
  count() %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc)) %>% 
  mutate(x = hsize)


lineWidth = 1
pointSize = 40
barplot <- ggplot(bar_data, aes(x = UP_or_DOWN, y = n ,fill = UP_or_DOWN)) + 
  geom_bar(position = 'dodge', stat = 'summary', fun = mean, width = 0.7, colour="black") +
  expand_limits(x = 0, y = 0) +
  theme(panel.spacing = unit(1, "lines")) +
  labs(title = paste0('DE TEs Pval05'),
       x = NULL, 
       y = c("Count")) +
  theme(text = element_text(size = pointSize, colour = "black"),
        rect = element_blank(),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(color="black", hjust =1, size=pointSize, face = "italic"),
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
  scale_fill_manual(values = c('dodgerblue3','firebrick'))+
  #scale_fill_brewer(palette = "Paired") + 
  scale_y_continuous(expand = c(0, 0, .05, 0))


pdf(file = paste0('plots/barplots/',file_prefix,'/',bar_prefix,'.pdf'), pointsize = 10, width = 5, height= 10)
print(barplot)
dev.off()



