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
setwd("/Users/kiliankleemann/sciebo - Kleemann, Kilian (kleemann@uni-bonn.de)@uni-bonn.sciebo.de/Immune_priming_TE/Chronic_vs_Acute/Chronic/IFET_vs_HET")
dir.create('results')
dir.create('plots')
dir.create('plots/PCA')
dir.create('plots/barplots')
dir.create('plots/heatmaps')
dir.create('plots/donut/')
dir.create('plots/volcano/')

###### ----- FILE PREPARATION -----#####
# Create a metadata map using your metadata file.
sample_data <- read.xlsx("metadata.xlsx", sheet = 1)
outliers = c('')
sample_data <- sample_data %>%            
  filter(!Sample_ID %in% outliers)

#Import TEtranscript counts
TEtranscript_multi_counts <- read.table(paste0('TETranscripts_multi/TEtranscripts_out.cntTable'),header=T,row.names=1)

colClean1 <- function(TEtranscript_multi_counts){ colnames(TEtranscript_multi_counts) <- gsub("Aligned.sortedByCoord.out.bam.T", "", colnames(TEtranscript_multi_counts)); TEtranscript_multi_counts } 
colClean2 <- function(TEtranscript_multi_counts){ colnames(TEtranscript_multi_counts) <- gsub("Aligned.sortedByCoord.out.bam.C", "", colnames(TEtranscript_multi_counts)); TEtranscript_multi_counts } 
colClean3 <- function(TEtranscript_multi_counts){ colnames(TEtranscript_multi_counts) <- gsub(".*?2023/", "", colnames(TEtranscript_multi_counts)); TEtranscript_multi_counts } 

TEtranscript_multi_counts <- colClean1(TEtranscript_multi_counts)
TEtranscript_multi_counts <- colClean2(TEtranscript_multi_counts)
TEtranscript_multi_counts <- colClean3(TEtranscript_multi_counts)


#Susetting analysis
#sample_data <- sample_data %>% filter(Condition_2 == "f")

#Arrange count_data columns as ordered in sample_data
TEtranscript_multi_counts <- TEtranscript_multi_counts %>% select(-contains(outliers))
TEtranscript_multi_counts_reordered <- TEtranscript_multi_counts %>% select(paste(sample_data$Sample_ID))

#Analyzis parameters
file_prefix = 'TEtranscript_multi_IFET_vs_HET'
reference_compare_to <- "HET"
experiment = TEtranscript_multi_counts_reordered
experiment
sample_data


###### ----- DESEQ2 ANALYSIS -----#####
#Designs 
dds <- DESeqDataSetFromMatrix(countData = experiment, colData = sample_data, design = ~ Genotype) 

#prefiltering on minimum of 5 reads (NOT REQUIRED AS DESEQ2 OPTIMIZES AND FILTERS AUTOMATICALLY)
dds <- estimateSizeFactors(dds)
# idx <- rowSums(counts(dds, normalized = TRUE)) >= 5
# dds <- dds[idx,]

#Setting the reference level (control group to compare against)
dds@colData@listData$Genotype <- relevel(dds@colData@listData$Genotype, ref = reference_compare_to)


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

pdf(file = paste0('plots/PCA/', file_prefix, '.pdf'), pointsize = 10)
ggplot(z$data,aes(x=z$data$PC1, y=z$data$PC2, )) +
  geom_point(aes(color = group)) +
  theme_PCA +
  labs(title = 'PCA (Top 200 variable genes)',   x=paste(z$labels$x), y=paste(z$labels$y))
dev.off()


pdf(file = paste0('plots/PCA/', file_prefix,'_labelled', '.pdf'), pointsize = 10)
ggplot(z$data,aes(x=z$data$PC1, y=z$data$PC2, label = z$data$name)) +
  geom_point(aes(color = group)) +
  theme_PCA + labs(title = 'PCA (Top 200 variable genes)',   x=paste(z$labels$x), y=paste(z$labels$y)) +
  geom_text_repel(aes(label = sample_data$Sample_ID),max.overlaps = Inf )
dev.off()


###### ----- SIG DATA & EXPORT -----#####
# Build significant gene table and extract list of sorted DE genes.
sig_res <- dds_result %>%
  data.frame() %>%
  rownames_to_column(var='gene') %>%  
  as_tibble %>%
  filter(pvalue < 0.05) %>%
  arrange(pvalue) 

sig_res_padj <- sig_res %>% 
  filter(padj < 0.05) %>%
  arrange(pvalue)

#!!!!! Export significant gene count data
name_list <- c('gene', (paste0(dds_run$Genotype, "_", dds_run$Sample_ID)))
name_list

sig_export <- DS_norm_counts %>%
  data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  as_tibble() %>%
  `colnames<-`(name_list) %>%
  semi_join(sig_res) %>%
  dplyr::slice(match(sig_res$gene, gene)) 

sig_export_padj <- sig_export %>% filter(gene %in% sig_res_padj$gene)

#Manual Reordering of Columns (if necessary)
reordered_index <- c(
  grep("IFET", names(sig_export), ignore.case = T),
  grep("HET", names(sig_export), ignore.case = T))

#Ordering
sig_export <- sig_export %>%
  select(c(1, reordered_index))

sig_export_padj <- sig_export_padj %>%
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
  `colnames<-`(check_cols) 

# Results
check_res <- dds_result %>%
  as.data.frame() %>%
  rownames_to_column('gene') 

# Output files
dir.create(paste0('results/', file_prefix))
#Statistics
write.xlsx(sig_res, file = paste0('results/', file_prefix, '/DEGene_statistics_pval05.xlsx'), overwrite = T)
write.xlsx(sig_res_padj, file = paste0('results/', file_prefix, '/DEGene_statistics_padj05.xlsx'), overwrite = T)
write.xlsx(check_res, file = paste0('results/', file_prefix, '/Statistics.xlsx'), overwrite = T)

#Counts
write.xlsx(sig_export, file = paste0('results/', file_prefix, '/DEGene_counts_pval05.xlsx'), overwrite = T)
write.xlsx(sig_export_padj, file = paste0('results/', file_prefix, '/DEGene_counts_padj05.xlsx'), overwrite = T)
write.xlsx(check_counts, file = paste0('results/', file_prefix, '/DS_counts_check.xlsx'), overwrite = T)


#TE results
Pval_cutoff <- 'pval01'

sig_res_TE <- sig_res %>% filter(pvalue<0.01) %>% filter(grepl(':',gene))
sig_counts_TE <- sig_export %>% filter(gene %in% sig_res_TE$gene) %>% filter(grepl(':',gene))


###### ----- BARPLOTS -----#####
goi <- sig_res_TE %>% head(40) %>% pull('gene')
goi_data <-  check_counts %>% filter(gene %in% goi)
goi_data2 <- rename(goi_data, 
                    'IFET-' = contains('IFET'),
                    'HET-' = contains('HET')
) 

dir.create(paste0('plots/barplots/',file_prefix,"_",Pval_cutoff))
barplot_dir <- paste0('plots/barplots/',file_prefix ,"_",Pval_cutoff,"/")


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
  
  pdf(file = paste0(barplot_dir, i,'.pdf'), pointsize = 10, width = 5, height= 10)
  print(barplot)
  dev.off()
}



###### ----- HEATMAPS -----#####
data_final <- sig_counts_TE %>% column_to_rownames('gene')
data_subset_50 <- data_final %>% head(n = 50)
data_subset_200 <- data_final %>% head(n = 200)
data_subset_400 <- data_final %>% head(n = 400)

#Heatmap setup
gaps = c(9)
clusters = 2
dir.create(paste0('plots/heatmaps/',file_prefix ,"_",Pval_cutoff))
heatmap_dir <- paste0('plots/heatmaps/',file_prefix ,"_",Pval_cutoff,"/")

# Create heatmap
phm_full <- pheatmap(data_final,
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
                     scale = 'row')

ggsave(paste0(heatmap_dir,'ALL_DETEs.pdf'),
       phm_full,
       dpi = 300)


phm_full <- pheatmap(data_final,
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
ggsave(paste0(heatmap_dir,'ALL_DETEs_labelled.pdf'),
       phm_full,
       height = 20,
       dpi = 300)


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
ggsave(paste0(heatmap_dir,'Top50_labelled.pdf'),
       phm_50,
       height = 15,
       dpi = 300)

#Export Heatmap per Family
data_final_family <- cSplit(sig_counts_TE, "gene", ":")
families_unique <- data_final_family$gene_3 %>% unique()

clusters <- "1"


for (TE in families_unique) {
  data_final_family_select <- data_final_family %>% 
    filter(gene_3 %in% TE) 
  data_final_family_select$gene <- paste0(data_final_family_select$gene_1,":", data_final_family_select$gene_2,":", data_final_family_select$gene_3)
  data_final_family_select <- data_final_family_select %>% 
    select(!c("gene_1","gene_2", "gene_3")) %>%
    column_to_rownames('gene')
  height_heatmap <- as.numeric(paste0(nrow(data_final_family_select))) *2
  phm_Fam <- pheatmap(data_final_family_select,
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
         units = 'in',
         dpi = 300)
}


###### ----- VOLCANO PLOT ----- ######
data_final_family <- cSplit(sig_res_TE, "gene", ":")

dir.create(paste0('plots/volcano/',file_prefix,"_",Pval_cutoff))
volcano_dir <- paste0('plots/volcano/',file_prefix,"_",Pval_cutoff,"/")

for (TE in families_unique) {
  data_final_family_select <- data_final_family %>% 
    filter(gene_3 %in% TE) 
  data_final_family_select$gene <- paste0(data_final_family_select$gene_1,":", data_final_family_select$gene_2,":", data_final_family_select$gene_3)
  data_final_family_select <- data_final_family_select %>% 
    select(!c("gene_1","gene_2", "gene_3")) 
  
  data_final_family_select$log_pval <- -log10(data_final_family_select$pvalue)
  data_final_family_select$log_baseMean <- log10(data_final_family_select$baseMean)
  col_max <- max(data_final_family_select$log2FoldChange)
  yinterct <- min(data_final_family_select$log_pval)
  
  scatter_volcano <- ggplot(data_final_family_select, aes(x= log2FoldChange, y= log_pval, label = gene)) +
  geom_point(data = data_final_family_select , aes(fill = log2FoldChange),size = 3, shape = 21, color = 'black') +
    scale_fill_gradient2(limits=c(-col_max, col_max), low="navyblue", mid="whitesmoke", high = "firebrick", na.value = 'navyblue') +
    theme(aspect.ratio = 1, 
        panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title  = element_text(size = 20, colour = "black"),
        axis.text.x=element_text(colour="black",size =20),
        axis.text.y=element_text(colour="black",size =20),
        axis.ticks=element_line(colour="black"),
        plot.margin=unit(c(1,1,1,1),"line"),
        plot.title = element_text(colour="black",size =20)) +
    geom_vline(xintercept = 0,colour="grey", linetype = "longdash") +
    geom_hline(yintercept = yinterct, colour="grey", linetype = "longdash") +
    labs(title = str_wrap(paste0(TE),60), x = "Log2FC", y = "-log(pvalue)") +
  geom_text_repel(data = data_final_family_select,
                  color = 'black', size = 2, fontface = 'italic',
                  min.segment.length = 0, max.overlaps = Inf)
  ggsave(paste0(volcano_dir,TE, '.pdf'),
       plot = scatter_volcano,
       device = NULL,
       path = NULL,
       #scale = 1,
       width = 20,
       height = 20,
       units = c("cm"),
       dpi = 600)
}

#All TEs
sig_res_TE$log_pval <- -log10(sig_res_TE$pvalue)
sig_res_TE$log_baseMean <- log10(sig_res_TE$baseMean)
#Color
col_max <- max(sig_res_TE$log2FoldChange)
yinterct <- min(sig_res_TE$log_pval)
  
scatter_volcano <- ggplot(sig_res_TE, aes(x= log2FoldChange, y= log_pval, label = gene)) +
  geom_point(data = sig_res_TE , aes(fill = log2FoldChange),size = 3, shape = 21, color = 'black') +
  scale_fill_gradient2(limits=c(-col_max, col_max), low="navyblue", mid="whitesmoke", high = "firebrick", na.value = 'navyblue') +
  theme(aspect.ratio = 1, 
        panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title  = element_text(size = 20, colour = "black"),
        axis.text.x=element_text(colour="black",size =20),
        axis.text.y=element_text(colour="black",size =20),
        axis.ticks=element_line(colour="black"),
        plot.margin=unit(c(1,1,1,1),"line"),
        plot.title = element_text(colour="black",size =20)) +
  geom_vline(xintercept = 0,colour="grey", linetype = "longdash") +
  geom_hline(yintercept = yinterct, colour="grey", linetype = "longdash") +
  labs(title = str_wrap(paste0("ALL_DETEs"),60), x = "Log2FC", y = "-log(pvalue)") 
ggsave(paste0(volcano_dir,"ALL_DETEs", '.pdf'),
       plot = scatter_volcano,
       device = NULL,
       path = NULL,
       #scale = 1,
       width = 20,
       height = 20,
       units = c("cm"),
       dpi = 600)


##### ----- DONUT PLOTS ----- #####
dir.create(paste0('plots/donut/',file_prefix,"_",Pval_cutoff))
donut_dir <- paste0('plots/donut/',file_prefix,"_",Pval_cutoff,"/")

sig_res_TE_dnt <- cSplit(sig_res_TE, "gene", ":")
sig_res_TE_dnt$UP_or_DOWN[sig_res_TE_dnt$log2FoldChange>0] <- 'UP'
sig_res_TE_dnt$UP_or_DOWN[sig_res_TE_dnt$log2FoldChange<0] <- 'DOWN'

sig_res_TE_UP <- sig_res_TE_dnt %>% filter(log2FoldChange > 0)
sig_res_TE_DOWN <- sig_res_TE_dnt %>% filter(log2FoldChange < 0)

# Piechart showing families UP 
donut_prefix <- 'sig_res_TE_UP'
donut_data <- sig_res_TE_UP

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
  labs(title = paste(paste(donut_prefix))) +
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

ggsave(paste0(donut_dir, donut_prefix,'_families', '.pdf'),
       plot = piechart,
       device = NULL,
       path = NULL,
       #scale = 1,
       width = 15,
       height = 15,
       units = c("cm"),
       dpi = 600)

#Piechart showing families DOWN 
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
  labs(title = paste(paste(donut_prefix))) +
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

ggsave(paste0(donut_dir, donut_prefix,'_families', '.pdf'),
       plot = piechart,
       device = NULL,
       path = NULL,
       #scale = 1,
       width = 15,
       height = 15,
       units = c("cm"),
       dpi = 600)

#Barplots of overall TE results
bar_prefix <- 'All_DETEs'
bar_data <- sig_res_TE_dnt
bar_data <- bar_data %>% 
  group_by(UP_or_DOWN) %>% 
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
  labs(title = paste0(Pval_cutoff),
       x = NULL, 
       y = c("Count")) +
  theme(text = element_text(size = pointSize, colour = "black"),
        rect = element_blank(),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(color="black", hjust =1, size=pointSize *0.7, face =),
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
  scale_fill_manual(values = c('blue3','red3'))+
  #scale_fill_brewer(palette = "Paired") + 
  scale_y_continuous(expand = c(0, 0, .05, 0))


ggsave(paste0(donut_dir, bar_prefix,'_UP_and_DOWN_families', '.pdf'),
       plot = barplot,
       device = NULL,
       path = NULL,
       #scale = 1,
       width = 15,
       height = 15,
       units = c("cm"),
       dpi = 600)

dev.off()

