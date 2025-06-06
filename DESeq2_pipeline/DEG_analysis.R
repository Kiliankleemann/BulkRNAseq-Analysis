#### -------- INSTALLING AND LOADING PACKAGES ------- ####
list.of.packages <- c("BiocGenerics","tximport","S4Vectors", "DESeq2", "biomaRt",'GseaVis',
                      "ggplot2", "ggsignif", "ggpubr", "sva", "devtools", "org.Hs.eg.db", 
                      "org.Mm.eg.db", "limma","stringr","KEGGREST","ggrepel", "openxlsx",
                      "fgsea","clusterProfiler","pheatmap","ggpubr","cowplot",'dplyr',
                      "RColorBrewer",'AnnotationDbi', 'tidyverse','pheatmap', 'dendextend')


new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
if(length(new.packages)) BiocManager::install(new.packages)

# Special packages
#devtools::install_github("junjunlab/GseaVis")

# Packages loading
invisible(lapply(list.of.packages, library, character.only = TRUE))

#### -------- SET WORKING DIRECTORY ------- ####
setwd("~/directory")


#### -------- CREATING DIRECTORIES ------- ####
dir.create('results')
dir.create('plots')
dir.create('plots/PCA')
dir.create('plots/heatmaps')
dir.create('plots/barplots')


#### -------- FILE PREPARATION ------- ####
# useMart is a function used to connect to the selected BioMart database and dataset (ListEnsemblArchives() to list host sites)
ensembl_mm_mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
# getBM is a function used to retrieve information from the BioMart database
ensembl_df <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id_version", "external_gene_name"), mart = ensembl_mm_mart) %>%
  as_tibble() %>%
  mutate(noVersion = sub('\\..*', '', ensembl_transcript_id_version))

# Build tx2gene table from ensembl_df
tx2gene <- ensembl_df %>%
  dplyr::select(noVersion, external_gene_name) %>%
  rename(c('Transcript' = "noVersion" , "Gene" = "external_gene_name" )) %>%
  data.frame()

# Create a metadata map using your metadata file.
sample_data <- read.xlsx("metadata.xlsx", sheet = 1)

# Outlier filtering  and subsetting(if necessary)
outliers = c('')
sample_data <- sample_data %>%            # Filter sample_data sheet
  filter(!Sample_ID %in% outliers)

# xxx <- sample_data %>% filter(condition_2 == 'xxx')
# xxx <- sample_data %>% filter(condition_2 == 'xxx')

#Setting experiment conditions 
file_prefix <- 'All_samples'
experiment <- sample_data
experiment


## List all directories containing data
#*** All the gene quantification files are in the folder "quant_files", make sure metafile sample order matches order of quant files.
all_files <- list.files(".//transcript_quant", full.names = T, pattern=NULL, all.files=FALSE)
quant_files <- file.path(all_files, "quant.sf")

position_list <- paste0(experiment$Sample_ID, collapse = '|')
position_list
sample_files <- grep(position_list, quant_files, value = TRUE)
sample_files

## QC check - Order of metafile samples must match quant file order!  Fix metadata file if not!
all(file.exists(quant_files))
sample_files

# Import Transcript quantification files
txi<-tximport(files=sample_files, type="salmon", tx2gene = tx2gene, ignoreTxVersion = T, importer=read.delim)

## Look at the counts and set the counts as a matrix file.
txicounts<-as.matrix(txi$counts)

# Write the counts to file, round them up, and convert to data.frame and remove txicounts file.
count_data <- txi$counts %>%
  round() %>%
  data.frame()

#### -------- ANALYSIS ------- ####
#Designs 
dds <- DESeqDataSetFromTximport(txi, colData = experiment, design = ~  Condition) 


#prefiltering on minimum of 5 reads (NOT REQUIRED AS DESEQ2 OPTIMIZES AND FILTERS AUTOMATICALLY)
dds <- estimateSizeFactors(dds)
# idx <- rowSums(counts(dds, normalized = TRUE)) >= 5
# dds <- dds[idx,]

#Setting the reference level (control group to compare against)
dds@colData@listData$condition_1 <- relevel(dds@colData@listData$condition_1, ref = "Ctrl")


## Run DESeq analysis to gather differential expression results
#Run DESeq (LRT)
dds_run <- DESeq(dds, test = 'LRT', reduced = ~ 1)
# Run DESeq (Wald)
#dds_run <- DESeq(dds, betaPrior = F)

# View names of estimated effects
results_names <- resultsNames(dds_run)

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
# plotDispEsts(dds_run)
# plotCounts(dds_run, 'Fus', intgroup = 'condition_1')
# 


####### ------ PCA ------- ####
# Transform counts for data visualization
vst <- vst(dds_run, blind=TRUE)

# Add nametags
z <- plotPCA(vst, intgroup=c("Disease" , "Genotype" ,"Sex"),ntop = 200)

#Adding extra group
z$data$group_2 <- paste0(z$data$Disease,z$data$Genotype)
#Calculating plotting ratio
process_PC <- function(PC) {
  PC <- gsub("^.*?\\s", "", PC)
  PC <- gsub("% variance", "", PC)
  PC <- as.numeric(PC)
  return(PC)
}
PC1 <- process_PC(z$labels$x)
PC2 <- process_PC(z$labels$y)
PCA_plot_ratio <- PC2/PC1

#Theme for PCA
theme_PCA <- theme(aspect.ratio = PCA_plot_ratio, 
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

# Assuming you have a DESeqTransform object called 'rld' (from rlog or vst transformation)
pca_data <- plotPCA(vst, intgroup = "Group", returnData = TRUE)

# Extract normalized data from the rlog or vst object
norm_data <- assay(vst)
# Remove genes with zero variance across samples
norm_data_filtered <- norm_data[rowVars(norm_data) > 0, ]

# Perform PCA
pca_result <- prcomp(t(norm_data_filtered), scale. = TRUE)

# Extract loadings
loadings <- pca_result$rotation

# Get top contributing genes for PC1 and PC2
top_genes_PC1 <- names(sort(abs(loadings[, 1]), decreasing = TRUE)[1:10])
top_genes_PC2 <- names(sort(abs(loadings[, 2]), decreasing = TRUE)[1:10])

# Combine top genes
top_genes <- unique(c(top_genes_PC1, top_genes_PC2))

# Prepare data for ggplot
pca_data <- as.data.frame(pca_result$x)
pca_data$sample <- rownames(pca_data)

# Add labels for top contributing genes
gene_labels <- as.data.frame(loadings[top_genes, c(1, 2)])
gene_labels$gene <- rownames(gene_labels)

# Subset PCA results for the top genes
loadings_top <- loadings[top_genes, ]

#Plot no labels
dir.create(paste0('plots/PCA/',file_prefix))
pdf(file = paste0('plots/PCA/', file_prefix, '/PCA_top200.pdf'), pointsize = 10)
ggplot(z$data,aes(x=z$data$PC1, y=z$data$PC2, )) +
  geom_point(aes(color = group_2, shape =Sex)) +
  theme_PCA +
  labs(title = 'PCA (Top 200 variable genes)',   x=paste(z$labels$x), y=paste(z$labels$y)) +
  stat_ellipse(aes(fill=group_2), geom = 'polygon',alpha = 0.1,level = 0.9)
dev.off()

# #second version
# pdf(file = paste0('plots/PCA/', file_prefix, '/PCA2_top200.pdf'), pointsize = 10)
# ggplot(z$data,aes(x=z$data$PC1, y=z$data$PC2, )) +
#   geom_point(aes(color = group_2, shape =Sex)) +
#   theme_PCA +
#   labs(title = 'PCA (Top 200 variable genes)',   x=paste(z$labels$x), y=paste(z$labels$y)) +
#   stat_ellipse(aes(fill=group_2), geom = 'polygon',alpha = 0.1,level = 0.6)
# dev.off()

#Plot with sample labels
pdf(file = paste0('plots/PCA/', file_prefix,'/PCA_top200_labelled', '.pdf'), pointsize = 10)
ggplot(z$data,aes(x=z$data$PC1, y=z$data$PC2, label = z$data$name)) +
  geom_point(aes(color = group_2, shape =Sex)) +
  theme_PCA + labs(title = 'PCA (Top 200 variable genes)',   x=paste(z$labels$x), y=paste(z$labels$y)) +
  geom_text_repel(aes(label = experiment$Sample_ID),max.overlaps = Inf, size = 1) +
  stat_ellipse(aes(fill=group_2), geom = 'polygon',alpha = 0.1,level = 0.9)
dev.off()

#### -------- OUTPUT TABLES & RESULTS ------- ####
#unwanted_genes = paste(c('^Gm', '^mt-', '^Vmn', '^Rpl', '^Rps', '^Olfr','Rik'), collapse = '|')
# Build significant gene table and extract list of sorted DE genes.
sig_res <- dds_result %>%
  data.frame() %>%
  rownames_to_column(var='gene') %>%  
  as_tibble %>%
  filter(pvalue < 0.05) %>%
  arrange(pvalue) #%>% 
#filter(!str_detect(gene, unwanted_genes))

sig_res_adj <- sig_res %>%
  filter(padj < 0.05) %>%
  arrange(padj) 

sorted_DEGenes <- sig_res$gene
sorted_DEGenes_padj <- sig_res_adj$gene

# Export significant gene count data
name_list <- c('gene', (paste0(dds_run$condition_1, "_",dds_run$condition_2,"_",dds_run$Sample_ID)))
name_list

sig_export <- DS_norm_counts %>%
  data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  as_tibble() %>%
  `colnames<-`(name_list) %>%
  semi_join(sig_res) %>%
  dplyr::slice(match(sorted_DEGenes, gene)) #%>%
#filter(!str_detect(gene, unwanted_genes)) 

sig_export_adj <- DS_norm_counts %>%
  data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  as_tibble() %>%
  `colnames<-`(name_list) %>%
  semi_join(sig_res) %>%
  dplyr::slice(match(sorted_DEGenes_padj, gene))


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
dir.create(paste0('results/', file_prefix))
#Statistics
write.xlsx(sig_res, file = paste0('results/', file_prefix, '/DEGene_statistics_pval05.xlsx'), overwrite = T)
write.xlsx(sig_res_adj, file = paste0('results/', file_prefix, '/DEGene_statistics_padj05.xlsx'), overwrite = T)
write.xlsx(check_res, file = paste0('results/', file_prefix, '/Statistics.xlsx'), overwrite = T)
#Counts
write.xlsx(sig_export, file = paste0('results/', file_prefix, '/DEGene_counts_pval05.xlsx'), overwrite = T)
write.xlsx(sig_export_adj, file = paste0('results/', file_prefix, '/DEGene_counts_padj05.xlsx'), overwrite = T)
write.xlsx(check_counts, file = paste0('results/', file_prefix, '/DS_counts.xlsx'), overwrite = T)


##### ------BARPLOTS ------- #######
goi <- sig_res %>% head(10) %>% pull('gene')
goi_data <-  check_counts %>% filter(gene %in% goi)
goi_data2 <- rename(goi_data, 
                    'control_PBS-' = contains('control_PBS'),
                    'selisistat_PBS-' = contains('selisistat_PBS'),
                    'resveratrol_PBS-' = contains('resveratrol_PBS'),
                    'control_LPS-' = contains('control_LPS'),
                    'selisistat_LPS-' = contains('selisistat_LPS'),
                    'resveratrol_LPS-' = contains('resveratrol_LPS')) 

#Settings for barplots
lineWidth = 1
pointSize = 25
dir.create(paste0('plots/barplots/',file_prefix))
#

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
    theme(
      text = element_text(size = pointSize, colour = "black"),
      rect = element_blank(),
      line = element_line(linewidth = lineWidth, colour = "black"),
      plot.title  = element_text(color="black", size=40, face="bold.italic"),
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
  
  pdf(file = paste0('plots/barplots/',file_prefix,'/',  i,'.pdf'), pointsize = 10, width = 5, height= 10)
  print(barplot)
  dev.off()
}





