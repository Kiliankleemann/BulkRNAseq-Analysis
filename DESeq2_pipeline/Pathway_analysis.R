#Setting up directories for analysis
dir.create(paste('Pathway_analysis/'))
dir.create(paste('Pathway_analysis/FGSEA_enrichment/'))
dir.create(paste('Pathway_analysis/KEGG_enrichment/'))
dir.create(paste('Pathway_analysis/GO_enrichment/'))

#Setting up organism data base for enrichment analysis
organism_name <-  'mmu' #human: hsa
OrgDb_name    <-  "org.Mm.eg.db" #human: 'org.Hs.eg.db'


##### --------- FGSEA analysis ####
pathway_data <- read.xlsx(paste0('results/',file_prefix, '/DEGene_statistics_padj05.xlsx')) %>% filter(!grepl(':',gene))

DEgenes <- pathway_data %>% filter(padj < 0.05)  %>% pull('gene')
#DEgenes <- mapIds(org.Mm.eg.db, DEgenes, 'ENTREZID', 'SYMBOL')

DEG_Log2FC <-  pathway_data %>% filter(padj < 0.05) %>%  pull('log2FoldChange')

names(DEG_Log2FC) = DEgenes
DEG_Log2FC <- sort(DEG_Log2FC, decreasing = TRUE)
DEG_Log2FC

#####CURATED pathways
pathways.all.curated <- gmtPathways("/media/kilian/OS/References/GSEA_database/m2.all.v2023.2.Mm.symbols.gmt") 
fgseaRes <- fgsea(pathways=pathways.all.curated, stats=DEG_Log2FC)
fgseaResTidy_curated <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  filter(pval<0.05)

dir.create(paste0('Pathway_analysis/FGSEA_enrichment/',file_prefix))
write.xlsx(fgseaResTidy_curated, paste0('Pathway_analysis/FGSEA_enrichment/',file_prefix,'/Curated_pathways.xlsx'))


#Select based on name or function 
# toMatch <- c("FOXP3", "IL17", "ANTIGEN")
# selected_data <- dplyr::filter(fgseaResTidy_curated,pathway %in% unique(grep(paste(toMatch,collapse="|"), fgseaResTidy_curated$pathway, value=TRUE))) 
# selected_data$log_pval <- -log10(selected_data$pval)

fgseaResTidy_curated$log_pval <- -log10(fgseaResTidy_curated$pval)
FGSEA_data <- fgseaResTidy_curated %>% head(15)

#Up
col_max <- 'firebrick1'
col_min <- 'navyblue'

#Fgsea plot theme
lineWidth = 0.5
pointSize = 20
theme_gsea_bar <- theme(
  text = element_text(size = pointSize, colour = "black"),
  rect = element_blank(),
  line = element_line(size = lineWidth, colour = "black"),
  plot.title  = element_text(color="black", size = pointSize),
  axis.title  = element_text(size = pointSize, colour = "black"),
  axis.text.x  = element_text(size = pointSize, colour = "black"),
  axis.text.y  = element_text(size = pointSize, colour = "black"),
  axis.ticks.x = element_line(size = lineWidth, colour = "black"),
  axis.ticks.y = element_line(size = lineWidth, colour = "black"),
  axis.line = element_line(size = lineWidth, colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.title = element_text(size = pointSize , colour = "black"),
  legend.text = element_text(size = pointSize, colour = "black"),
  legend.key.height = unit(0.5, "cm"),
  legend.key.width = unit(0.5, "cm"))

#Barplot
barplot_curated <- ggplot(FGSEA_data, aes(NES,reorder(pathway, NES))) +
  geom_col(aes(fill= NES), width = 0.5,color = 'black') +
  geom_point(aes(size = log_pval))+
  theme_gsea_bar +
  scale_fill_gradient2(low = col_min, high = col_max) +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title= str_wrap("Curated pathways NES from GSEA",width =30)) + 
  geom_vline(xintercept = 0, color ='black') +
  scale_x_continuous(expand = c(0, 0, 0.1, 0))

barplot_curated

dev.off()

dir.create(paste0('plots/pathways/'))
dir.create(paste0('plots/pathways/FGSEA_pathways/'))
dir.create(paste0('plots/pathways/FGSEA_pathways/',file_prefix))
ggsave(paste0('plots/pathways/FGSEA_pathways/', file_prefix,'/Barplot_curated_pathways.pdf'),
       barplot_curated,
       height = 20,
       width = 30,
       dpi = 300)


#### ---- KEGG ENRICHMENT ---- #####
pathway_data <- read.xlsx(paste0('results/',file_prefix, '/DEGene_statistics_pval05.xlsx'))

DEgenes <- pathway_data %>% filter(pvalue < 0.01)  %>% pull('gene')
DEgenes <- mapIds(org.Mm.eg.db, DEgenes, 'ENTREZID', 'SYMBOL')

DEG_Log2FC <-  pathway_data %>% filter(pvalue < 0.01) %>%  pull('log2FoldChange')

names(DEG_Log2FC) = DEgenes
DEG_Log2FC <- sort(DEG_Log2FC, decreasing = TRUE)
DEG_Log2FC

#Run KEGG 
kegg_result <- gseKEGG(geneList = DEG_Log2FC,
                       organism     = organism_name,
                       pAdjustMethod = "BH",
                       minGSSize = 15,
                       maxGSSize = 500,
                       nPermSimple = 1000,
                       pvalueCutoff = 0.05)

kegg_result <- setReadable(kegg_result,OrgDb= 'org.Mm.eg.db', keyType = 'ENTREZID')
#head(kegg_result@geneSets)

#GSEA plot
# gseaNb(object = kk2,
#        geneSetID = 'hsa05169',
#        subPlot = 2,
#        addGene = T,
#        kegg = T,
#        markTopgene = T)
# 
# ggsave(paste0('plots/pathways/KEGG_pathways/', file_export,'.pdf'),
#        gsea_KEGG,
#        dpi = 300)


#Barplot
KEGG_data <- as.data.frame(kegg_result) %>% filter(NES > 0.5 | NES < -1.5)
KEGG_data$log_pval <- log(KEGG_data$pvalue,10)

#Filter names 
KEGG_data$Description <- gsub('Mus','', KEGG_data$Description)
KEGG_data$Description <- gsub('musculus','', KEGG_data$Description)
KEGG_data$Description <- gsub('house','', KEGG_data$Description)
KEGG_data$Description <- gsub('mouse','', KEGG_data$Description)
KEGG_data$Description <- gsub('\\( )','', KEGG_data$Description)
KEGG_data$Description <- gsub('\\-','', KEGG_data$Description)
KEGG_data$Description <- gsub('    ','', KEGG_data$Description)

barplot_Kegg <- ggplot(KEGG_data, aes(NES,reorder(Description, NES))) +
  geom_col(aes(fill= NES), width = 0.5,color = 'black') +
  geom_point(aes(size = log_pval))+
  theme_gsea_bar +
  scale_fill_gradient2(low = col_min, high = col_max) +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title= str_wrap("Curated pathways NES from GSEA",width =30)) + 
  geom_vline(xintercept = 0, color ='black') +
  scale_x_continuous(expand = c(0, 0, 0.1, 0))

barplot_Kegg

dev.off()

dir.create(paste0('plots/pathways/'))
dir.create(paste0('plots/pathways/KEGG_pathways/'))
dir.create(paste0('plots/pathways/KEGG_pathways/',file_prefix))
ggsave(paste0('plots/pathways/KEGG_pathways/', file_prefix,'/Barplot_pathways.pdf'),
       barplot_Kegg,
       height = 20,
       width = 30,
       dpi = 300)
dir.create(paste0('Pathway_analysis/KEGG_enrichment/',file_prefix))
write.xlsx(KEGG_data, paste0('Pathway_analysis/KEGG_enrichment/',file_prefix,'/KEGG_pathways.xlsx'))


#### ---- EXTRACTING SPECIFIC PATHWAYS - HEATMAP ---- ####
toMatch <- c("Antigen")
pathway_data <- dplyr::filter(kegg_result,Description %in% unique(grep(paste(toMatch,collapse="|"), kegg_result$Description, value=TRUE))) 
pathway_ID <- pathway_data@result$ID


#Extract genes and make heatmap for pathways matching term
cluster_no = 2
gaps = c(4)
data <- read.xlsx(paste0('results/',file_prefix, '/DEGene_counts_pval05.xlsx'))

remove_species_info <- function(input_strings) {
  cleaned_strings <- gsub('Mus','', input_strings)
  cleaned_strings <- gsub('musculus','', cleaned_strings)
  cleaned_strings <- gsub('house','', cleaned_strings)
  cleaned_strings <- gsub('mouse','', cleaned_strings)
  cleaned_strings <- gsub('Mus','', cleaned_strings)
  cleaned_strings <- gsub('\\( )','', cleaned_strings)
  cleaned_strings <- gsub('\\-','', cleaned_strings)
  cleaned_strings <- gsub('    ','', cleaned_strings)
  return(cleaned_strings)
}

for (id in pathway_ID) {
  pathway_data_genes <- pathway_data[[id]]
  pathway_name <- pathway_data[id] %>% pull(Description)
  pathway_name <- remove_species_info(pathway_name)
  pathway_data_selected <- data %>% filter(gene %in% pathway_data_genes) %>% column_to_rownames('gene')
  phm_kegg <- pheatmap(pathway_data_selected,
                       color = colorRampPalette(color_set)(40),
                       breaks = seq(-2, 2, by = 0.1),
                       kmeans_k = NA,
                       cluster_rows = F,
                       cutree_row = cluster_no,
                       cluster_cols = F,
                       gaps_col = gaps,             
                       legend = TRUE,
                       show_rownames = T,
                       border_color = 'NA',
                       fontsize = 10,
                       treeheight_col = 50,
                       treeheight_row = 50,
                       cellwidth = 10,
                       scale = 'row',
                       cellheight = 10,
                       main = pathway_name
  )
  
  ggsave(paste0('plots/pathways/KEGG_pathways/',file_prefix, '/Heatmap_',pathway_name , '_', 'pval05.pdf'),
         phm_kegg,
         height = 10,
         width = 5,
         dpi = 300)
  
}





#### ---- GO ENRICHMENT ---- #####
pathway_data <- read.xlsx(paste0('results/',file_prefix, '/DEGene_statistics_pval05.xlsx'))

DEgenes <- pathway_data %>% filter(pvalue < 0.05)  %>% pull('gene')
DEgenes <- mapIds(org.Mm.eg.db, DEgenes, 'ENTREZID', 'SYMBOL')

DEG_Log2FC <-  pathway_data %>% filter(pvalue < 0.05) %>%  pull('log2FoldChange')

names(DEG_Log2FC) = DEgenes
DEG_Log2FC <- sort(DEG_Log2FC, decreasing = TRUE)
DEG_Log2FC

#Run GO 
go_result <- gseGO(geneList = DEG_Log2FC,
                       OrgDb     = OrgDb_name,
                       pAdjustMethod = "BH",
                       minGSSize = 15,
                       maxGSSize = 500,
                       nPermSimple = 1000,
                       pvalueCutoff = 0.05)

go_result <- setReadable(go_result, OrgDb= OrgDb_name, keyType = 'ENTREZID')


#Barplot
GO_data <- as.data.frame(go_result) %>% filter(NES > 1.8 | NES < -2.4)
GO_data$log_pval <- log(GO_data$pvalue,10)

barplot_Go <- ggplot(GO_data, aes(NES,reorder(Description, NES))) +
  geom_col(aes(fill= NES), width = 0.5,color = 'black') +
  geom_point(aes(size = log_pval))+
  theme_gsea_bar +
  scale_fill_gradient2(low = col_min, high = col_max) +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title= str_wrap("Curated pathways NES from GSEA",width =30)) + 
  geom_vline(xintercept = 0, color ='black') +
  scale_x_continuous(expand = c(0, 0, 0.1, 0))

barplot_Go

dev.off()
dir.create(paste0('plots/pathways/GO_pathways/'))
dir.create(paste0('plots/pathways/GO_pathways/',file_prefix))
ggsave(paste0('plots/pathways/GO_pathways/', file_prefix,'/Barplot_pathways.pdf'),
       barplot_Go,
       height = 20,
       width = 30,
       dpi = 300)

dir.create(paste0('Pathway_analysis/GO_enrichment/',file_prefix))
write.xlsx(GO_data, paste0('Pathway_analysis/GO_enrichment/',file_prefix,'/GO_pathways.xlsx'))

#### ---- EXTRACTING SPECIFIC PATHWAYS  ---- ####
toMatch <- c("immune")
go_result
pathway_data <- dplyr::filter(go_result,Description %in% unique(grep(paste(toMatch,collapse="|"), go_result$Description, value=TRUE))) 
pathway_ID <- pathway_data@result$ID
pathway_ID

cluster_no = 2
gaps = c(4)
data <- read.xlsx(paste0('results/',file_prefix, '/DEGene_counts_pval05.xlsx'))


for (id in pathway_ID) {
  pathway_data_genes <- pathway_data[[id]]
  pathway_name <- pathway_data[id] %>% pull(Description)
  pathway_data_selected <- data %>% filter(gene %in% pathway_data_genes) %>% column_to_rownames('gene')
  phm_kegg <- pheatmap(pathway_data_selected,
                       color = colorRampPalette(color_set)(40),
                       breaks = seq(-2, 2, by = 0.1),
                       kmeans_k = NA,
                       cluster_rows = F,
                       cutree_row = cluster_no,
                       cluster_cols = F,
                       gaps_col = gaps,             
                       legend = TRUE,
                       show_rownames = T,
                       border_color = 'NA',
                       fontsize = 10,
                       treeheight_col = 50,
                       treeheight_row = 50,
                       cellwidth = 10,
                       scale = 'row',
                       cellheight = 10,
                       main = pathway_name
  )
  
  ggsave(paste0('plots/pathways/GO_pathways/',file_prefix, '/Heatmap_',pathway_name , '_', 'pval05.pdf'),
         phm_kegg,
         height = 30,
         units = 'cm',
         width = 12,
         dpi = 300)
  
}







#### ---- HEATMAP CLUSTER ENRICHMENT ANALYSIS ---- ####
data <- read.xlsx(paste0('results/', file_prefix, 'clustergenes.xlsx'))
clusters <- unique(data$cluster)
for (i in clusters) {
  cluster_genes <- data %>% filter(cluster ==i) %>% pull(gene)
  cluster_genes <- mapIds(org.Mm.eg.db, cluster_genes, 'ENTREZID', 'SYMBOL')
  cluster_go <- enrichGO(gene = cluster_genes,
                         'org.Mm.eg.db',
                         keyType = "ENTREZID",
                         ont = "BP",
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH")
  p1 <- as.data.frame(cluster_go@result) %>% filter(pvalue < 0.05)
  p1$log_p <- -log10(p1$pvalue)
  color_threshold <- max(p1$log_p)+1
  dotplot_GO <- ggplot(p1[c(1,2,3,4,5,6,7,8,9,10,11,12),],  # you can replace the numbers to the row number of pathway of your interest
                       aes(x = Count, y = reorder(Description, log_p))) + 
    geom_point(aes(size = Count, color = log_p)) +
    theme_bw(base_size = 14) +
    theme +
    scale_colour_gradient2(limits=c(0, color_threshold), low="white", mid="white", high = "purple") +
    ylab(NULL) +
    xlab('Count of Genes') +
    labs(title = str_wrap(paste0('GO - Biological Pathway Cluster',i), 20))
  dotplot_GO
  
  ggsave(paste0('plots/',file_prefix,'_cluster',i,'_GO_pathway','.pdf'),
         dotplot_GO,
         dpi = 300)
}



