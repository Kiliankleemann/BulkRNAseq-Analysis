# Define a normalization function to calculate Z scores.
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# Import in data and subset it for how many genes you want to see.
data <- read.xlsx(paste0('results/', file_prefix, '/DEGene_counts_padj05.xlsx'), rowNames = T)
file_prefix

#Data subsetting and gene filtering (if necessary)
#unwanted_genes = paste(c('^Gm', '^mt-', '^Vmn', '^Rpl', '^Rps', '^Olfr','Rik$'), collapse = '|')

data_final <- data %>%
  rownames_to_column('gene') %>%
  #filter(!str_detect(gene, c(unwanted_genes))) %>%
  #anti_join(unwanted_clusters, by = 'gene') %>%
  column_to_rownames('gene')
data_subset_50 <- data_final %>% head(n = 50)
data_subset_200 <- data_final %>% head(n = 200)
data_subset_400 <- data_final %>% head(n = 400)


#Grouped data final
data_grouped <- data.frame(
  control = rowMeans(select(data_final, contains(("ctrl")))),
  Sirt1_KO = rowMeans(select(data_final, contains(("Sirt"))))
  )


# Normalize your data according to Z score using cal_z_score function.
data_norm <- t(apply(data_final, 1, cal_z_score))
data_subset_200 <- t(apply(data_subset_200, 1, cal_z_score))
data_subset_400 <- t(apply(data_subset_400, 1, cal_z_score))
data_grouped_norm <- t(apply(data_grouped, 1, cal_z_score))

#Export directory
dir.create(paste0('plots/heatmaps/'))
dir.create(paste0('plots/heatmaps/',file_prefix))

#Setting up Heatmap parameters
gaps = c(4)
clusters = 2
color_set <- c("navy", "white", "firebrick3")
all_degs_export_height <- nrow(data_final)

## ------------- DENDROGRAM AND CLUSTER MANIPULATION (ADVANCED) -------------##
# Create heatmap
dev.off()
pdf(file = paste0('plots/heatmaps/',file_prefix,'/','All_DEGs_pval05.pdf'), pointsize = 10)
phm_full <- pheatmap(data_norm,
                     color = colorRampPalette(color_set)(41),
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


pdf(file = paste0('plots/heatmaps/', file_prefix,'/', 'All_DEGs_labelled_pval05.pdf'), pointsize = 10, height = all_degs_export_height/4)
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

pdf(file = paste0('plots/heatmaps/', file_prefix,'/', 'Top_50.pdf'), pointsize = 10, height = 10)
phm_50 <- pheatmap(data_subset_50,
                    color = colorRampPalette(color_set)(41),
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

pdf(file = paste0('plots/heatmaps/', file_prefix,'/', 'Top_200.pdf'), pointsize = 10, height = 35)
phm_200 <- pheatmap(data_subset_200,
                    color = colorRampPalette(color_set)(41),
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

pdf(file = paste0('plots/heatmaps/', file_prefix,'/', 'Top_400.pdf'), pointsize = 10, height = 60)
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


### GOURPED
pdf(file = paste0('plots/heatmaps/', file_prefix,'/', 'All_DEGs_grouped.pdf'), pointsize = 10)
phm_full <- pheatmap(data_grouped_norm,
                     color = colorRampPalette(c(color_set))(41),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = F,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     #cellwidth = 20,
                     #cellheight = 9,
                     scale = 'row')
dev.off()


#Grouped labelled
pdf(file = paste0('plots/heatmaps/', file_prefix, '/','All_DEGs_labelled_grouped.pdf'), pointsize = 10, height = all_degs_export_height/4)
phm_full <- pheatmap(data_grouped_norm,
                     color = colorRampPalette(c(color_set))(41),
                     kmeans_k = NA,
                     cluster_rows = T,
                     cutree_row = clusters,
                     cluster_cols = F,
                     #cutree_cols = 3,
                     gaps_col = c(),             
                     legend = TRUE,
                     show_rownames = T,
                     border_color = 'NA',
                     fontsize = 8,
                     treeheight_col = 0,
                     treeheight_row = 0,
                     cellwidth = 20,
                     cellheight = 9,
                     scale = 'row')
dev.off()





####################################### CLUSTER EXTRACTION ####################################33
# Build cluster tables from dendrogram based on k
gene_clusters <- data.frame(sort(cutree(phm_full$tree_row, k=4))) %>%
  rownames_to_column(var="gene") 

temp = gene_clusters$sort.cutree.phm_full.tree_row..k...4..

gene_clusters <- gene_clusters %>% 
  mutate(cluster = temp) %>%
  select(gene, cluster)
gene_clusters

write.xlsx(gene_clusters, file=paste0('results/', file_prefix, 'clustergenes.xlsx'))


# Subset and organize clusters (iterative process)
gene_clusters_clean <- gene_clusters %>%
  arrange(match(cluster, c(1,2,3,4,5,6)))


dev.off()
                     
###################################### HEATMAP CUSTOMIZATION ##################################
# Subset and organize clusters (iterative process)
gene_clusters_clean <- gene_clusters %>%
  arrange(match(cluster, c(1,2,3,4,5,6)))

# Reorder transcript count data to match custom order
#data_ordered <- data_final[gene_clusters_clean$gene, ]    
data_ordered <- data_final[gene_clusters_clean$gene, ]
rownames(gene_clusters) <- rownames(data_ordered)
gene_clusters

# Build annotation dataframes metadata
gene_clusters_clean <- gene_clusters_clean %>% mutate(primary_effect = ifelse(cluster == 1, 'cluster_1',
                        ifelse(cluster == 2, 'cluster_2',
                        ifelse(cluster == 3, 'cluster_3',
                        ifelse(cluster == 4, 'cluster_4',          
                        ifelse(cluster == 5, 'cluster_5', 
                         'cluster_6' ))))))
                        
# Subset out annotation data for plotting
cluster_anno <- gene_clusters_clean %>% select(primary_effect)
rownames(cluster_anno) <- rownames(data_ordered)


# creat colours for each group
newCols <- colorRampPalette(grDevices::rainbow(length(unique(cluster_anno$primary_effect))))
annoCol <- newCols(length(unique(cluster_anno$primary_effect)))
names(annoCol) <- unique(cluster_anno$primary_effect)
annoCol <- list(category = annoCol)

# make the heatmap
pheatmap(zz, scale = "row", cluster_cols = F, cluster_rows = F, annotation_row = anno, annotation_colors = annoCol, cellheight = 10, cellwidth = 10, file = "expr_group.png")

f <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)
paired <- as.list((cols <- f("Paired")))

pdf(file = paste0('plots/heatmap_cluster_annotation_', file_prefix,'.pdf'), pointsize = 10, height = 10)
phm_clusters <- pheatmap(data_ordered,
                         color = colorRampPalette(c("navy", "white", "firebrick3"))(41),
                         breaks = seq(-2, 2, by = 0.1),
                         kmeans_k = NA,
                         cluster_rows = T,
                         cutree_row = 6,
                         cluster_cols = F,
                         #cutree_cols = 4,
                         gaps_col = c(4,8,12),
                         annotation_row = cluster_anno,
                         annotation_colors = annoCol,
                         legend = TRUE,
                         show_rownames = F,
                         border_color = 'NA',
                         fontsize = 8,
                         treeheight_col = 0,
                         treeheight_row = 10,
                         #cellwidth = 25,
                         scale = 'row')
dev.off()

# Write cluster metadata to file
write.xlsx(gene_clusters_clean, file = paste0('results/', file_prefix, 'DE_clusters.xlsx'), overwrite = T)
