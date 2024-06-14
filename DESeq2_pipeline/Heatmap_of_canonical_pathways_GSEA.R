library(ComplecHeatmap)
library(circlise)


#Download canonical gene set from GSEA (m2.cp.v2022.1.Mm.entrez.gmt)

#Canonical gene sets
pathways.canoical <- gmtPathways("gsea_pathways/m2.cp.v2022.1.Mm.entrez.gmt") 
rm(canonical_data_total)
canonical_data_total <- data.frame()
for (i in SubCellTypes_no_macrophage){
  pathway_analysis <- pathway_data %>% filter(cluster == paste0(i))
  DEgenes <- pathway_analysis %>% filter(p_val < 0.05) %>% pull('gene')
  DEgenes <- mapIds(org.Mm.eg.db, DEgenes, 'ENTREZID', 'SYMBOL')
  DEG_Log2FC <-  pathway_analysis %>% filter(p_val < 0.05) %>%  pull('avg_log2FC')
  names(DEG_Log2FC) = DEgenes
  DEG_Log2FC <- sort(DEG_Log2FC, decreasing = TRUE)
  fgseaRes <- fgsea(pathways=pathways.canoical, stats=DEG_Log2FC, scoreType ='pos')
  fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  filter(pval<0.05)
  write.xlsx(fgseaResTidy, paste0('results/Hallmark_pathways_',i, 'Hallmark_pathways.xlsx'))
  fgseaResTidy$log_p <- -log10(fgseaResTidy$pval)
  fgseaResTidy$group <- paste0(i)
  canonical_data_total <- rbind(canonical_data_total, fgseaResTidy)
}


## Heatmap
canonical_selected_pathways <- read.xlsx(paste0('results/canonical_data_total.xlsx'), sheet = 2)%>% pull(pathway)
canonical_data_selected <- canonical_data_total %>% filter(pathway %in% canonical_selected_pathways)

canonical_data_selected_pval <- canonical_data_selected %>% select(c('pathway','group', 'log_p'))
canonical_data_selected_NES <- canonical_data_selected %>% select(c('pathway','group', 'NES'))


canonical_data_selected_pval_mt <- canonical_data_selected_pval %>% 
  as_tibble() %>% pivot_wider(id_cols = pathway,
                              names_from = group,
                              values_from = log_p,
                              values_fill = 0)%>%
  column_to_rownames('pathway')

canonical_data_selected_NES_mt <- canonical_data_selected_NES %>% 
  as_tibble() %>% pivot_wider(id_cols = pathway,
                              names_from = group,
                              values_from = NES,
                              values_fill = 0)%>%
  column_to_rownames('pathway')

#Change Column order
mt_colum_order <- c( "M0","Ribosome-Intermediate","HSP-Intermediate","Cytokine-Intermediate","Interferon","Cycling-G2M","Cycling-S","MGnD","MGnD-Antigen")
canonical_data_selected_NES_mt <- canonical_data_selected_NES_mt[,mt_colum_order]
canonical_data_selected_pval_mt <- canonical_data_selected_pval_mt[,mt_colum_order]

mt_row_order <- read.xlsx(paste0('results/canonical_data_total.xlsx'), sheet = 3)%>% pull(ordered_pathways)
canonical_data_selected_NES_mt <- canonical_data_selected_NES_mt[mt_row_order,]
canonical_data_selected_pval_mt <- canonical_data_selected_pval_mt[mt_row_order,]

#Edit rownames
rownames(canonical_data_selected_pval_mt) <- sub('REACTOME_','',rownames(canonical_data_selected_pval_mt))
rownames(canonical_data_selected_pval_mt) <- sub('BIOCARTA_','',rownames(canonical_data_selected_pval_mt))
rownames(canonical_data_selected_pval_mt) <- sub('WP_','',rownames(canonical_data_selected_pval_mt))
rownames(canonical_data_selected_pval_mt) <- gsub('_',' ',rownames(canonical_data_selected_pval_mt))

rownames(canonical_data_selected_NES_mt) <- sub('REACTOME_','',rownames(canonical_data_selected_NES_mt))
rownames(canonical_data_selected_NES_mt) <- sub('BIOCARTA_','',rownames(canonical_data_selected_NES_mt))
rownames(canonical_data_selected_NES_mt) <- sub('WP_','',rownames(canonical_data_selected_NES_mt))
rownames(canonical_data_selected_NES_mt) <- gsub('_',' ',rownames(canonical_data_selected_NES_mt))

canonical_data_selected_pval_mt <- as.matrix(canonical_data_selected_pval_mt)
canonical_data_selected_NES_mt <- as.matrix(canonical_data_selected_NES_mt)


#Set up legend
col_fun = colorRamp2(c(1.3, 2, 2.6), c("white", "orange", "red"))
col_fun(seq(1.3, 2.6))
lgd = Legend(col_fun = col_fun, title = "NES")

htmp_pval <- Heatmap(canonical_data_selected_NES_mt,
        col = col_fun,
        rect_gp = gpar(col = "black", lwd = 0.1),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        width = ncol(canonical_data_selected_pval_mt)*unit(5, "mm"), 
        height = nrow(canonical_data_selected_pval_mt)*unit(5, "mm"),
        show_heatmap_legend = F,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(canonical_data_selected_pval_mt[i, j] > 2) {
            grid.text("***", x, y)
          } else if(canonical_data_selected_pval_mt[i, j] > 1.6) {
            grid.text("**", x, y)
          } else if(canonical_data_selected_pval_mt[i, j] > 1.3 ) {
            grid.text("*", x, y)
          }else if(canonical_data_selected_pval_mt[i, j] < 1.3) {
            grid.text("", x, y)
          }
        }) 

pdf(file = paste0('plots/heatmap_canonical_pathways_', 'SubCellType', '.pdf'), width = 10, height = 12)
htmp_pval
draw(lgd, x = unit(0.5, "cm"), just = c("left"))
dev.off()
