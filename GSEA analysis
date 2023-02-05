######--------GSEA analysis--------########

library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(clusterProfiler)
library(limma)
library(openxlsx)
library(dplyr)
library(tidyverse)
library(data.table)
library(stringr)
library(fgsea)
library(enrichplot)


file_prefix <- 'Wald_test/Control_data_APOE34_vs_APOE33_'
file_export <- 'Control_data_APOE34_vs_APOE33_'
pathway_data <- read.xlsx(paste0('results/',file_prefix, 'DEGene_statistics_pval05.xlsx'))
DEgenes <- pathway_data %>% filter(pvalue < 0.05)  %>% pull('gene')
DEgenes <- mapIds(org.Hs.eg.db, DEgenes, 'ENTREZID', 'SYMBOL')
DEG_Log2FC <-  pathway_data %>% filter(pvalue < 0.05) %>%  pull('log2FoldChange')

names(DEG_Log2FC) = DEgenes
DEG_Log2FC <- sort(DEG_Log2FC, decreasing = TRUE)
DEG_Log2FC

#KEGG enrichment analsis
kegg_result <- gseKEGG(geneList = DEG_Log2FC,
                       organism     = 'hsa',
                       pAdjustMethod = "BH",
                       # minGSSize = 15, 
                       # maxGSSize = 500,
                       nPermSimple = 10000,
                       pvalueCutoff = 1)

head(kegg_result)
kegg_result <- setReadable(kegg_result, 'org.Hs.eg.db', 'ENTREZID')
kegg_result_dataframe <- as.data.frame(kegg_result)
write.xlsx(kegg_result_dataframe, paste0('results/Pathway_analysis/KEGG_enrichment/',file_prefix, '.xlsx'))

#BAR plot
gsea_KEGG <- gseaplot2(kegg_result, geneSetID = c(1,5,7,10), pvalue_table = T,ES_geom = "line", color = scale_fill_brewer(palette = "Paired"))
gsea_KEGG
ggsave(paste0('plots/pathways/KEGG_pathways/', file_export,'.pdf'),
       gsea_KEGG,
       dpi = 300)


# fgseaResTidy_curated_pval05 <- kegg_result_dataframe %>% filter(pathway %in% c('Neutrophil extracellular trap formation',
#                                                                                ))
# fgseaResTidy_curated_pval05$log_pval <- -log10(fgseaResTidy_curated_pval05$pval)
# barplot_curated <- ggplot(fgseaResTidy_curated_pval05, aes(reorder(pathway, NES), NES)) +
#   geom_col(aes(fill= NES), color = 'black') +
#   coord_flip() +
#   theme_gsea +
#   scale_fill_gradient2(low = 'blue', mid = "white",
#                        high = 'red', midpoint = 0, space = "rgb",
#                        na.value = "grey50", guide = "colourbar") +
#   labs(x="Pathway", y="Normalized Enrichment Score",
#        title="Curated pathways NES from GSEA") + 
#   geom_vline(xintercept = 0, color ='black')



############ --------- GO enrichment analsis CELLULAR COMPONENT ----------############
go_result <- gseGO(geneList = DEG_Log2FC,
                   OrgDb    = org.Hs.eg.db,
                   ont      = "CC",
                   # minGSSize = 15, 
                   # maxGSSize = 500,
                   nPermSimple = 10000,
                   pvalueCutoff = 1)

head(go_result)
go_result <- setReadable(go_result, 'org.Hs.eg.db', 'ENTREZID')
go_result_dataframe <- as.data.frame(go_result)
write.xlsx(kegg_result_dataframe, paste0('results/Pathway_analysis/GO_enrichment/',file_prefix, 'cellular_componenet.xlsx'))

#Gsea plot
dotplot_GO <- gseaplot2(go_result, geneSetID = c(4), pvalue_table = T,ES_geom = "line", color = scale_fill_brewer(palette = "Paired"))
dotplot_GO
ggsave(paste0('plots/pathways/GO_pathways/', file_prefix,'cellular_componenet.pdf'),
       dotplot_GO,
       dpi = 300)

#GO enrichment analsis BIOLOGICAL PROCESS
go_result <- gseGO(geneList = DEG_Log2FC,
                   OrgDb    = org.Hs.eg.db,
                   ont      = "BP",
                   # minGSSize = 15, 
                   # maxGSSize = 500,
                   nPermSimple = 10000,
                   pvalueCutoff = 1)

go_result <- setReadable(go_result, 'org.Hs.eg.db', 'ENTREZID')
go_result_dataframe <- as.data.frame(go_result)
write.xlsx(kegg_result_dataframe, paste0('results/Pathway_analysis/GO_enrichment/',file_prefix, 'biological_process.xlsx'))

#Gsea plot
#select_result <- go_result %>% filter(ID %in% c('GO:0002366','GO:0001818','GO:0051606','GO:0007186'))
dotplot_GO <- gseaplot2(go_result, geneSetID = c(1,2,3,4), pvalue_table = T,ES_geom = "line", color = scale_fill_brewer(palette = "Paired"))
dotplot_GO
ggsave(paste0('plots/pathways/GO_pathways/', file_export,'biological_process.pdf'),
       dotplot_GO,
       dpi = 300)

#GO enrichment analsis BIOLOGICAL PROCESS
go_result <- gseGO(geneList = DEG_Log2FC,
                   OrgDb    = org.Hs.eg.db,
                   ont      = "MF",
                   # minGSSize = 15, 
                   # maxGSSize = 500,
                   nPermSimple = 10000,
                   pvalueCutoff = 1)

head(go_result)
go_result <- setReadable(go_result, 'org.Hs.eg.db', 'ENTREZID')
go_result_dataframe <- as.data.frame(go_result)
write.xlsx(kegg_result_dataframe, paste0('results/Pathway_analysis/GO_enrichment/',file_prefix, 'molecular_function.xlsx'))

#Gsea plot
dotplot_GO <- gseaplot2(go_result, geneSetID = c(1,2,3), pvalue_table = T,ES_geom = "line", color = scale_fill_brewer(palette = "Paired"))
dotplot_GO
ggsave(paste0('plots/GO_pathways/', file_prefix,'cellular_componenet.pdf'),
       dotplot_GO,
       dpi = 300)




#fgsea analysis
#Fgsea plot theme
theme_gsea <- theme(
  text = element_text(size = pointSize, colour = "black"),
  rect = element_blank(),
  line = element_line(size = lineWidth, colour = "black"),
  plot.title  = element_text(color="black", size = pointSize),
  axis.title  = element_text(size = pointSize * 0.8, colour = "black"),
  axis.text.x  = element_text(size = pointSize * 0.8, colour = "black"),
  axis.text.y  = element_text(size = pointSize * 0.8, colour = "black"),
  axis.ticks.x = element_line(size = 0.5, colour = "black"),
  axis.ticks.y = element_line(size = 0.5, colour = "black"),
  axis.line = element_line(size = lineWidth, colour = "black"),
  axis.line.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.title = element_text(size = pointSize * 0.7, colour = "black"),
  legend.text = element_text(size = pointSize * 0.4, colour = "black"),
  legend.key.height = unit(0.2, "cm"),
  legend.key.width = unit(0.2, "cm"))

#####CURATED pathways
pathways.all.curated <- gmtPathways("GSEA_database/c2.all.v2022.1.Hs.entrez.gmt") 
fgseaRes <- fgsea(pathways=pathways.all.curated, stats=DEG_Log2FC)

fgseaResTidy_curated <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  filter(pval<0.1)

write.xlsx(fgseaResTidy_curated, paste0('results/Pathway_analysis/FGSEA_enrichment/',file_export, 'Curated_pathways.xlsx'))

#GSEA plot
fgseaResTidy_curated_pval05 <- fgseaResTidy_curated %>% filter(pval < 0.05)
fgseaResTidy_curated_pval05 <- fgseaResTidy_curated_pval05 %>% filter(pathway %in% c('REACTOME_NEUTROPHIL_DEGRANULATION','LIAN_NEUTROPHIL_GRANULE_CONSTITUENTS','MARTINELLI_IMMATURE_NEUTROPHIL_UP',
                                                                                'MAHAJAN_RESPONSE_TO_IL1A_UP', 'REACTOME_FCGR3A_MEDIATED_IL10_SYNTHESIS', 'KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY',
                                                                                'BIOCARTA_IL17_PATHWAY','WP_IL17_SIGNALING_PATHWAY','MAHAJAN_RESPONSE_TO_IL1A_UP','REACTOME_COLLAGEN_DEGRADATION',
                                                                                'COULOUARN_TEMPORAL_TGFB1_SIGNATURE_UP','WP_ALLOGRAFT_REJECTION','REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX',
                                                                                'VERRECCHIA_RESPONSE_TO_TGFB1_C2', 'REACTOME_REGULATION_OF_IFNA_IFNB_SIGNALING','MARZEC_IL2_SIGNALING_UP',
                                                                                'GILMORE_CORE_NFKB_PATHWAY','KEGG_LYSOSOME','BIOCARTA_NEUTROPHIL_PATHWAY','REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS',
                                                                                'WP_DEGRADATION_PATHWAY_OF_SPHINGOLIPIDS_INCLUDING_DISEASES','LIU_IL13_PRIMING_MODEL',
                                                                                'WP_IL4_SIGNALING_PATHWAY', 'REACTOME_SIGNALING_BY_TGFB_FAMILY_MEMBERS','REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX'))
fgseaResTidy_curated_pval05$log_pval <- -log10(fgseaResTidy_curated_pval05$pval)
barplot_curated <- ggplot(fgseaResTidy_curated_pval05, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES), color = 'black') +
  coord_flip() +
  theme_gsea +
  scale_fill_gradient2(low = 'blue', mid = "white",
                       high = 'red', midpoint = 0, space = "rgb",
                       na.value = "grey50", guide = "colourbar") +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Curated pathways NES from GSEA") + 
  geom_vline(xintercept = 0, color ='black')
 
barplot_curated

ggsave(paste0('plots/pathways/FGSEA_pathways/', file_export,'curated_pathways.pdf'),
       barplot_curated,
       dpi = 300)



#Immunological pathways
pathways.immuno <- gmtPathways("GSEA_database/immunologic_sets_c7.all.v2022.1.Hs.entrez.gmt") 
fgseaRes <- fgsea(pathways=pathways.immuno, stats=DEG_Log2FC)

fgseaResTid_immuno <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  filter(pval<0.1)

write.xlsx(fgseaResTid_immuno, paste0('results/Pathway_analysis/FGSEA_enrichment/',file_export, 'Immunological_pathways.xlsx'))

#Gsea plot
fgseaResTid_immuno$log_pval <- -log10(fgseaResTid_immuno$pval)
barplot_immuno <- ggplot(fgseaResTid_immuno, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES), color = 'black') +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  geom_vline(xintercept = 0, color ='black')+
  theme_gsea
barplot_immuno


#Hallmark Pathways
pathways.hallmark <- gmtPathways("GSEA_database/hallmarks.all.v2022.1.Hs.entrez.gmt") 
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=DEG_Log2FC)

fgseaResTidy_hallmark <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  filter(pval<0.05)

write.xlsx(fgseaResTidy, paste0('results/Pathway_analysis/FGSEA_enrichment/',file_export, 'Hallmark_pathways.xlsx'))

#Gsea plot
fgseaResTidy_hallmark$log_pval <- -log10(fgseaResTidy_hallmark$pval)
barplot_hallmarks <- ggplot(fgseaResTidy_hallmark, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES), color = 'black') +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  geom_vline(xintercept = 0, color ='black')+
  theme_gsea+
  scale_fill_gradient2(low = 'blue', mid = "white",
                        high = 'red', midpoint = 0, space = "rgb",
                        na.value = "grey50", guide = "colourbar")
barplot_hallmarks

ggsave(paste0('plots/pathways/FGSEA_pathways/', file_export,'Hallmark_pathways.pdf'),
       barplot_hallmarks,
       dpi = 300)









