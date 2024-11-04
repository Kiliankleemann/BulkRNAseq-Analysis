Combined significant data
#MGnD paradigm
IL17RA_KO_results <- read.xlsx('/Users/kiliankleemann/Dropbox/Neutrophils Projects/SK-5R3J_GSK_Neta_Neutro_Jul_2023/results/MG_AD_IL17_KO_data_DEGene_statistics_pval05.xlsx')
IL17RA_KO_DEGenes <- IL17RA_KO_results %>% pull(gene)

#Gene list 1
data_result1 <- read.xlsx('results/data_result1.xlsx') %>% pull(gene)

#Gene list 2
data_result2 <- read.xlsx('results/data_result2.xlsx') %>% pull(gene)

#All data APOE comparisons 
overlap_genes <-  Reduce(intersect,list(data_result2, data_result1))
overlap_genes

#Calculate overlap significance - Chi squared test
nrow(data_result1)
nrow(overlap_genes)
x = 60      #Overlap 
m = 785     #data_result1
n = 11878   #Total number of instances tested i.e. total RNA editing sites
k = 728     #data_result2
hyper_overlap_stat <- phyper(x, m, n, k, lower.tail = FALSE)
hyper_overlap_stat


# Make venn diagram for overlap
# Chart
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
venn <- venn.diagram(
  x = list(data_result1,data_result2),
  category.names = c('data_result1_name', 'data_result2_name'),
  #filename = paste0('plots/overlap_comparison/',gene_set_comparison_name,'.pdf'),
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
  col=c("#2A2AE1","#FF0090" ), #'orchid3'), 
  fill = c(alpha("#2A2AE1" ,0.5),alpha("#FF0090",0.5)), #alpha('orchid3',0.5)),#
  # Numbers
  cex = 5,
  #fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 0.5,
  cat.default.pos = "outer",
  cat.pos = c(1,1,1),
  cat.dist = c(0.3,0.3,0.5),
  #cat.fontface = "bold",
  cat.col = c("#2A2AE1","#FF0090"), #"orchid3"),
  cat.fontfamily = "sans")
grid.draw(venn)


pdf(file=paste0('plots/overlap_comparison/','Gene_set_comparison','.pdf'))
grid.draw(venn)
dev.off()