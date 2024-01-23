#Barplots of results from DESEQ - Adding DESeq2 pvalue (for one-way comparison)
file_prefix <- ''
check_counts <- read.xlsx('results/', file_prefix, '/DS_counts_check.xlsx')
sig_res <- read.xlsx('results/',file_prefix,  '/DEGene_statistics_pval05.xlsx')

goi <- c('')
#goi <- sig_res %>% head(20) %>% pull('gene')
goi_data <-  check_counts %>% filter(gene %in% goi)
goi_data2 <- rename(goi_data, 
                    'WT_' = contains('WT'),
                    'KI_' = contains('KI')) 

lineWidth = 1
pointSize = 40
for (i in goi) { 
  goi_data3 <- goi_data2 %>% filter(gene == i) %>% pivot_longer(
    cols = -(1:1),
    values_to = c("DS_Counts"),
    names_to = c("condition", "replicate"),
    names_sep  = "_")
  
  pvalue_data <- round(sig_res %>% filter(gene == i) %>% pull(padj),digits = 10)
  
  goi_data3$condition <- factor(goi_data3$condition, levels=unique(goi_data3$condition))
  
  comparisons <- compare_means(
    data = goi_data3,
    formula = DS_Counts ~ condition,
    method = "t.test",
    p.adjust.method = "BH") #%>% filter(p<0.05)
  
  comparison_list_sign <- comparisons %>% mutate(comparison_list = map2(group1, group2,c)) %>% pull(comparison_list)
  
  barplot <- ggplot(goi_data3, aes(x = condition, y = DS_Counts, fill = condition)) + 
    geom_bar(position = 'dodge', stat = 'summary', fun = mean, width = 0.7, colour="black") +
    geom_errorbar(stat = 'summary', position = 'dodge', width = 0.4) +
    geom_point(aes(x = condition),size=4,shape=21, color ='white', position = position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width=0.9)) +
    expand_limits(x = 0, y = 0) +
    theme(panel.spacing = unit(1, "lines")) +
    labs(title = paste0(i),
         subtitle = paste0("Padj = ", pvalue_data),
         x = NULL, 
         y = c("Normalized Counts ±SE")) +
    ggtitle(paste0(i))+
    theme(
      text = element_text(size = pointSize, colour = "black"),
      rect = element_blank(),
      line = element_line(size = lineWidth, colour = "black"),
      plot.title  = element_text(color="black", size=pointSize),
      plot.subtitle= element_text(size=pointSize *0.8, hjust= 1.5, face="italic", color="black"),
      axis.title  = element_text(size = pointSize, colour = "black"),
      axis.text.x  = element_text(size = pointSize , colour = "black", angle = 45, hjust = 1),
      axis.text.y  = element_text(size = pointSize , colour = "black"),
      axis.ticks.x = element_line(size = lineWidth, colour = "black"),
      axis.ticks.y = element_line(size = lineWidth, colour = "black"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = pointSize , colour = "black"),
      # legend.key.height = unit(0.1, "cm"),
      # legend.key.width = unit(0.2, "cm"),
      axis.line = element_line(size = lineWidth, colour = "black"),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
    #stat_compare_means(comparisons = comparison_list_sign , size = 3) +
    #scale_fill_brewer(palette = "Paired") + 
    scale_fill_manual(values = c('black',"red3")) +
    scale_y_continuous(expand = c(0, 0, 0.1, 0))
  
  pdf(file = paste0('plots/barplots/',file_prefix, i,'.pdf'), pointsize = 10, width = 5, height= 10)
  print(barplot)
  dev.off()
}
