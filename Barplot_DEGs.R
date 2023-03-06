#Barplots
BiocManager::install('ggsignif')
BiocManager::install('ggpubr')
library(ggsignif)
library(ggpubr)
library(cowplot)
library(RColorBrewer)

comparisons <- goi_data3 %>%
  group_by(condition) %>%
  distinct(condition) %>% pull(condition)

goi <- c('Esr1','Cx3cr1','Stat1','Clec7a','Axl')
goi_data <-  sig_export %>% filter(gene %in% c(goi))
goi_data2 <- rename(goi_data, 
                    'NonPhag-WT-IgG_' = contains('Non_phag_WT_IgG'),
                    'NonPhag-WT-Blocker_' = contains('Non_phag_WT_Blocker'),
                    'NonPhag-miR155cKO-IgG_' = contains('Non_phag_miR155cKO_IgG'),
                    'NonPhag-miR155cKO-Blocker_' = contains('Non_phag_miR155cKO_Blocker'),
                    'Phagocytic-WT-IgG_' = contains('Phagocytic_WT_IgG'),
                    'Phagocytic-WT-Blocker_' = contains('Phagocytic_WT_Blocker'),
                    'Phagocytic-miR155cKO-IgG_' = contains('Phagocytic_miR155cKO_IgG'),
                    'Phagocytic-miR155cKO-Blocker_' = contains('Phagocytic_miR155cKO_Blocker')
                   ) 

for (i in goi) { 
  goi_data3 <- goi_data2 %>% filter(gene == i) %>% pivot_longer(
  cols = -(1:1),
  values_to = c("DS_Counts"),
  names_to = c("condition", "replicate"),
  names_sep  = "_")
  
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
  geom_point(aes(x = condition), position = position_jitterdodge(jitter.width = 0.3, jitter.height=0.4, dodge.width=0.9)) +
  expand_limits(x = 0, y = 0) +
  theme(panel.spacing = unit(1, "lines")) +
  labs(x = NULL, y = c("Normalized Counts ±SE")) +
  theme(
    text = element_text(size = pointSize, colour = "black"),
    rect = element_blank(),
    line = element_line(size = lineWidth, colour = "black"),
    plot.title  = element_text(color="black", size=1, face="bold.italic"),
    axis.title  = element_text(size = pointSize * 0.8, colour = "black"),
    axis.text.x  = element_text(size = pointSize * 0.6, colour = "black", angle = 45, hjust = 1),
    axis.text.y  = element_text(size = pointSize * 0.6, colour = "black"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = pointSize * 0.6, colour = "black"),
    legend.key.height = unit(0.1, "cm"),
    legend.key.width = unit(0.2, "cm"),
    axis.line = element_line(size = lineWidth, colour = "black"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
  stat_compare_means(comparisons = comparison_list_sign, method = 't.test', label = "p.signif", ) +
  scale_fill_brewer(palette = "Paired") + scale_y_continuous(expand = c(0, 0, .05, 0))
  
  pdf(file = paste0('plots/barplot_', i,'_',file_prefix,'.pdf'), pointsize = 10)
  print(barplot)
  dev.off()
}









