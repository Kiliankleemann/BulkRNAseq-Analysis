#Volcano
volcano_data <- read.xlsx(paste0(DRIRECTORY/NAME_OF_FILE))
volcano_data$log_pval <- -log10(volcano_data$p_val)

UP_data <- volcano_data %>% filter(p_val < 0.05 & avg_log2FC > 0)
DOWN_data <- volcano_data %>% filter(p_val < 0.05 & avg_log2FC < 0)

UP_DEGs <- UP_data %>% pull('gene')
UP_DEGs <- mapIds(org.Mm.eg.db, UP_DEGs, 'ENTREZID', 'SYMBOL')

DOWN_DEGs <- DOWN_data %>% pull('gene')
DOWN_DEGs <- mapIds(org.Mm.eg.db, DOWN_DEGs, 'ENTREZID', 'SYMBOL')

#Pathways
go <- enrichGO(gene = DOWN_DEGs,
               'org.Mm.eg.db',
               keyType = "ENTREZID",
               ont = "BP",
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH")
p1 <- as.data.frame(go@result) %>% filter(pvalue < 0.06)
p1$log_p <- -log10(p1$pvalue)

pathway_data <- AnnotationDbi::select(org.Mm.eg.db, keytype="GOALL", keys='GO:0043161', columns = 'SYMBOL')
pathway_genes <- pathway_data %>% pull (SYMBOL)
Proteasome_data <- DOWN_data %>% filter(gene %in% c(pathway_genes))
Antigen_data <- UP_data %>% filter(gene %in% c(pathway_genes))

#Theme
theme <- theme(aspect.ratio = 1, 
               panel.background = element_blank(),
               panel.border=element_rect(fill=NA),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               strip.background=element_blank(),
               axis.title = element_text(size = 20, colour = "black"),
               axis.text.x=element_text(colour="black",size =20),
               axis.text.y=element_text(colour="black",size =20),
               axis.ticks=element_line(colour="black"),
               plot.margin=unit(c(1,1,1,1),"line"),
               plot.title = element_text(colour="black",size =20))


scatter <- ggplot(volcano_data, aes(x= avg_log2FC, y= log_pval, label = gene)) +
  #geom_smooth(data = MGnD_data, method = "lm", se=TRUE, color="red", formula = y ~ x) +
  #stat_cor(data = MGnD_data, color = 'red', label.x.npc = "left", label.y.npc = "top") +
  #geom_smooth(data = M0_data, method = "lm", se=TRUE, color="dodgerblue", formula = y ~ x) +
  #stat_cor(data = M0_data, color = 'dodgerblue', label.x.npc = "center", label.y.npc = "bottom") +
  geom_point(size=0.5, color = 'grey90') +
  geom_point(data = UP_data, color = 'red', size = 1)+
  geom_point(data = DOWN_data, color = 'dodgerblue3', size = 1) +
  #geom_point(data = Antigen_data, color = 'chartreuse4', size = 1)+
  theme +
  xlim(-5,5)+
  ylim(0,20)+
  geom_vline(xintercept = 0.15,colour="grey", linetype = "longdash") +
  geom_vline(xintercept = -0.15,colour="grey", linetype = "longdash") +
  geom_hline(yintercept = 1.3, colour="grey", linetype = "longdash") +
  labs(title = str_wrap("B1 (2 Months) APOE4 vs APOE3",60),
       x = "Log2FC",
       y = "-log(pvalue)") +
  #geom_label(data = data_selected, color = 'orange', size = 5, label.padding = unit(0.1,"lines")) +
  #geom_label(data = Antigen_data, color = 'chartreuse4', size = 2, label.padding = unit(0.1, "lines")) 
  geom_label_repel(data = Antigen_data,
                   color = 'orangered', size = 3,
                   min.segment.length = 0,fontface = 'italic',
                   label.padding = unit(0.2,"lines"),
                   nudge_y = 1,
                   nudge_x = 2.5,
                   segment.angle = 20,
                   max.overlaps = 15) +
  annotate("text", x = 2.5, y = 19, label= str_wrap("Proteasome processing",10), size = 6) + 
  geom_label_repel(data = Proteasome_data,
                   color = 'dodgerblue1',size = 3, fontface = 'italic',
                   min.segment.length =0,
                   nudge_y = 2,
                   nudge_x = -2.5,
                   label.padding = unit(0.2,"lines"),
                   max.overlaps = 20) +
  annotate("text", x = -2.5, y = 19, label= str_wrap("Antigen presentation",10), size = 6)  
scatter 


ggsave(
  paste0('plots/Volcano_B1_2Months_APOE4vsAPOE3', '.pdf'),
  plot = scatter,
  device = NULL,
  path = NULL,
  scale = 1,
  width = NA,
  height = NA,
  units = c("mm"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)
